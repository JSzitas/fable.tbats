// Tiny QR solver, header only library
//
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// Copyright (C) 2023- Juraj Szitas
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef TINYQR_H_
#define TINYQR_H_

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>
// inlining
#define INLINE_THIS
#if defined(__clang__) || defined(__GNUC__)
#undef INLINE_THIS
#define INLINE_THIS __attribute__((always_inline))
#elif defined(_MSC_VER)
#undef INLINE_THIS
#define INLINE_THIS __forceinline
#endif
// restrict
#define RESTRICT_THIS
#if defined(__clang__) || defined(__GNUC__)
#undef RESTRICT_THIS
#define RESTRICT_THIS __restrict__
#elif defined(_MSC_VER)
#undef RESTRICT_THIS
#define RESTRICT_THIS __restrict
#endif

// vectorized math macros
#if !defined(NO_MANUAL_VECTORIZATION) && defined(__GNUC__) && \
    (__GNUC__ > 6) && defined(__AVX512F__)
#define USE_AVX512
#endif
#if !defined(NO_MANUAL_VECTORIZATION) && defined(__AVX__) && \
    defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#define USE_AVX
#endif

#if defined(USE_AVX) || defined(USE_AVX512)
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__)
#include <immintrin.h>  //<x86intrin.h>
#endif
#endif

namespace tinyqr {
namespace internal {

inline void print_rec_mat(const double *x, const size_t n = 4,
                          const size_t p = 4) {
  std::cout << "\n";
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < p; j++) {
      std::cout << *(x + j * n + i) << ',';
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
}

// the dumbest way to use an arena; all things allocated once at the start
template <typename scalar_t>
struct arena {
  scalar_t *storage;
  explicit arena(const size_t n_entries) {
    // use calloc for zeroes since we technicaly do need them
    storage = static_cast<scalar_t *>(calloc(n_entries, sizeof(scalar_t)));
  }
  scalar_t *operator()() { return storage; }
  scalar_t *operator+(const size_t offset) { return storage + offset; }
  ~arena() { free(storage); }
};

// since we are using 1/sqrt(x) here - default
template <typename scalar_t>
inline INLINE_THIS scalar_t inv_sqrt(scalar_t x) {
  return 1.0 / std::sqrt(x);
}
// fast inverse square root - might buy us a tiny bit, and I have been looking
// for forever to use this :)
template <>
[[maybe_unused]] inline INLINE_THIS float inv_sqrt(float x) {
  // technically omits any refinement via Newton-Raphson, but that is not
  // necessary for our purposes
  long i = *reinterpret_cast<long *>(&x);  // NOLINT [runtime/int]
  i = 0x5f3759df - (i >> 1);
  return *reinterpret_cast<float *>(&i);
}
[[maybe_unused]] inline INLINE_THIS double inv_sqrt(double x) {
  return 1.0 / std::sqrt(x);
  // taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC7827340/
  /*double halfx=0.5*x;
  long long i=*reinterpret_cast<long long*>(&x);
  i=0x5FE6ED2102DCBFDA - (i>>1);
  auto y = *reinterpret_cast<double*>(&i);
  return y = 1.50087895511633457 - halfx*y*y;
  y* = 1.50000057967625766 - halfx*y*y;
  y* = 1.5000000000002520 - halfx*y*y;
  y* = 1.5000000000000000 - halfx*y*y;*/
  // return y;
}
// |a| < |b|; the compiler emits the sign-bit mask and one compare for this,
// so a bit-pattern comparison buys nothing and would read the wrong width on
// platforms where the integer type is narrower than the floating point one
inline bool compare_abs_a_lt_b(double a, double b) {
  return std::abs(a) < std::abs(b);
}
inline bool compare_abs_a_lt_b(float a, float b) {
  return std::abs(a) < std::abs(b);
}
template <typename scalar_t>
inline INLINE_THIS std::pair<scalar_t, scalar_t> givens_rotation(
    const scalar_t a, const scalar_t b) {
  if (compare_abs_a_lt_b(a, b)) {
    const scalar_t r = a / b;
    const auto s = static_cast<scalar_t>(inv_sqrt((r * r) + 1.0));
    return {s * r, s};
  }
  const scalar_t r = b / a;
  const auto c = static_cast<scalar_t>(inv_sqrt((r * r) + 1.0));
  return {c, c * r};
}
template <typename scalar_t>
[[maybe_unused]] inline INLINE_THIS void tA_matmul_B_to_C(
    std::vector<scalar_t> &RESTRICT_THIS A,
    std::vector<scalar_t> &RESTRICT_THIS B,
    std::vector<scalar_t> &RESTRICT_THIS C, const size_t ncol) {
  for (size_t k = 0; k < ncol; k++) {
    for (size_t l = 0; l < ncol; l++) {
      scalar_t accumulator = 0.0;
      for (size_t m = 0; m < ncol; m++) {  // room for optimization; the bellow
                                           // should be easy to vectorize
        accumulator += A[k * ncol + m] * B[l * ncol + m];
      }
      C[l * ncol + k] = accumulator;
    }
  }
}
// transpose a square matrix in place
template <typename scalar_t>
inline INLINE_THIS void transpose_square(scalar_t *X, const size_t p) {
  for (size_t i = 0; i < p; i++) {
    for (size_t j = i + 1; j < p; j++) {
      // note the dereferences! we want to swap the actual data,
      // not where the pointers are pointing to;
      // simply swapping the pointers would not improve cache locality
      std::swap(*(X + (j * p) + i), *(X + (i * p) + j));
    }
  }
}
template <typename scalar_t>
inline INLINE_THIS void rotate_matrix(scalar_t *RESTRICT_THIS lower,
                                      scalar_t *RESTRICT_THIS upper,
                                      const scalar_t c, const scalar_t s,
                                      size_t p) {
  for (; p > 0; --p) {
    const scalar_t temp_1 = *lower;
    const scalar_t temp_2 = *upper;
    *lower = c * temp_1 + s * temp_2;
    *upper = -s * temp_1 + c * temp_2;
    ++lower;
    ++upper;
  }
}
#ifdef USE_AVX
template <>
inline INLINE_THIS void rotate_matrix(float *RESTRICT_THIS lower,
                                      float *RESTRICT_THIS upper, const float c,
                                      const float s, size_t p) {
  if (p > 7) {
    const __m256 c_ = _mm256_set1_ps(c);
    const __m256 s_ = _mm256_set1_ps(s);
    __m256 res = _mm256_setzero_ps();
    for (; p > 7; p -= 8) {
      // set current register values
      const auto lower_ = _mm256_loadu_ps(lower);
      const auto upper_ = _mm256_loadu_ps(upper);
      // updates on lower
      res = _mm256_add_ps(_mm256_mul_ps(c_, lower_), _mm256_mul_ps(s_, upper_));
      // store in lower
      _mm256_storeu_ps(lower, res);
      // updates on upper
      res = _mm256_sub_ps(_mm256_mul_ps(c_, upper_), _mm256_mul_ps(s_, lower_));
      // store in upper
      _mm256_storeu_ps(upper, res);
      lower += 8;
      upper += 8;
    }
  }
  for (; p > 0; --p) {
    const float temp_1 = *lower;
    const float temp_2 = *upper;
    *lower = c * temp_1 + s * temp_2;
    *upper = -s * temp_1 + c * temp_2;
    ++lower;
    ++upper;
  }
}
#endif
#ifdef USE_AVX_512
template <>
inline INLINE_THIS void rotate_matrix(float *RESTRICT_THIS lower,
                                      float *RESTRICT_THIS upper, const float c,
                                      const float s, size_t p) {
  if (p > 15) {
    const __m512 c_ = _mm512_set1_ps(c);
    const __m512 s_ = _mm512_set1_ps(s);
    __m256 res = _mm512_setzero_ps();
    for (; p > 15; p -= 16) {
      // set current register values
      const auto lower_ = _mm512_loadu_ps(lower);
      const auto upper_ = _mm512_loadu_ps(upper);
      // updates on lower
      res = _mm512_add_ps(_mm512_mul_ps(c_, lower_), _mm512_mul_ps(s_, upper_));
      // store in lower
      _mm512_storeu_ps(lower, res);
      // updates on upper
      res = _mm512_sub_ps(_mm512_mul_ps(c_, upper_), _mm512_mul_ps(s_, lower_));
      // store in upper
      _mm512_storeu_ps(upper, res);
      lower += 16;
      upper += 16;
    }
  }
  for (; p > 0; --p) {
    const float temp_1 = *lower;
    const float temp_2 = *upper;
    *lower = c * temp_1 + s * temp_2;
    *upper = -s * temp_1 + c * temp_2;
    ++lower;
    ++upper;
  }
}
#endif
template <typename scalar_t>
inline INLINE_THIS void make_identity(scalar_t *arena, const size_t n) {
  for (size_t i = 0; i < n; i++)
    *(arena + (i * n + i)) = static_cast<scalar_t>(1.0);
}
template <typename scalar_t, const bool report_success = true,
          const size_t tune = 1000>
[[maybe_unused]] void validate_qr(const std::vector<scalar_t> &RESTRICT_THIS X,
                                  const std::vector<scalar_t> &RESTRICT_THIS Q,
                                  const std::vector<scalar_t> &RESTRICT_THIS R,
                                  const size_t n, const size_t p) {
  // constant factor here added since epsilon is too small otherwise
  constexpr auto eps =
      std::numeric_limits<scalar_t>::epsilon() * static_cast<scalar_t>(tune);
  // this trick is done since some third-party, limited precision floats
  // do not provide an impl of the constexpr numeric limits epsilon function
  // const auto eps = static_cast<scalar_t>(eps_);
  // Matrix multiplication QR
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) {
      auto tmp = scalar_t(0.0);
      for (size_t k = 0; k < p; ++k) {
        tmp += Q[k * n + i] * R[j * p + k];
      }
      // Compare to original matrix X
      if (std::abs(X[j * n + i] - tmp) > eps) {
        std::cout << "Error in {validate_qr}, " << tmp << " != " << X[i * p + j]
                  << " diff: " << std::abs(X[j * n + i] - tmp)
                  << " eps: " << eps << "\n";
        std::cout << "Failed to recreate input from QR matrices for size " << n
                  << ", " << p << "\n";
        return;
      }
    }
  }
  if constexpr (report_success) {
    std::cout << "Validation of QR successful for size " << n << ", " << p
              << std::endl;
  }
}
template <typename scalar_t, const bool cleanup = false>
void qr_impl(scalar_t *RESTRICT_THIS Q_data, scalar_t *RESTRICT_THIS R_data,
             const size_t n, const size_t p, const scalar_t tol) {
  // for 'ridge' all we really need to is carry out the rotations that 'would'
  // be on the diagonal rotations only need to annihilate entries on that
  // diagonal (which are **all**) equal to lambda (or rather, sqrt(lambda), to
  // be more accurate); since we know these rotations only serve to annihilate
  // these diagonal entries
  constexpr auto min_ = std::numeric_limits<scalar_t>::min();
  // the key to optimizing this is probably to take R as R transposed - most
  // likely a lot of work is done just in the k loops, which is probably a good
  // place to optimize
  // the first iteration has j == 0; we can leverage this for **slightly**
  // better code
  for (size_t i = n - 1; i > 0; --i) {
    const auto i_p = i * p;
    const auto b = *(R_data + i_p);
    // no rotation
    if (std::abs(b) <= min_) continue;
    const auto i_p_1 = i_p - p;
    const auto [c, s] = givens_rotation(*(R_data + i_p_1), b);
    // you can make the matrix multiplication implicit, as the given's rotation
    // only impacts a moving 2x2 block
    // R is transposed
    const auto i_n = i * n;
    rotate_matrix(R_data + i_p_1, R_data + i_p, c, s, p);
    rotate_matrix(Q_data + i_n - n, Q_data + i_n, c, s, n);
  }
  for (size_t j = 1; j < p; j++) {
    for (size_t i = n - 1; i > j; --i) {
      const auto i_p = i * p;
      const auto i_p_1 = i_p - p;
      // no rotation
      if (std::abs(*(R_data + i_p + j)) <= min_) continue;
      const auto [c, s] =
          givens_rotation(*(R_data + i_p_1 + j), *(R_data + i_p + j));
      // you can make the matrix multiplication implicit, as the given's
      // rotation only impacts a moving 2x2 block R is transposed
      const auto i_n = i * n;
      rotate_matrix(R_data + i_p_1, R_data + i_p, c, s, p);
      rotate_matrix(Q_data + i_n - n, Q_data + i_n, c, s, n);
    }
  }
  // clean up R - particularly under the diagonal - only useful if you are
  // interested in the actual decomposition
  if constexpr (cleanup) {
    for (size_t j = 0; j < p; j++) {
      for (size_t i = 0; i < n; i++) {
        auto *val_ptr = (R_data + i * p + j);
        auto val = *val_ptr;
        val = std::abs(val) < tol ? 0.0 : val;
        *val_ptr = val;
      }
    }
  }
}
// Function for back substitution using QR decomposition result
// this does not require any temporaries
template <const bool is_ridge = false, typename scalar_t>
std::vector<scalar_t> back_solve(const scalar_t *RESTRICT_THIS Q,
                                 const scalar_t *RESTRICT_THIS R,
                                 const scalar_t *RESTRICT_THIS y,
                                 const size_t nrow, const size_t ncol) {
  // possibly replacing this with something else is actually not more efficient
  std::vector<scalar_t> result(ncol, 0.0);
  // offsetting in Q controlled this way
  const size_t Q_size = nrow + (is_ridge * ncol);
  for (size_t i = ncol; i-- > 0;) {
    scalar_t temp = 0.0;
    // this is actually the product with R' since that has a much nicer
    // access pattern
    // TODO(JSzitas): Technically permits FMA
    const scalar_t *view_r = R + (i * ncol);
    for (size_t j = i + 1; j < ncol; ++j) {
      temp += *(view_r + j) * result[j];
    }
    // do a multiplication by -1;
    // the result is that you save on a temporary since the next computation
    // can be directly accumulated here; possibly replacing the increments
    // above by subtractions is faster, but I assume adds/mults are a bit
    // cheaper than subtractions
    temp *= -1.;
    const scalar_t *view_q = Q + (i * Q_size);
    // TODO(JSzitas): Technically permits FMA
    for (size_t j = 0; j < nrow; ++j) {
      // product Q'y need not be computed ahead of time; we can lazily compute
      // coefficient by coefficient, requiring only one temporary
      // stack-allocated variable
      temp += *(view_q + j) * *(y + j);
    }
    result[i] = temp / *(R + (i * ncol + i));
  }
  return result;
}
// constexpr implementations of the above, plus basic functions like pow
template <typename T>
constexpr T cabs(T x) {
  return x < 0 ? -x : x;
}
template <typename T>
constexpr T cpow(T base, int exp) {
  if (exp == 0) return 1.0;  // base case for zero exponents
  if (exp < 0) return 1.0 / cpow(base, -exp);
  T result = 1.0;
  while (exp > 0) {
    if (exp % 2 == 1) {
      result *= base;
    }
    base *= base;
    exp /= 2;
  }
  return result;
}
template <typename T>
constexpr T csqrt(T x, T curr = 1.0, T prev = 0.0) {
  return curr == prev ? curr : csqrt(x, 0.5 * (curr + x / curr), curr);
}
template <typename T>
constexpr T cinv_sqrt(T x, T curr = 1.0, T prev = 0.0) {
  return 1. / csqrt(x, curr, prev);
}
template <typename T>
constexpr void cswap(T &a, T &b) {
  T temp = a;
  a = b;
  b = temp;
}
template <typename scalar_t>
constexpr std::pair<scalar_t, scalar_t> cgivens_rotation(const scalar_t a,
                                                         const scalar_t b) {
  if (cabs(b) > cabs(a)) {
    const scalar_t r = a / b;
    const auto s = static_cast<scalar_t>(cinv_sqrt(cpow(r, 2) + 1.0));
    return {s * r, s};
  }
  const scalar_t r = b / a;
  const auto c = static_cast<scalar_t>(cinv_sqrt(cpow(r, 2) + 1.0));
  return {c, c * r};
}
template <typename scalar_t>
void constexpr crotate_matrix(scalar_t *lower, scalar_t *upper,
                              const scalar_t c, const scalar_t s, size_t p) {
  for (; p > 0; --p) {
    const scalar_t temp_1 = *lower;
    const scalar_t temp_2 = *upper;
    *lower = c * temp_1 + s * temp_2;
    *upper = -s * temp_1 + c * temp_2;
    ++lower;
    ++upper;
  }
}
template <typename scalar_t, const size_t n>
constexpr std::array<scalar_t, n * n> cmake_identity() {
  std::array<scalar_t, n *n> result = {};
  for (size_t i = 0; i < n; i++) result[i * n + i] = static_cast<scalar_t>(1.0);
  return result;
}
template <typename scalar_t, const size_t n, const size_t p>
constexpr void ctranspose_rect(std::array<scalar_t, n * p> &X) {
  for (size_t i = 0; i < p; i++) {
    for (size_t j = i + 1; j < p; j++) {
      cswap(X[(j * p) + i], X[(i * p) + j]);
    }
  }
}
template <typename scalar_t, const size_t n, const size_t p,
          const bool cleanup = false>
constexpr void cqr_impl(std::array<scalar_t, n * n> &Q,
                        std::array<scalar_t, n * p> &R, const scalar_t tol) {
  // the key to optimizing this is probably to take R as R transposed - most
  // likely a lot of work is done just in the k loops, which is probably a good
  // place to optimize
  constexpr auto min_ = std::numeric_limits<scalar_t>::min();
  for (size_t j = 0; j < p; j++) {
    for (size_t i = n - 1; i > j; --i) {
      const auto i_p = i * p;
      const auto i_p_1 = i_p - p;
      // no rotation
      if (std::abs(R[i_p + j]) <= min_) continue;
      // using tuples and structured bindings should make this fairly ok
      // performance wise
      // check if R[j * n + i] - is not zero; otherwise we can skip this
      // iteration
      // if (std::abs(R[i * p + j]) <= std::numeric_limits<scalar_t>::min())
      // continue;
      const auto [c, s] = cgivens_rotation(R[i_p_1 + j], R[i_p + j]);
      // you can make the matrix multiplication implicit, as the given's
      // rotation only impacts a moving 2x2 block R is transposed
      const auto i_n = i * n;
      crotate_matrix(R.data() + i_p_1, R.data() + i_p, c, s, p);
      crotate_matrix(Q.data() + i_n - n, Q.data() + i_n, c, s, n);
    }
  }
  // clean up R - particularly under the diagonal - only useful if you are
  // interested in the actual decomposition
  if constexpr (cleanup) {
    for (auto &val : R) {
      val = cabs(val) < tol ? 0.0 : val;
    }
  }
}
template <typename scalar_t, const size_t n, const size_t p>
constexpr std::array<scalar_t, p> cback_solve(
    const std::array<scalar_t, n * p> &Q, const std::array<scalar_t, p * p> &R,
    const std::array<scalar_t, n> &y) {
  std::array<scalar_t, p> result = {};
  for (size_t i = p; i-- > 0;) {
    scalar_t temp = 0.0;
    // this might benefit from transposes
    for (size_t j = i + 1; j < p; ++j) {
      temp += R[j * p + i] * result[j];
    }
    scalar_t y_tmp = 0;
    for (size_t j = 0; j < n; ++j) {
      // product Q'y need not be computed ahead of time; we can lazily compute
      // coefficient by coefficient, requiring only one temporary
      // stack-allocated variable
      y_tmp += Q[i * n + j] * y[j];
    }
    result[i] = (y_tmp - temp) / R[i * p + i];
  }
  return result;
}
}  // namespace internal
template <typename scalar_t>
struct QR {
  std::vector<scalar_t> Q;
  std::vector<scalar_t> R;
};
template <const bool is_ridge = false, typename scalar_t>
[[maybe_unused]] void qr_decomposition(
    const scalar_t *X, tinyqr::internal::arena<scalar_t> &memory,
    const size_t n, const size_t p, scalar_t ridge_lambda = 0.,
    const scalar_t tol = 1e-8) {
  scalar_t *storage_arena = memory();
  size_t Q_offset, n_, p_;
  if constexpr (is_ridge) {
    n_ = n + p;
    p_ = p;
    Q_offset = n_ * n_;
    // initialize Q and R
    internal::make_identity<scalar_t>(storage_arena, n_);
    // initialize R as transposed
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < p; j++) {
        *(storage_arena + Q_offset + (i * p) + j) = *(X + (j * n) + i);
      }
    }
    for (size_t i = n; i < n_; i++) {
      for (size_t j = 0; j < p; j++) {
        *(storage_arena + Q_offset + ((n + j) * p) + j) = ridge_lambda;
      }
    }

  } else {
    n_ = n;
    p_ = p;
    // initialize Q and R
    internal::make_identity<scalar_t>(storage_arena, n);
    // initialize R as transposed
    Q_offset = n * n;
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < p; j++) {
        *(storage_arena + Q_offset + (i * p) + j) = *(X + (j * n) + i);
      }
    }
  }
  // technically returns transposed Q and transposed R; both should have
  // better cache locality in lm
  internal::qr_impl<scalar_t, true>(storage_arena, storage_arena + Q_offset, n_,
                                    p_, tol);
}
// this gets called by users, not by the implementation; thus it actually
// returns
template <typename scalar_t>
[[maybe_unused]] QR<scalar_t> qr_decomposition(const scalar_t *X,
                                               const size_t n, const size_t p,
                                               const scalar_t tol = 1e-8) {
  // the arena manages the memory, we just take a pointer
  internal::arena<scalar_t> arena(n * n + n * p);
  // no dangling
  scalar_t *storage = arena();
  // initialize Q and R
  internal::make_identity<scalar_t>(storage, n);
  // initialize R as transposed
  size_t Q_offset = n * n;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < p; j++) {
      *(storage + Q_offset + (i * p) + j) = *(X + (j * n) + i);
    }
  }
  // technically returns transposed Q and transposed R; both should have
  // better cache locality in lm
  internal::qr_impl<scalar_t, true>(storage, storage + Q_offset, n, p, tol);
  // transpose to return R, not R'
  internal::transpose_square(storage + Q_offset, p);
  return {
      std::vector<scalar_t>{storage, storage + (n * p)},
      std::vector<scalar_t>{storage + Q_offset, storage + Q_offset + (p * p)}};
}
template <typename scalar_t>
struct eigendecomposition {
  std::vector<scalar_t> eigenvals;
  std::vector<scalar_t> eigenvecs;
};
// Unshifted QR iteration A_{k+1} = R_k Q_k, reading eigenvalues off the
// diagonal. It converges only when every eigenvalue is real and the moduli
// are distinct, so it suits symmetric (e.g. covariance) matrices; a complex
// conjugate pair leaves a 2 by 2 block whose diagonal is not an eigenvalue.
// For a general real matrix use eigenvalues() below.
template <typename scalar_t>
[[maybe_unused]] eigendecomposition<scalar_t> qr_algorithm(
    const std::vector<scalar_t> &A, const size_t max_iter = 25,
    const scalar_t tol = 1e-8) {
  // TODO(JSzitas): Could the arena be a bit cheaper to use here than vectors?
  auto Ak = A;
  // A must-be square for this to work
  const size_t n = std::sqrt(A.size());
  std::vector<scalar_t> QQ(n * n, 0.0);
  for (size_t i = 0; i < n; i++) QQ[i * n + i] = 1.0;
  // initialize Q and R
  std::vector<scalar_t> Q(n * n, 0.0);
  for (size_t i = 0; i < n; i++) Q[i * n + i] = 1.0;
  std::vector<scalar_t> R(n * n);  // = Ak;
  std::vector<scalar_t> temp(Q.size());
  auto Q_data = Q.data();
  auto R_data = R.data();
  for (size_t i = 0; i < max_iter; i++) {
    // reset Q and R, G gets reset inside qr_impl
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < n; k++) {
        // probably a decent way to reset to a diagonal matrix
        Q[j * n + k] = static_cast<scalar_t>(k == j);
        // qr_impl takes the matrix transposed (row-major), as
        // qr_decomposition hands it over; a plain copy would factor Ak'
        R[j * n + k] = Ak[k * n + j];
      }
    }
    // call QR decomposition
    internal::qr_impl<scalar_t, false>(Q_data, R_data, n, n, tol);
    // note QR decomposition returns Rt, not R!!!
    internal::tA_matmul_B_to_C<scalar_t>(R, Q, Ak, n);
    // overwrite QQ in place
    size_t p = 0;
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < n; k++) {
        temp[p] = 0;
        for (size_t l = 0; l < n; l++) {
          // TODO(JSzitas): if QQ is transposed elsewhere we can do FMA here
          temp[p] += QQ[l * n + k] * Q[j * n + l];
        }
        p++;
      }
    }
    // write to A directly
    for (size_t k = 0; k < QQ.size(); k++) {
      QQ[k] = temp[k];
    }
  }
  // diagonal elements of Ak are eigenvalues - we can just shuffle elements of A
  // and resize
  for (size_t i = 1; i < n; i++) {
    Ak[i] = Ak[i * n + i];
  }
  Ak.resize(n);
  return {Ak, QQ};
}
template <typename scalar_t>
class [[maybe_unused]] QRSolver {
  const size_t n;
  std::vector<scalar_t> Ak, QQ, Q, R, temp, eigval;

 public:
  [[maybe_unused]] explicit QRSolver<scalar_t>(const size_t n) : n(n) {
    // TODO(JSzitas): replace with arena
    this->Ak = std::vector<scalar_t>(n * n);
    this->QQ = std::vector<scalar_t>(n * n, 0.0);
    for (size_t i = 0; i < n; i++) this->QQ[i * n + i] = 1.0;
    // initialize Q and R
    this->Q = std::vector<scalar_t>(n * n, 0.0);
    for (size_t i = 0; i < n; i++) this->Q[i * n + i] = 1.0;
    this->R = std::vector<scalar_t>(n * n);
    this->temp = std::vector<scalar_t>(n * n);
    this->eigval = std::vector<scalar_t>(n);
  }
  [[maybe_unused]] void solve(const std::vector<scalar_t> &A,
                              const size_t max_iter = 25,
                              const scalar_t tol = 1e-8) {
    this->Ak = A;
    // in case we need to reset QQ
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        this->QQ[i * n + j] = static_cast<scalar_t>(i == j);
      }
    }
    auto Q_data = Q.data();
    auto R_data = R.data();
    for (size_t i = 0; i < max_iter; i++) {
      // reset Q and R, G gets reset inside qr_impl
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < n; k++) {
          // probably a decent way to reset to a diagonal matrix
          Q[j * n + k] = static_cast<scalar_t>(k == j);
          // transposed, as qr_impl expects (see qr_algorithm)
          R[j * n + k] = Ak[k * n + j];
        }
      }
      // call QR decomposition
      internal::qr_impl<scalar_t, false>(Q_data, R_data, n, n, tol);
      // note QR decomposition returns Qt, not Q!!!
      internal::tA_matmul_B_to_C<scalar_t>(R, Q, Ak, n);
      // overwrite QQ in place
      size_t p = 0;
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < n; k++) {
          temp[p] = 0;
          for (size_t l = 0; l < n; l++) {
            temp[p] += QQ[l * n + k] * Q[j * n + l];
          }
          p++;
        }
      }
      // write to A colwise - i.e. directly
      for (size_t k = 0; k < QQ.size(); k++) {
        QQ[k] = temp[k];
      }
    }
    for (size_t i = 0; i < n; i++) {
      eigval[i] = Ak[i * n + i];
    }
  }
  [[maybe_unused]] const std::vector<scalar_t> &eigenvalues() const {
    return eigval;
  }
  [[maybe_unused]] const std::vector<scalar_t> &eigenvectors() const {
    return this->QQ;
  }
};
template <typename scalar_t>
[[maybe_unused]] std::vector<scalar_t> lm(
    const scalar_t *RESTRICT_THIS X, const scalar_t *RESTRICT_THIS y,
    const size_t nrow, const size_t ncol,
    const scalar_t ridge_lambda = std::numeric_limits<scalar_t>::min(),
    const scalar_t tol = 1e-7) {
  // lives until this function goes out of scope; don't need to manage lifetime
  // include ridge?
  const bool is_ridge =
      ridge_lambda > 10e7 * std::numeric_limits<scalar_t>::min();
  // allocates auxiliary memory for Q and R
  const size_t Q_offset =
      is_ridge ? (nrow + ncol) * (nrow + ncol) : nrow * nrow;
  size_t mem_size =
      is_ridge ? Q_offset + ((nrow * ncol) + (ncol * ncol))
               :                      // Q + R (by brackets)
          nrow * nrow + nrow * ncol;  // or both, without the ridge penalty
  internal::arena<scalar_t> arena(mem_size);
  is_ridge ? qr_decomposition<true>(X, arena, nrow, ncol, ridge_lambda, tol)
           : qr_decomposition<false>(X, arena, nrow, ncol, ridge_lambda, tol);
  // does this actually just work for ridge?
  // solve R'x = Q'y - actually still using the transposed R
  return is_ridge ? internal::back_solve<true>(arena(), (arena + Q_offset), y,
                                               nrow, ncol)
                  : internal::back_solve<false>(arena(), (arena + Q_offset), y,
                                                nrow, ncol);
}
template <typename scalar_t>
[[maybe_unused]] std::vector<scalar_t> lm(
    const std::vector<scalar_t> &RESTRICT_THIS X,
    const std::vector<scalar_t> &RESTRICT_THIS y,
    scalar_t ridge_lambda = std::numeric_limits<scalar_t>::min(),
    const scalar_t tol = 1e-7) {
  const size_t nrow = y.size();
  if (nrow == 0) return {};  // no observations: nothing to fit
  const size_t ncol = X.size() / nrow;
  return lm(X.data(), y.data(), nrow, ncol, ridge_lambda, tol);
}
// Eigenvalues of a general real square matrix.
//
// The matrix is reduced to upper Hessenberg form by Givens similarity
// transforms, then the Francis double-shift QR iteration with deflation (the
// EISPACK hqr algorithm) drives it to real Schur form; real eigenvalues are
// read off 1 by 1 diagonal blocks and complex conjugate pairs off the 2 by 2
// blocks that remain.
// [[CITATION]] Golub and Van Loan, Matrix Computations, 4th ed., sections 7.4
// (Hessenberg reduction) and 7.5 (the practical QR algorithm).
namespace internal {
// In place similarity reduction of the column-major n by n matrix `a` to
// upper Hessenberg form (zero below the first subdiagonal). Within column k
// the entries are annihilated from the bottom up by rotations in the plane
// (i-1, i); each rotation is applied to the rows and to the columns, so the
// spectrum is unchanged. Rows above k+1 and columns up to k hold zeros already
// and are not touched.
template <typename scalar_t>
void hessenberg_reduce(scalar_t *a, const size_t n) {
  constexpr auto min_ = std::numeric_limits<scalar_t>::min();
  auto at = [a, n](const size_t i, const size_t j) -> scalar_t & {
    return a[j * n + i];
  };
  for (size_t k = 0; k + 2 < n; ++k) {
    for (size_t i = n - 1; i > k + 1; --i) {
      const scalar_t b = at(i, k);
      if (std::abs(b) <= min_) continue;
      const auto [c, s] = givens_rotation(at(i - 1, k), b);
      for (size_t j = k; j < n; ++j) {
        const scalar_t lower = at(i - 1, j);
        const scalar_t upper = at(i, j);
        at(i - 1, j) = c * lower + s * upper;
        at(i, j) = -s * lower + c * upper;
      }
      at(i, k) = 0;
      for (size_t r = 0; r < n; ++r) {
        const scalar_t left = at(r, i - 1);
        const scalar_t right = at(r, i);
        at(r, i - 1) = c * left + s * right;
        at(r, i) = -s * left + c * right;
      }
    }
  }
}
// Eigenvalues of the upper Hessenberg column-major n by n matrix `a`, which
// is overwritten. Throws when an eigenvalue has not deflated after `max_iter`
// double steps; exceptional shifts at steps 10 and 20 break the rare cycles
// that stall the standard shift.
template <typename scalar_t>
void hessenberg_eigenvalues(scalar_t *a, const size_t n, const size_t max_iter,
                            std::vector<std::complex<scalar_t>> &result) {
  auto at = [a, n](const size_t i, const size_t j) -> scalar_t & {
    return a[j * n + i];
  };
  result.assign(n, std::complex<scalar_t>(0, 0));
  // scale of the matrix, used in the negligibility test when the two
  // diagonal entries next to a subdiagonal entry are both zero
  scalar_t anorm = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = (i == 0 ? 0 : i - 1); j < n; ++j) {
      anorm += std::abs(at(i, j));
    }
  }
  scalar_t t = 0;  // sum of the exceptional shifts applied so far
  // rows and columns 0 .. n_active - 1 still hold undeflated eigenvalues
  size_t n_active = n;
  while (n_active > 0) {
    const size_t nn = n_active - 1;
    size_t its = 0;
    scalar_t x = 0, y = 0, w = 0;
    while (true) {
      // l is the start of the trailing unreduced block: the lowest row whose
      // subdiagonal entry is negligible, or 0
      size_t l = 0;
      for (size_t ll = nn; ll > 0; --ll) {
        scalar_t s = std::abs(at(ll - 1, ll - 1)) + std::abs(at(ll, ll));
        if (s == 0) s = anorm;
        if (std::abs(at(ll, ll - 1)) + s == s) {
          at(ll, ll - 1) = 0;
          l = ll;
          break;
        }
      }
      x = at(nn, nn);
      if (l == nn) {  // a 1 by 1 block: one real eigenvalue
        result[nn] = {x + t, 0};
        n_active -= 1;
        break;
      }
      y = at(nn - 1, nn - 1);
      w = at(nn, nn - 1) * at(nn - 1, nn);
      if (l + 1 == nn) {  // a 2 by 2 block: a real or a complex pair
        const scalar_t p = static_cast<scalar_t>(0.5) * (y - x);
        const scalar_t q = p * p + w;
        scalar_t z = std::sqrt(std::abs(q));
        x += t;
        if (q >= 0) {
          z = p + std::copysign(z, p);
          result[nn - 1] = result[nn] = {x + z, 0};
          if (z != 0) result[nn] = {x - w / z, 0};
        } else {
          result[nn - 1] = {x + p, -z};
          result[nn] = {x + p, z};
        }
        n_active -= 2;
        break;
      }
      if (its == max_iter) {
        throw std::runtime_error(
            "tinyqr::eigenvalues: no convergence after max_iter QR steps");
      }
      if (its == 10 || its == 20) {  // exceptional shift
        t += x;
        for (size_t i = 0; i <= nn; ++i) at(i, i) -= x;
        const scalar_t s =
            std::abs(at(nn, nn - 1)) + std::abs(at(nn - 1, nn - 2));
        y = x = static_cast<scalar_t>(0.75) * s;
        w = static_cast<scalar_t>(-0.4375) * s * s;
      }
      ++its;
      // form the double shift and look for two consecutive small subdiagonal
      // entries, which let the step start below row l
      size_t m = nn - 2;
      scalar_t p = 0, q = 0, r = 0, z = 0;
      while (true) {
        z = at(m, m);
        r = x - z;
        scalar_t s = y - z;
        p = (r * s - w) / at(m + 1, m) + at(m, m + 1);
        q = at(m + 1, m + 1) - z - r - s;
        r = at(m + 2, m + 1);
        s = std::abs(p) + std::abs(q) + std::abs(r);
        p /= s;
        q /= s;
        r /= s;
        if (m == l) break;
        const scalar_t u = std::abs(at(m, m - 1)) * (std::abs(q) + std::abs(r));
        const scalar_t v =
            std::abs(p) * (std::abs(at(m - 1, m - 1)) + std::abs(z) +
                           std::abs(at(m + 1, m + 1)));
        if (u + v == v) break;
        --m;
      }
      for (size_t i = m + 2; i <= nn; ++i) {
        at(i, i - 2) = 0;
        if (i != m + 2) at(i, i - 3) = 0;
      }
      // the double QR step on rows l..nn and columns m..nn, as a sequence of
      // Householder reflections of order 3 (order 2 for the last one)
      for (size_t k = m; k < nn; ++k) {
        if (k != m) {
          p = at(k, k - 1);
          q = at(k + 1, k - 1);
          r = 0;
          if (k + 1 != nn) r = at(k + 2, k - 1);
          x = std::abs(p) + std::abs(q) + std::abs(r);
          if (x != 0) {
            p /= x;
            q /= x;
            r /= x;
          }
        }
        const scalar_t s = std::copysign(std::sqrt(p * p + q * q + r * r), p);
        if (s == 0) continue;
        if (k == m) {
          if (l != m) at(k, k - 1) = -at(k, k - 1);
        } else {
          at(k, k - 1) = -s * x;
        }
        p += s;
        x = p / s;
        y = q / s;
        z = r / s;
        q /= p;
        r /= p;
        for (size_t j = k; j <= nn; ++j) {  // rows k, k+1, k+2
          p = at(k, j) + q * at(k + 1, j);
          if (k + 1 != nn) {
            p += r * at(k + 2, j);
            at(k + 2, j) -= p * z;
          }
          at(k + 1, j) -= p * y;
          at(k, j) -= p * x;
        }
        const size_t mmin = nn < k + 3 ? nn : k + 3;
        for (size_t i = l; i <= mmin; ++i) {  // columns k, k+1, k+2
          p = x * at(i, k) + y * at(i, k + 1);
          if (k + 1 != nn) {
            p += z * at(i, k + 2);
            at(i, k + 2) -= p * r;
          }
          at(i, k + 1) -= p * q;
          at(i, k) -= p;
        }
      }
    }
  }
}
}  // namespace internal
// Eigenvalues of the real column-major n by n matrix `a`, in no particular
// order, written to `out`; `a` is overwritten and nothing is allocated once
// `out` has capacity n. `max_iter` bounds the QR steps spent on any one
// eigenvalue; the algorithm needs a handful per eigenvalue, so exceeding it
// signals a matrix the method cannot handle (NaN or Inf entries, typically).
template <typename scalar_t>
[[maybe_unused]] void eigenvalues(scalar_t *a, const size_t n,
                                  std::vector<std::complex<scalar_t>> &out,
                                  const size_t max_iter = 30) {
  internal::hessenberg_reduce(a, n);
  internal::hessenberg_eigenvalues(a, n, max_iter, out);
}
// The same on a copy of A, which is left untouched.
template <typename scalar_t>
[[maybe_unused]] std::vector<std::complex<scalar_t>> eigenvalues(
    const scalar_t *A, const size_t n, const size_t max_iter = 30) {
  std::vector<scalar_t> a(A, A + n * n);
  std::vector<std::complex<scalar_t>> out;
  eigenvalues(a.data(), n, out, max_iter);
  return out;
}
template <typename scalar_t>
[[maybe_unused]] std::vector<std::complex<scalar_t>> eigenvalues(
    const std::vector<scalar_t> &A, const size_t max_iter = 30) {
  const auto n = static_cast<size_t>(std::sqrt(A.size()));
  return eigenvalues(A.data(), n, max_iter);
}
template <typename scalar_t, const size_t n, const size_t p>
struct cQR {
  std::array<scalar_t, n * p> Q;
  std::array<scalar_t, p * p> R;
};
template <typename scalar_t, const size_t n, const size_t p>
constexpr cQR<scalar_t, n, p> cqr_decomposition(
    const std::array<scalar_t, n * p> &X, const scalar_t tol = 1e-8) {
  // initialize Q and R
  std::array<scalar_t, n *n> Q = internal::cmake_identity<scalar_t, n>();
  // initialize R as transposed
  std::array<scalar_t, n *p> R = {};
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < p; j++) {
      R[i * p + j] = X[j * n + i];
    }
  }
  internal::cqr_impl<scalar_t, n, p, true>(Q, R, tol);

  // we do not need to manipulate more than pxp block of R
  internal::ctranspose_rect<scalar_t, n, p>(R);
  size_t k = 0;
  std::array<scalar_t, n *p> Q_ = {};
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < p; j++) {
      Q_[i * p + j] = Q[k++];
    }
  }
  k = 0;
  std::array<scalar_t, p *p> R_ = {};
  for (size_t i = 0; i < p; i++) {
    for (size_t j = 0; j < p; j++) {
      R_[i * p + j] = R[k++];
    }
  }
  return {Q_, R_};
}
template <typename scalar_t, const size_t n, const size_t p>
[[maybe_unused]] constexpr std::array<scalar_t, p> clm(
    const std::array<scalar_t, n * p> &X, const std::array<scalar_t, n> &y,
    const scalar_t tol = 1e-12) {
  // compute QR decomposition
  const auto qr = tinyqr::cqr_decomposition<scalar_t, n, p>(X, tol);
  // solve Rx = Q'y
  return internal::cback_solve<scalar_t, n, p>(qr.Q, qr.R, y);
}
}  // namespace tinyqr
#endif  // TINYQR_H_"
