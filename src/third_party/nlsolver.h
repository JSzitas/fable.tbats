// Nonlinear optimization in C++
// https://github.com/JSzitas/nlsolver
//
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// Copyright (c) 2023- Juraj Szitas
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

#ifndef NLSOLVER_H_
#define NLSOLVER_H_

#if defined(__clang__)
#pragma clang diagnostic push
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "./tinyqr.h"

// TODO(JSzitas): delete
template <typename T>
void print_vector(T x) {
  for (auto &val : x) {
    std::cout << val << ",";
  }
  std::cout << "\n";
}
// mostly dot products and other fun vector math stuff
namespace nlsolver::math {
template <typename T>
[[maybe_unused]] inline T dot(const T *x, const T *y, int f) {
  T s = 0;
  for (int z = 0; z < f; z++) {
    s += (*x) * (*y);
    x++;
    y++;
  }
  return s;
}

template <typename T>
[[maybe_unused]] inline T fast_sum(const T *x, int size) {
  T result = 0;
  for (int i = 0; i < size; i++) {
    result += (*x);
    x++;
  }
  return result;
}

template <typename T>
[[maybe_unused]] inline T vec_scalar_mult(const T *vec, const T *scalar,
                                          int f) {
  T result = 0;
  // load single scalar
  // Don't forget the remaining values.
  for (int i = 0; i < f; i++) {
    result += *vec * *scalar;
    vec++;
  }
  return result;
}
template <typename T>
[[maybe_unused]] inline T norm(const T *x, int f) {
  T s = 0;
  for (int z = 0; z < f; z++) {
    s += (*x) * (*x);
    x++;
  }
  return std::sqrt(s);
}
template <typename T>
[[maybe_unused]] inline T norm_diff(const T *x, const T *y, int f) {
  T s = 0;
  for (int z = 0; z < f; z++) {
    T tmp = *x - *y;
    s += tmp * tmp;
    x++;
    y++;
  }
  return std::sqrt(s);
}
template <typename T>
[[maybe_unused]] inline void a_plus_b(T *a, const T *b, int f) {
  for (int i = 0; i < f; i++) {
    *a += *b;
    a++;
    b++;
  }
}
template <typename T>
[[maybe_unused]] inline void a_minus_b(T *a, const T *b, int f) {
  for (int i = 0; i < f; i++) {
    *a -= *b;
    a++;
    b++;
  }
}
template <typename T>
[[maybe_unused]] inline void a_plus_b_to_c(const T *a, const T *b, T *c,
                                           int f) {
  for (int i = 0; i < f; i++) {
    *c = *a + *b;
    a++;
    b++;
    c++;
  }
}

template <typename T>
[[maybe_unused]] inline void a_minus_b_to_c(const T *a, const T *b, T *c,
                                            int f) {
  for (int i = 0; i < f; i++) {
    *c = *a - *b;
    a++;
    b++;
    c++;
  }
}

template <typename T>
[[maybe_unused]] inline T sum_a_plus_b_times_c(const T *a, const T *b,
                                               const T *c, int f) {
  T result = 0;
  for (int i = 0; i < f; i++) {
    result += (*a + *b) * (*c);
    a++;
    b++;
    c++;
  }
  return result;
}

template <typename T>
[[maybe_unused]] inline T sum_a_minus_b_times_c(const T *a, const T *b,
                                                const T *c, int f) {
  T result = 0;
  for (int i = 0; i < f; i++) {
    result += (*a - *b) * (*c);
    a++;
    b++;
    c++;
  }
  return result;
}
template <typename T>
[[maybe_unused]] inline void a_plus_scalar_to_b(const T *a, const T scalar,
                                                T *b, int f) {
  for (int i = 0; i < f; i++) {
    *b = (*a + scalar);
    a++;
    b++;
  }
}

template <typename T>
[[maybe_unused]] inline void a_minus_scalar_to_b(const T *a, const T scalar,
                                                 T *b, int f) {
  for (int i = 0; i < f; i++) {
    *b = (*a - scalar);
    a++;
    b++;
  }
}
template <typename T>
[[maybe_unused]] inline void a_mult_scalar_to_b(const T *a, const T scalar,
                                                T *b, int f) {
  for (int i = 0; i < f; i++) {
    *b = (*a * scalar);
    a++;
    b++;
  }
}
template <typename T>
[[maybe_unused]] inline void a_mult_scalar_add_b(const T *a, const T scalar,
                                                 T *b, int f) {
  for (int i = 0; i < f; i++) {
    *b += (*a * scalar);
    a++;
    b++;
  }
}
template <typename T>
[[maybe_unused]] inline void a_minus_b_mult_scalar_add_c(const T *a, const T *b,
                                                         const T scalar, T *c,
                                                         int f) {
  for (int i = 0; i < f; i++) {
    *c += (*a - *b) * scalar;
    a++;
    b++;
    c++;
  }
}
template <typename T>
[[maybe_unused]] inline void a_mul_scalar(T *a, const T scalar, int f) {
  for (int i = 0; i < f; i++) {
    *a *= scalar;
    a++;
  }
}
// this is best defined here even though it is not technically a function
// we probably want to reuse much
template <typename T>
[[maybe_unused]] inline void hessian_update_inner_loop(
    T *inv_hessian, const T *step, const T *grad_diff_inv_hess, const T rho,
    const T denom, const int n_dim) {
  for (int j = 0; j < n_dim; j++) {
    for (int i = 0; i < n_dim; i++) {
      // do not replace this with -= or the whole thing falls apart
      // because of operator order precedence - e.g. whole rhs would
      // get evaluated before -=, whereas we want to do inv_hessian - first part
      // + second part
      *(inv_hessian + j * n_dim + i) =
          *(inv_hessian + j * n_dim + i) -
          // first part
          rho * (*(step + i) * *(grad_diff_inv_hess + j) +
                 *(grad_diff_inv_hess + i) * *(step + j) +
                 // second part ->  multiply(step[i], denom * step[j])
                 denom * *(step + i) * *(step + j));
    }
  }
}
template <typename scalar_t>
void cholesky(std::vector<scalar_t> &A) {
  const auto n = static_cast<size_t>(std::sqrt(A.size()));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      scalar_t sum = 0;
      for (size_t k = 0; k < j; ++k) {
        sum += A[i * n + k] * A[j * n + k];
      }
      A[i * n + j] = (1.0 / A[j * n + j] * (A[i * n + j] - sum));
    }
    // diagonal elements only
    scalar_t sum = 0;
    for (size_t k = 0; k < i; ++k) {
      sum += A[i * n + k] * A[i * n + k];
    }
    A[i * n + i] = sqrt(A[i * n + i] - sum);
  }
}
template <typename scalar_t>
void backsolve_inplace_t(std::vector<scalar_t> &U, std::vector<scalar_t> &b,
                         const size_t n) {
  int i = static_cast<int>(n) - 1;
  for (; i >= 0; --i) {
    scalar_t sum = 0.0;
    for (size_t j = i + 1; j < n; ++j) {
      sum += U[j * n + i] * b[j];
    }
    b[i] = (b[i] - sum) / U[i * n + i];
  }
}
template <typename scalar_t = double>
void forwardsolve_inplace(std::vector<scalar_t> &update,
                          const std::vector<scalar_t> &L,
                          const std::vector<scalar_t> &b, const size_t n) {
  std::fill(update.begin(), update.end(), 0.0);
  for (size_t i = 0; i < n; ++i) {
    scalar_t sum = 0.0;
    for (size_t j = 0; j < i; ++j) {
      sum += L[i * n + j] * update[j];
    }
    update[i] = (b[i] - sum) / L[i + i * n];
  }
}
template <typename T>
bool is_diagonal(const std::vector<T> &A) {
  const size_t size = static_cast<size_t>(std::sqrt(A.size()));
  // tolerance relative to the largest-magnitude diagonal entry (the previous
  // check used a one-sided `>` that ignored negative off-diagonals, with an
  // absolute threshold of ~2.2e-4 that was far too loose)
  T scale = 0;
  for (size_t i = 0; i < size; ++i) {
    scale = std::max(scale, std::abs(A[i * size + i]));
  }
  const T tol =
      (scale > 0 ? scale : T(1)) * std::numeric_limits<T>::epsilon() * 100;
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      if (i != j && std::abs(A[i * size + j]) > tol) {
        return false;
      }
    }
  }
  return true;
}

template <typename scalar_t>
void get_update_with_hessian(std::vector<scalar_t> &update,
                             std::vector<scalar_t> &hess,
                             std::vector<scalar_t> &grad) {
  // if hessian is diagonal we have a much faster update - since the inverse
  // is just 1/diagonal entry, and the products are just simple multiplications
  const size_t n = grad.size();
  if (is_diagonal(hess)) {
    // much faster update loop
    for (size_t i = 0; i < n; i++) {
      update[i] = grad[i] / hess[i * n + i];
    }
    return;
  }
  // note that this overwrites the hessian due to the cholesky decomposition
  // being in-place
  // compute cholesky of hessian
  cholesky(hess);
  // use cholesky to solve system of equations
  forwardsolve_inplace(update, hess, grad, n);
  backsolve_inplace_t(hess, update, n);
}

// inspired by annoylib, see
// https://github.com/spotify/annoy/blob/main/src/annoylib.h
#if !defined(NO_MANUAL_VECTORIZATION) && defined(__GNUC__) && \
    (__GNUC__ > 6) && defined(__AVX512F__)
#define DOT_USE_AVX512
#endif
#if !defined(NO_MANUAL_VECTORIZATION) && defined(__AVX__) && \
    defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#define DOT_USE_AVX
#endif

#if defined(DOT_USE_AVX) || defined(DOT_USE_AVX512)
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__)
#include <immintrin.h>  //<x86intrin.h>
#endif
#endif

#ifdef DOT_USE_AVX
// Horizontal single sum of 256bit vector.
inline float hsum256_ps_avx(__m256 v) {
  const __m128 x128 =
      _mm_add_ps(_mm256_extractf128_ps(v, 1), _mm256_castps256_ps128(v));
  const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
  const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
  return _mm_cvtss_f32(x32);
}

inline double hsum256_pd_avx(__m256d v) {
  __m128d vlow = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);  // high 128
  vlow = _mm_add_pd(vlow, vhigh);               // reduce down to 128
  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

template <>
[[maybe_unused]] inline float dot<float>(const float *x, const float *y,
                                         int f) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    acc = _mm256_add_ps(
        acc, _mm256_mul_ps(_mm256_loadu_ps(x + i), _mm256_loadu_ps(y + i)));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < f; i++) result += x[i] * y[i];
  return result;
}

template <>
[[maybe_unused]] inline float fast_sum<float>(const float *x, int size) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  for (; i + 8 <= size; i += 8) {
    acc = _mm256_add_ps(acc, _mm256_loadu_ps(x + i));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < size; i++) result += x[i];
  return result;
}

template <>
[[maybe_unused]] inline float vec_scalar_mult<float>(const float *vec,
                                                     const float *scalar,
                                                     int f) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  const __m256 S = _mm256_set1_ps(*scalar);
  for (; i + 8 <= f; i += 8) {
    acc = _mm256_add_ps(acc, _mm256_mul_ps(_mm256_loadu_ps(vec + i), S));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < f; i++) result += vec[i] * (*scalar);
  return result;
}

template <>
[[maybe_unused]] inline float norm<float>(const float *x, int f) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    acc = _mm256_add_ps(
        acc, _mm256_mul_ps(_mm256_loadu_ps(x + i), _mm256_loadu_ps(x + i)));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < f; i++) result += x[i] * x[i];
  return std::sqrt(result);
}

template <>
[[maybe_unused]] inline float norm_diff<float>(const float *x, const float *y,
                                               int f) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    const __m256 d =
        _mm256_sub_ps(_mm256_loadu_ps(x + i), _mm256_loadu_ps(y + i));
    acc = _mm256_add_ps(acc, _mm256_mul_ps(d, d));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < f; i++) result += (x[i] - y[i]) * (x[i] - y[i]);
  return std::sqrt(result);
}

template <>
[[maybe_unused]] inline float sum_a_plus_b_times_c<float>(const float *a,
                                                          const float *b,
                                                          const float *c,
                                                          int f) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    acc = _mm256_add_ps(
        acc, _mm256_mul_ps(
                 _mm256_add_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)),
                 _mm256_loadu_ps(c + i)));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < f; i++) result += (a[i] + b[i]) * c[i];
  return result;
}

template <>
[[maybe_unused]] inline float sum_a_minus_b_times_c<float>(const float *a,
                                                           const float *b,
                                                           const float *c,
                                                           int f) {
  __m256 acc = _mm256_setzero_ps();
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    acc = _mm256_add_ps(
        acc, _mm256_mul_ps(
                 _mm256_sub_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)),
                 _mm256_loadu_ps(c + i)));
  }
  float result = hsum256_ps_avx(acc);
  for (; i < f; i++) result += (a[i] - b[i]) * c[i];
  return result;
}

template <>
[[maybe_unused]] inline void a_plus_b(float *a, const float *b, int f) {
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(
        a + i, _mm256_add_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)));
  }
  for (; i < f; i++) a[i] += b[i];
}

template <>
[[maybe_unused]] inline void a_minus_b(float *a, const float *b, int f) {
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(
        a + i, _mm256_sub_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)));
  }
  for (; i < f; i++) a[i] -= b[i];
}

template <>
[[maybe_unused]] inline void a_plus_b_to_c(const float *a, const float *b,
                                           float *c, int f) {
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(
        c + i, _mm256_add_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)));
  }
  for (; i < f; i++) c[i] = a[i] + b[i];
}

template <>
[[maybe_unused]] inline void a_minus_b_to_c(const float *a, const float *b,
                                            float *c, int f) {
  int i = 0;
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(
        c + i, _mm256_sub_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)));
  }
  for (; i < f; i++) c[i] = a[i] - b[i];
}

template <>
[[maybe_unused]] inline void a_plus_scalar_to_b(const float *a,
                                                const float scalar, float *b,
                                                int f) {
  int i = 0;
  const __m256 S = _mm256_set1_ps(scalar);
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(b + i, _mm256_add_ps(_mm256_loadu_ps(a + i), S));
  }
  for (; i < f; i++) b[i] = a[i] + scalar;
}

template <>
[[maybe_unused]] inline void a_minus_scalar_to_b(const float *a,
                                                 const float scalar, float *b,
                                                 int f) {
  int i = 0;
  const __m256 S = _mm256_set1_ps(scalar);
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(b + i, _mm256_sub_ps(_mm256_loadu_ps(a + i), S));
  }
  for (; i < f; i++) b[i] = a[i] - scalar;
}

template <>
[[maybe_unused]] inline void a_mult_scalar_to_b(const float *a,
                                                const float scalar, float *b,
                                                int f) {
  int i = 0;
  const __m256 S = _mm256_set1_ps(scalar);
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(b + i, _mm256_mul_ps(_mm256_loadu_ps(a + i), S));
  }
  for (; i < f; i++) b[i] = a[i] * scalar;
}

template <>
[[maybe_unused]] inline void a_mult_scalar_add_b(const float *a,
                                                 const float scalar, float *b,
                                                 int f) {
  int i = 0;
  const __m256 S = _mm256_set1_ps(scalar);
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(b + i,
                     _mm256_add_ps(_mm256_loadu_ps(b + i),
                                   _mm256_mul_ps(_mm256_loadu_ps(a + i), S)));
  }
  for (; i < f; i++) b[i] += a[i] * scalar;
}

template <>
[[maybe_unused]] inline void a_minus_b_mult_scalar_add_c(const float *a,
                                                         const float *b,
                                                         const float scalar,
                                                         float *c, int f) {
  int i = 0;
  const __m256 S = _mm256_set1_ps(scalar);
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(
        c + i,
        _mm256_add_ps(_mm256_loadu_ps(c + i),
                      _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(a + i),
                                                  _mm256_loadu_ps(b + i)),
                                    S)));
  }
  for (; i < f; i++) c[i] += (a[i] - b[i]) * scalar;
}

template <>
[[maybe_unused]] inline void a_mul_scalar(float *a, const float scalar, int f) {
  int i = 0;
  const __m256 S = _mm256_set1_ps(scalar);
  for (; i + 8 <= f; i += 8) {
    _mm256_storeu_ps(a + i, _mm256_mul_ps(_mm256_loadu_ps(a + i), S));
  }
  for (; i < f; i++) a[i] *= scalar;
}

template <>
[[maybe_unused]] inline double dot<double>(const double *x, const double *y,
                                           int f) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    acc = _mm256_add_pd(
        acc, _mm256_mul_pd(_mm256_loadu_pd(x + i), _mm256_loadu_pd(y + i)));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < f; i++) result += x[i] * y[i];
  return result;
}

template <>
[[maybe_unused]] inline double fast_sum<double>(const double *x, int size) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= size; i += 4) {
    acc = _mm256_add_pd(acc, _mm256_loadu_pd(x + i));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < size; i++) result += x[i];
  return result;
}

template <>
[[maybe_unused]] inline double vec_scalar_mult<double>(const double *vec,
                                                       const double *scalar,
                                                       int f) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  const __m256d S = _mm256_set1_pd(*scalar);
  for (; i + 4 <= f; i += 4) {
    acc = _mm256_add_pd(acc, _mm256_mul_pd(_mm256_loadu_pd(vec + i), S));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < f; i++) result += vec[i] * (*scalar);
  return result;
}

template <>
[[maybe_unused]] inline double norm<double>(const double *x, int f) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    acc = _mm256_add_pd(
        acc, _mm256_mul_pd(_mm256_loadu_pd(x + i), _mm256_loadu_pd(x + i)));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < f; i++) result += x[i] * x[i];
  return std::sqrt(result);
}

template <>
[[maybe_unused]] inline double norm_diff<double>(const double *x,
                                                 const double *y, int f) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    const __m256d d =
        _mm256_sub_pd(_mm256_loadu_pd(x + i), _mm256_loadu_pd(y + i));
    acc = _mm256_add_pd(acc, _mm256_mul_pd(d, d));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < f; i++) result += (x[i] - y[i]) * (x[i] - y[i]);
  return std::sqrt(result);
}

template <>
[[maybe_unused]] inline double sum_a_plus_b_times_c<double>(const double *a,
                                                            const double *b,
                                                            const double *c,
                                                            int f) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    acc = _mm256_add_pd(
        acc, _mm256_mul_pd(
                 _mm256_add_pd(_mm256_loadu_pd(a + i), _mm256_loadu_pd(b + i)),
                 _mm256_loadu_pd(c + i)));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < f; i++) result += (a[i] + b[i]) * c[i];
  return result;
}

template <>
[[maybe_unused]] inline double sum_a_minus_b_times_c<double>(const double *a,
                                                             const double *b,
                                                             const double *c,
                                                             int f) {
  __m256d acc = _mm256_setzero_pd();
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    acc = _mm256_add_pd(
        acc, _mm256_mul_pd(
                 _mm256_sub_pd(_mm256_loadu_pd(a + i), _mm256_loadu_pd(b + i)),
                 _mm256_loadu_pd(c + i)));
  }
  double result = hsum256_pd_avx(acc);
  for (; i < f; i++) result += (a[i] - b[i]) * c[i];
  return result;
}

template <>
[[maybe_unused]] inline void a_plus_b(double *a, const double *b, int f) {
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(
        a + i, _mm256_add_pd(_mm256_loadu_pd(a + i), _mm256_loadu_pd(b + i)));
  }
  for (; i < f; i++) a[i] += b[i];
}

template <>
[[maybe_unused]] inline void a_minus_b(double *a, const double *b, int f) {
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(
        a + i, _mm256_sub_pd(_mm256_loadu_pd(a + i), _mm256_loadu_pd(b + i)));
  }
  for (; i < f; i++) a[i] -= b[i];
}

template <>
[[maybe_unused]] inline void a_plus_b_to_c(const double *a, const double *b,
                                           double *c, int f) {
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(
        c + i, _mm256_add_pd(_mm256_loadu_pd(a + i), _mm256_loadu_pd(b + i)));
  }
  for (; i < f; i++) c[i] = a[i] + b[i];
}

template <>
[[maybe_unused]] inline void a_minus_b_to_c(const double *a, const double *b,
                                            double *c, int f) {
  int i = 0;
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(
        c + i, _mm256_sub_pd(_mm256_loadu_pd(a + i), _mm256_loadu_pd(b + i)));
  }
  for (; i < f; i++) c[i] = a[i] - b[i];
}

template <>
[[maybe_unused]] inline void a_plus_scalar_to_b(const double *a,
                                                const double scalar, double *b,
                                                int f) {
  int i = 0;
  const __m256d S = _mm256_set1_pd(scalar);
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(b + i, _mm256_add_pd(_mm256_loadu_pd(a + i), S));
  }
  for (; i < f; i++) b[i] = a[i] + scalar;
}

template <>
[[maybe_unused]] inline void a_minus_scalar_to_b(const double *a,
                                                 const double scalar, double *b,
                                                 int f) {
  int i = 0;
  const __m256d S = _mm256_set1_pd(scalar);
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(b + i, _mm256_sub_pd(_mm256_loadu_pd(a + i), S));
  }
  for (; i < f; i++) b[i] = a[i] - scalar;
}

template <>
[[maybe_unused]] inline void a_mult_scalar_to_b(const double *a,
                                                const double scalar, double *b,
                                                int f) {
  int i = 0;
  const __m256d S = _mm256_set1_pd(scalar);
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(b + i, _mm256_mul_pd(_mm256_loadu_pd(a + i), S));
  }
  for (; i < f; i++) b[i] = a[i] * scalar;
}

template <>
[[maybe_unused]] inline void a_mult_scalar_add_b(const double *a,
                                                 const double scalar, double *b,
                                                 int f) {
  int i = 0;
  const __m256d S = _mm256_set1_pd(scalar);
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(b + i,
                     _mm256_add_pd(_mm256_loadu_pd(b + i),
                                   _mm256_mul_pd(_mm256_loadu_pd(a + i), S)));
  }
  for (; i < f; i++) b[i] += a[i] * scalar;
}

template <>
[[maybe_unused]] inline void a_minus_b_mult_scalar_add_c(const double *a,
                                                         const double *b,
                                                         const double scalar,
                                                         double *c, int f) {
  int i = 0;
  const __m256d S = _mm256_set1_pd(scalar);
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(
        c + i,
        _mm256_add_pd(_mm256_loadu_pd(c + i),
                      _mm256_mul_pd(_mm256_sub_pd(_mm256_loadu_pd(a + i),
                                                  _mm256_loadu_pd(b + i)),
                                    S)));
  }
  for (; i < f; i++) c[i] += (a[i] - b[i]) * scalar;
}

template <>
[[maybe_unused]] inline void a_mul_scalar(double *a, const double scalar,
                                          int f) {
  int i = 0;
  const __m256d S = _mm256_set1_pd(scalar);
  for (; i + 4 <= f; i += 4) {
    _mm256_storeu_pd(a + i, _mm256_mul_pd(_mm256_loadu_pd(a + i), S));
  }
  for (; i < f; i++) a[i] *= scalar;
}

#endif
}  // namespace nlsolver::math
namespace nlsolver::rng {
#define MAX_SIZE_64_BIT_UINT (18446744073709551615U)
template <typename scalar_t = float>
struct [[maybe_unused]] halton {
  explicit halton<scalar_t>(const scalar_t base = 2)
      : b(base), y(1), n(0), d(1), x(1) {}
  scalar_t yield() {
    x = d - n;
    if (x == 1) {
      n = 1;
      d *= b;
    } else {
      y = d;
      while (x <= y) {
        y /= b;
        n = (b + 1) * y - x;
      }
    }
    return (scalar_t)(n / d);
  }
  scalar_t operator()() { return this->yield(); }
  [[maybe_unused]] void reset() {
    b = 2;
    y = 1;
    n = 0;
    d = 1;
    x = 1;
  }
  [[maybe_unused]] std::vector<scalar_t> get_state() const {
    std::vector<scalar_t> result(5);
    result[0] = b;
    result[1] = y;
    result[2] = n;
    result[3] = d;
    result[4] = x;
    return result;
  }
  [[maybe_unused]] void set_state(const scalar_t b_, const scalar_t y_,
                                  const scalar_t n_, const scalar_t d_,
                                  const scalar_t x_) {
    this->b = b_;
    this->y = y_;
    this->n = n_;
    this->d = d_;
    this->x = x_;
  }

 private:
  scalar_t b, y, n, d, x;
};

template <typename scalar_t = float>
struct [[maybe_unused]] recurrent {
  recurrent<scalar_t>() : seed_(0.5), alpha_(0.618034), z_(alpha_ + seed_) {
    this->z_ -= static_cast<scalar_t>(static_cast<uint64_t>(this->z_));
  }
  [[maybe_unused]] explicit recurrent(scalar_t seed)
      : seed_(seed), alpha_(0.618034), z_(alpha_ + seed_) {
    this->z_ -= static_cast<scalar_t>(static_cast<uint64_t>(this->z_));
  }
  scalar_t yield() {
    this->z_ += this->alpha_;
    // a slightly evil way to do z % 1 with floats
    this->z_ -= static_cast<scalar_t>(static_cast<uint64_t>(this->z_));
    return this->z_;
  }
  scalar_t operator()() { return this->yield(); }
  [[maybe_unused]] void reset() {
    this->alpha_ = 0.618034;
    this->seed_ = 0.5;
    this->z_ = 0;
  }
  [[maybe_unused]] std::vector<scalar_t> get_state() const {
    std::vector<scalar_t> result(2);
    result[0] = this->alpha_;
    result[1] = this->z_;
    return result;
  }
  [[maybe_unused]] void set_state(scalar_t alpha = 0.618034, scalar_t z = 0) {
    this->alpha_ = alpha;
    this->z_ = z;
  }

 private:
  scalar_t alpha_ = 0.618034, seed_ = 0.5, z_ = 0;
};

template <typename scalar_t = float>
struct splitmix {
  explicit splitmix<scalar_t>() : s(12374563468) {}
  scalar_t yield() {
    uint64_t result = (s += 0x9E3779B97f4A7C15);
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    // map to the half-open interval (0, 1): take the top 53 bits and offset by
    // half a ULP, so we never return exactly 0 (would break rnorm's log) or 1
    // (would break floor(gen*max) indexing).
    return static_cast<scalar_t>(
        (static_cast<double>((result ^ (result >> 31)) >> 11) + 0.5) *
        (1.0 / 9007199254740992.0));
  }
  scalar_t operator()() { return this->yield(); }
  uint64_t yield_init() {
    uint64_t result = (s += 0x9E3779B97f4A7C15);
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return result ^ (result >> 31);
  }
  [[maybe_unused]] void set_state(uint64_t seed) { this->s = seed; }
  [[maybe_unused]] std::vector<scalar_t> get_state() const {
    std::vector<scalar_t> result(1);
    result[0] = this->s;
    return result;
  }

 private:
  uint64_t s;
};
template <typename scalar_t = float>
struct xoshiro {
  xoshiro<scalar_t>() {  // NOLINT
    splitmix<scalar_t> gn;
    // seed each 64-bit word from an independent splitmix draw; the previous
    // code derived words via >>32 (half the bits zero) and assigned a (0,1)
    // double to s[2] which truncated to 0, leaving s[2]=s[3]=0 (degenerate)
    s[0] = gn.yield_init();
    s[1] = gn.yield_init();
    s[2] = gn.yield_init();
    s[3] = gn.yield_init();
  }
  scalar_t yield() {
    uint64_t const result = s[0] + s[3];
    uint64_t const t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = bitwise_rotate(s[3], 64, 45);

    // half-open (0, 1), see splitmix::yield
    return static_cast<scalar_t>((static_cast<double>(result >> 11) + 0.5) *
                                 (1.0 / 9007199254740992.0));
  }
  uint64_t bitwise_rotate(uint64_t x, int bits, int rotate_bits) {
    return (x << rotate_bits) | (x >> (bits - rotate_bits));
  }
  scalar_t operator()() { return this->yield(); }
  [[maybe_unused]] void reset() {
    splitmix<scalar_t> gn;
    s[0] = gn.yield_init();
    s[1] = gn.yield_init();
    s[2] = gn.yield_init();
    s[3] = gn.yield_init();
  }
  [[maybe_unused]] void set_state(uint64_t x, uint64_t y, uint64_t z,
                                  uint64_t t) {
    this->s[0] = x;
    this->s[1] = y;
    this->s[2] = z;
    this->s[3] = t;
  }
  [[maybe_unused]] std::vector<scalar_t> get_state() const {
    std::vector<scalar_t> result(4);
    for (size_t i = 0; i < 4; i++) {
      result[i] = this->s[i];
    }
    return result;
  }

 private:
  uint64_t s[4];
};

template <typename scalar_t = float>
struct xorshift {
  xorshift<scalar_t>() {  // NOLINT
    splitmix<scalar_t> gn;
    // independent draws (was x[1] = x[0] >> 32, leaving half the bits zero)
    x[0] = gn.yield_init();
    x[1] = gn.yield_init();
  }
  scalar_t yield() {
    uint64_t t = x[0];
    uint64_t const s = x[1];
    x[0] = s;
    t ^= t << 23;  // a
    t ^= t >> 18;  // b -- Again, the shifts and the multipliers are tunable
    t ^= s ^ (s >> 5);  // c
    x[1] = t;
    // half-open (0, 1), see splitmix::yield
    return static_cast<scalar_t>((static_cast<double>((t + s) >> 11) + 0.5) *
                                 (1.0 / 9007199254740992.0));
  }
  scalar_t operator()() { return this->yield(); }
  [[maybe_unused]] void reset() {
    splitmix<scalar_t> gn;
    x[0] = gn.yield_init();
    x[1] = gn.yield_init();
  }
  [[maybe_unused]] void set_state(uint64_t y, uint64_t z) {
    x[0] = y;
    x[1] = z;
  }
  [[maybe_unused]] std::vector<scalar_t> get_state() const {
    std::vector<scalar_t> result(2);
    for (size_t i = 0; i < 2; i++) {
      result[i] = this->x[i];
    }
    return result;
  }

 private:
  uint64_t x[2]{};
};
}  // namespace nlsolver::rng
namespace nlsolver::finite_difference {
// The 'accuracy' can be 0, 1, 2, 3.
template <typename Callable, typename scalar_t, const size_t accuracy = 0>
void finite_difference_gradient(Callable &f, std::vector<scalar_t> &x,
                                std::vector<scalar_t> &grad) {
  // all constexpr values - this is why accuracy is a template parameter
  // base relative step (~ sqrt(machine epsilon)); scaled per component below
  // Central stencil of truncation order p = 2 (accuracy + 1). Balancing the
  // O(h^p) truncation error against the O(eps / h) rounding error gives a step
  // h ~ eps^(1 / (p + 1)): eps^(1/3) for the two-point stencil, eps^(1/5) for
  // the four-point one, with a gradient error of O(eps^(p / (p + 1))). A
  // sqrt(eps)-sized step, the forward-difference rule, leaves either stencil
  // rounding-bound at O(1e-8) and wastes the extra points.
  constexpr int truncation_order = 2 * (accuracy + 1);
  const scalar_t eps =
      std::pow(std::numeric_limits<scalar_t>::epsilon(),
               static_cast<scalar_t>(1) / (truncation_order + 1));
  constexpr std::array<scalar_t, 20> coeff = {1,   -1,   1,   -8,   8,  -1, -1,
                                              9,   -45,  45,  -9,   1,  3,  -32,
                                              168, -672, 672, -168, 32, -3};
  constexpr std::array<scalar_t, 20> coeff2 = {
      1, -1, -2, -1, 1, 2, -3, -2, -1, 1, 2, 3, -4, -3, -2, -1, 1, 2, 3, 4};
  constexpr std::array<scalar_t, 4> dd = {2, 12, 60, 840};
  constexpr int innerSteps = 2 * (accuracy + 1);
  constexpr std::array<size_t, 4> offset_index = {0, 2, 6, 12};
  constexpr size_t offset = offset_index[accuracy];
  // actual stuff that should exist at runtime
  const size_t x_size = x.size();
  std::fill(grad.begin(), grad.end(), 0.0);
  for (size_t d = 0; d < x_size; d++) {
    const scalar_t h = eps * std::max(static_cast<scalar_t>(1), std::abs(x[d]));
    for (size_t s = 0; s < innerSteps; ++s) {
      scalar_t tmp = x[d];
      x[d] += coeff2[offset + s] * h;
      grad[d] += coeff[offset + s] * f(x);
      x[d] = tmp;
    }
    grad[d] /= (dd[accuracy] * h);
  }
}

template <typename Callable, typename scalar_t, const size_t accuracy = 0>
void finite_difference_hessian(Callable &f, std::vector<scalar_t> &x,  // NOLINT
                               std::vector<scalar_t> &hess) {
  constexpr auto eps_mult = static_cast<scalar_t>(1.0 / 4.0);
  const scalar_t eps =
      std::pow(std::numeric_limits<scalar_t>::epsilon(), eps_mult);
  // std::cout << "Epsilon: "<< eps;
  const size_t p = x.size();
  if constexpr (accuracy == 0) {
    // symmetric central differences (the previous one-sided forward stencil was
    // only first-order accurate and biased). Per-dimension steps scaled by |x|.
    const scalar_t f0 = f(x);
    for (size_t i = 0; i < p; i++) {
      const scalar_t temp_i = x[i];
      const scalar_t hi =
          eps * std::max(static_cast<scalar_t>(1), std::abs(temp_i));
      for (size_t j = i; j < p; j++) {
        const scalar_t temp_j = x[j];
        scalar_t value;
        if (i == j) {
          // diagonal: (f(x+h) - 2 f(x) + f(x-h)) / h^2
          x[i] = temp_i + hi;
          const scalar_t fp = f(x);
          x[i] = temp_i - hi;
          const scalar_t fm = f(x);
          x[i] = temp_i;
          value = (fp - 2.0 * f0 + fm) / (hi * hi);
        } else {
          // mixed: (f++ - f+- - f-+ + f--) / (4 hi hj)
          const scalar_t hj =
              eps * std::max(static_cast<scalar_t>(1), std::abs(temp_j));
          x[i] = temp_i + hi;
          x[j] = temp_j + hj;
          const scalar_t fpp = f(x);
          x[j] = temp_j - hj;
          const scalar_t fpm = f(x);
          x[i] = temp_i - hi;
          const scalar_t fmm = f(x);
          x[j] = temp_j + hj;
          const scalar_t fmp = f(x);
          x[i] = temp_i;
          x[j] = temp_j;
          value = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
        }
        hess[i * p + j] = value;
        hess[j * p + i] = value;
      }
    }
  } else {
    const scalar_t denom = (600.0 * eps * eps), two_eps = 2 * eps,
                   three_eps = 3 * eps, four_eps = 4 * eps;
    for (size_t i = 0; i < p; i++) {
      const scalar_t temp_i = x[i];
      for (size_t j = 0; j < p; j++) {
        scalar_t result = 0.0, temp = 0.0;
        scalar_t temp_j = x[j];
        x[i] += eps;      // x_i + eps
        x[j] -= two_eps;  // x_j - 2 * eps
        temp += f(x);
        x[i] += eps;  // x_i + 2 * eps
        x[j] += eps;  // x_j - eps
        temp += f(x);
        x[i] -= four_eps;  // x_i - 2 * eps
        x[j] += two_eps;   // x_j + eps
        temp += f(x);
        x[i] += eps;  // x_i - eps
        x[j] += eps;  // x_j + 2 * eps
        temp += f(x);
        result -= 63 * temp;
        temp = 0.0;
        // x_i remains at (x_i - eps)
        x[j] -= four_eps;  // x_j - 2 * eps
        temp += f(x);
        x[i] -= eps;  // x_i - 2 * eps
        x[j] += eps;  // x_j - eps
        temp += f(x);
        x[i] += three_eps;  // x_i + eps
        x[j] += three_eps;  // x_j + 2 * eps
        temp += f(x);
        x[i] += eps;  // x_i + 2 * eps
        x[j] -= eps;  // x_j + eps
        temp += f(x);
        result += 63 * temp;
        temp = 0.0;
        // x_i remains at (x_i + 2 * eps)
        x[j] -= three_eps;  // x_j -2 * eps
        temp += f(x);
        x[i] -= four_eps;  // x_i - 2 * eps
        x[j] += four_eps;  // x_j + 2 * eps
        temp += f(x);
        // x_i remains at (x_i - 2 * eps)
        x[j] -= four_eps;  // x_j - 2 * eps
        temp -= f(x);
        x[i] += four_eps;  // x_i + 2 * eps
        x[j] += four_eps;  // x_j + 2 * eps
        temp -= f(x);
        result += 44 * temp;
        temp = 0.0;
        x[i] -= three_eps;  // x_i - eps
        x[j] -= three_eps;  // x_j - eps
        temp += f(x);
        x[i] += two_eps;  // x_i + eps
        x[j] += two_eps;  // x_j + eps
        temp += f(x);
        // x_i remains at (x_i + eps)
        x[j] -= two_eps;  // x_j - eps
        temp -= f(x);
        x[i] -= two_eps;  // x_i - eps
        x[j] += two_eps;  // x_j + eps
        temp -= f(x);
        result += 74 * temp;
        // reset
        x[i] = temp_i;
        x[j] = temp_j;
        hess[i * p + j] = result / denom;
      }
    }
  }
}
}  // namespace nlsolver::finite_difference
namespace nlsolver::linesearch {
template <typename scalar_t>
scalar_t max_abs(scalar_t x, scalar_t y, scalar_t z) {
  return std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
}
// this is in a bit of a sad state since the result gets discarded
// and everything is done by reference - potentially might also
// benefit from branch elimination
template <typename scalar_t>
[[maybe_unused]] static int cstep(scalar_t &stx, scalar_t &fx,
                                  scalar_t &dx,   // NOLINT
                                  scalar_t &sty,  // NOLINT
                                  scalar_t &fy, scalar_t &dy,
                                  scalar_t &stp,  // NOLINT
                                  scalar_t &fp,   // NOLINT
                                  scalar_t &dp, bool &brackt,
                                  scalar_t &stpmin,               // NOLINT
                                  scalar_t &stpmax, int &info) {  // NOLINT
  info = 0;
  bool bound;

  // Check the input parameters for errors.
  if ((brackt & ((stp <= std::min<scalar_t>(stx, sty)) ||
                 (stp >= std::max<scalar_t>(stx, sty)))) ||
      (dx * (stp - stx) >= 0.0) || (stpmax < stpmin)) {
    return -1;
  }

  scalar_t sgnd = dp * (dx / fabs(dx));
  scalar_t stpf = 0;
  scalar_t stpc;
  scalar_t stpq;

  if (fp > fx) {
    info = 1;
    bound = true;
    scalar_t theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
    scalar_t s = max_abs(theta, dx, dp);
    scalar_t gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp < stx) gamma = -gamma;
    scalar_t p = (gamma - dx) + theta;
    scalar_t q = ((gamma - dx) + gamma) + dp;
    scalar_t r = p / q;
    stpc = stx + r * (stp - stx);
    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
    if (fabs(stpc - stx) < fabs(stpq - stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq - stpc) / 2;
    brackt = true;
  } else if (sgnd < 0.0) {
    info = 2;
    bound = false;
    scalar_t theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    scalar_t s = max_abs(theta, dx, dp);
    scalar_t gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp > stx) gamma = -gamma;

    scalar_t p = (gamma - dp) + theta;
    scalar_t q = ((gamma - dp) + gamma) + dx;
    scalar_t r = p / q;
    stpc = stp + r * (stx - stp);
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (fabs(stpc - stp) > fabs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    brackt = true;
  } else if (fabs(dp) < fabs(dx)) {
    info = 3;
    bound = true;
    scalar_t theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    scalar_t s = max_abs(theta, dx, dp);
    scalar_t gamma = s * sqrt(std::max<scalar_t>(
                             static_cast<scalar_t>(0.),
                             (theta / s) * (theta / s) - (dx / s) * (dp / s)));
    if (stp > stx) gamma = -gamma;
    scalar_t p = (gamma - dp) + theta;
    scalar_t q = (gamma + (dx - dp)) + gamma;
    scalar_t r = p / q;
    if ((r < 0.0) & (gamma != 0.0)) {
      stpc = stp + r * (stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (brackt) {
      if (fabs(stp - stpc) < fabs(stp - stpq)) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
    } else {
      if (fabs(stp - stpc) > fabs(stp - stpq)) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
    }
  } else {
    info = 4;
    bound = false;
    if (brackt) {
      scalar_t theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
      scalar_t s = max_abs(theta, dy, dp);
      scalar_t gamma =
          s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
      if (stp > sty) gamma = -gamma;

      scalar_t p = (gamma - dp) + theta;
      scalar_t q = ((gamma - dp) + gamma) + dy;
      scalar_t r = p / q;
      stpc = stp + r * (sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }

    stx = stp;
    fx = fp;
    dx = dp;
  }

  stpf = std::clamp(stpf, stpmin, stpmax);
  stp = stpf;

  if (brackt & bound) {
    if (sty > stx) {
      stp = std::min<scalar_t>(stx + static_cast<scalar_t>(0.66) * (sty - stx),
                               stp);
    } else {
      stp = std::max<scalar_t>(stx + static_cast<scalar_t>(0.66) * (sty - stx),
                               stp);
    }
  }
  return 0;
}

template <typename Callable, typename Grad, typename scalar_t = double>
[[maybe_unused]] static int cvsrch(
    Callable &f, std::vector<scalar_t> &x, scalar_t current_f_value,
    std::vector<scalar_t> &gradient, scalar_t *stp,
    const std::vector<scalar_t> &search_direction,
    std::vector<scalar_t> &linesearch_temp, Grad &g) {
  // we rewrite this from MIN-LAPACK and some MATLAB code
  int info = 0;
  int infoc = 1;
  constexpr scalar_t xtol = 1e-15;
  constexpr scalar_t ftol = 1e-4;
  constexpr scalar_t gtol = 1e-2;
  constexpr scalar_t stpmin = 1e-15;
  constexpr scalar_t stpmax = 1e15;
  constexpr scalar_t xtrapf = 4;
  constexpr int maxfev = 20;
  int nfev = 0;

  scalar_t dginit =
      nlsolver::math::dot(gradient.data(), search_direction.data(), x.size());
  if (dginit >= 0.0) {
    return -1;
  }

  bool brackt = false;
  bool stage1 = true;

  scalar_t finit = current_f_value;
  scalar_t dgtest = ftol * dginit;
  scalar_t width = stpmax - stpmin;
  scalar_t width1 = 2 * width;
  // vector_t wa = x->eval();

  scalar_t stx = 0.0;
  scalar_t fx = finit;
  scalar_t dgx = dginit;
  scalar_t sty = 0.0;
  scalar_t fy = finit;
  scalar_t dgy = dginit;

  scalar_t stmin;
  scalar_t stmax;

  while (true) {
    // Make sure we stay in the interval when setting min/max-step-width.
    if (brackt) {
      stmin = std::min<scalar_t>(stx, sty);
      stmax = std::max<scalar_t>(stx, sty);
    } else {
      stmin = stx;
      stmax = *stp + xtrapf * (*stp - stx);
    }

    // Force the step to be within the bounds stpmax and stpmin.
    *stp = std::clamp(*stp, stpmin, stpmax);

    // Oops, let us return the last reliable values.
    if ((brackt && ((*stp <= stmin) || (*stp >= stmax))) ||
        (nfev >= maxfev - 1) || (infoc == 0) ||
        (brackt && ((stmax - stmin) <= (xtol * stmax)))) {
      *stp = stx;
    }

    // Test new point.
    for (size_t i = 0; i < x.size(); i++) {
      linesearch_temp[i] = x[i] + *stp * search_direction[i];
    }
    current_f_value = f(linesearch_temp);
    g(linesearch_temp, gradient);
    nfev++;
    scalar_t dg =
        nlsolver::math::dot(gradient.data(), search_direction.data(), x.size());
    scalar_t ftest1 = finit + *stp * dgtest;

    // All possible convergence tests.
    if ((brackt & ((*stp <= stmin) | (*stp >= stmax))) | (infoc == 0)) info = 6;
    if ((*stp == stpmax) & (current_f_value <= ftest1) & (dg <= dgtest))
      info = 5;
    if ((*stp == stpmin) & ((current_f_value > ftest1) | (dg >= dgtest)))
      info = 4;
    if (nfev >= maxfev) info = 3;
    if (brackt & (stmax - stmin <= xtol * stmax)) info = 2;
    if ((current_f_value <= ftest1) & (fabs(dg) <= gtol * (-dginit))) info = 1;
    // Terminate when convergence reached.
    if (info != 0) return -1;

    if (stage1 & (current_f_value <= ftest1) &
        (dg >= std::min<scalar_t>(ftol, gtol) * dginit))
      stage1 = false;

    if (stage1 & (current_f_value <= fx) & (current_f_value > ftest1)) {
      scalar_t fm = current_f_value - *stp * dgtest;
      scalar_t fxm = fx - stx * dgtest;
      scalar_t fym = fy - sty * dgtest;
      scalar_t dgm = dg - dgtest;
      scalar_t dgxm = dgx - dgtest;
      scalar_t dgym = dgy - dgtest;

      cstep(stx, fxm, dgxm, sty, fym, dgym, *stp, fm, dgm, brackt, stmin, stmax,
            infoc);

      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    } else {
      // This is ugly and some variables should be moved to the class scope.
      cstep(stx, fx, dgx, sty, fy, dgy, *stp, current_f_value, dg, brackt,
            stmin, stmax, infoc);
    }

    if (brackt) {
      if (fabs(sty - stx) >= 0.66 * width1) {
        *stp = stx + 0.5 * (sty - stx);
      }
      width1 = width;
      width = fabs(sty - stx);
    }
  }
  return 0;
}
template <typename Callable, typename scalar_t>
[[nodiscard]] inline scalar_t line_at_alpha(
    Callable &f, std::vector<scalar_t> &linesearch_temp,
    std::vector<scalar_t> &x, const std::vector<scalar_t> &search_direction,
    const scalar_t alpha = 1.0) {
  const size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    linesearch_temp[i] = x[i] + alpha * search_direction[i];
  }
  return f(linesearch_temp);
}
template <typename Callable, typename scalar_t = double>
[[nodiscard]] [[maybe_unused]] scalar_t armijo_search(
    Callable &f, const scalar_t current_f_val, std::vector<scalar_t> &x,
    std::vector<scalar_t> &gradient,
    const std::vector<scalar_t> &search_direction, scalar_t alpha = 1.0) {
  constexpr scalar_t c = 0.2, rho = 0.9;
  const size_t x_size = x.size();
  scalar_t limit = nlsolver::math::dot(gradient.data(), search_direction.data(),
                                       static_cast<int>(x_size)) *
                   c;

  std::vector<scalar_t> linesearch_temp(x_size);
  scalar_t search_val = line_at_alpha<Callable, scalar_t>(
      f, linesearch_temp, x, search_direction, alpha);
  while (search_val > (current_f_val + alpha * limit)) {
    alpha *= rho;
    search_val = line_at_alpha<Callable, scalar_t>(f, linesearch_temp, x,
                                                   search_direction, alpha);
  }
  return alpha;
}
// this is just an overload which allows us to pass a temporary
template <typename Callable, typename scalar_t = double>
[[maybe_unused]] scalar_t armijo_search(
    Callable &f, const scalar_t current_f_val, std::vector<scalar_t> &x,
    std::vector<scalar_t> &gradient,
    const std::vector<scalar_t> &search_direction,
    std::vector<scalar_t> &linesearch_temp, scalar_t alpha = 1.0) {
  constexpr scalar_t c = 0.2, rho = 0.9;
  const size_t x_size = x.size();
  scalar_t limit = nlsolver::math::dot(gradient.data(), search_direction.data(),
                                       static_cast<int>(x_size)) *
                   c;
  scalar_t search_val =
      line_at_alpha(f, linesearch_temp, x, search_direction, alpha);
  while (search_val > (current_f_val + alpha * limit)) {
    alpha *= rho;
    search_val = line_at_alpha(f, linesearch_temp, x, search_direction, alpha);
  }
  return alpha;
}
// this is just an overload which allows us to pass a temporary
template <typename Callable, typename scalar_t = double>
[[maybe_unused]] scalar_t armijo_search(
    Callable &f, std::vector<scalar_t> &x, std::vector<scalar_t> &gradient,
    const std::vector<scalar_t> &search_direction,
    std::vector<scalar_t> &linesearch_temp, scalar_t alpha = 1.0) {
  constexpr scalar_t c = 0.2, rho = 0.9;
  const scalar_t current_f_val = f(x);
  const size_t x_size = x.size();
  scalar_t limit = nlsolver::math::dot(gradient.data(), search_direction.data(),
                                       static_cast<int>(x_size)) *
                   c;
  scalar_t search_val =
      line_at_alpha(f, linesearch_temp, x, search_direction, alpha);
  while (search_val > (current_f_val + alpha * limit)) {
    alpha *= rho;
    search_val = line_at_alpha(f, linesearch_temp, x, search_direction, alpha);
  }
  return alpha;
}

//
template <typename Callable, typename Grad, typename scalar_t = double>
scalar_t more_thuente_search(Callable &f, const scalar_t current_f_val,
                             std::vector<scalar_t> &x,
                             std::vector<scalar_t> &gradient,
                             const std::vector<scalar_t> &search_direction,
                             std::vector<scalar_t> &linesearch_temp,
                             scalar_t alpha, Grad g) {
  scalar_t alpha_ = alpha;
  cvsrch(f, x, current_f_val, gradient, &alpha_, search_direction,
         linesearch_temp, g);
  return alpha_;
}
//
template <typename Callable, typename Grad, typename scalar_t = double>
[[maybe_unused]] scalar_t more_thuente_search(
    Callable &f, std::vector<scalar_t> &x, std::vector<scalar_t> &gradient,
    const std::vector<scalar_t> &search_direction,
    std::vector<scalar_t> &linesearch_temp, scalar_t alpha, Grad g) {
  const scalar_t current_f_val = f(x);
  scalar_t alpha_ = alpha;
  cvsrch(f, x, current_f_val, gradient, &alpha_, search_direction,
         linesearch_temp, g);
  return alpha_;
}
}  // namespace nlsolver::linesearch
namespace nlsolver {
template <typename scalar_t = double>
inline scalar_t max_abs_vec(const std::vector<scalar_t> &x) {
  auto result = std::abs(x[0]);
  for (size_t i = 1; i < x.size(); i++) {
    scalar_t temp = std::abs(x[i]);
    if (result < temp) {
      result = temp;
    }
  }
  return result;
}
template <typename scalar_t = double>
struct simplex {
  explicit simplex<scalar_t>(const size_t i = 0) {
    this->vals = std::vector<std::vector<scalar_t>>(i + 1);
  }
  explicit simplex<scalar_t>(const std::vector<scalar_t> &x,
                             const scalar_t step = -1) {
    std::vector<std::vector<scalar_t>> init_simplex(x.size() + 1);
    // init_simplex[0] = x;
    //  this follows Gao and Han, see:
    //  'Proper initialization is crucial for the Nelder–Mead simplex search.'
    //  (2019), Wessing, S.  Optimization Letters 13, p. 847–856
    //  (also at https://link.springer.com/article/10.1007/s11590-018-1284-4)
    //  default initialization
    if (step < 0) {
      // get infinity norm of initial vector
      scalar_t x_inf_norm = max_abs_vec(x);
      // if smaller than 1, set to 1
      scalar_t a = x_inf_norm < 1.0 ? 1.0 : x_inf_norm;
      // if larger than 10, set to 10
      scalar_t scale = a < 10 ? a : 10;
      for (auto &vertex : init_simplex) {
        vertex = x;
      }
      for (size_t i = 1; i < init_simplex.size(); i++) {
        init_simplex[i][i - 1] += scale;
      }
      // update first simplex point
      auto n = static_cast<scalar_t>(x.size());
      for (size_t i = 0; i < x.size(); i++) {
        init_simplex[0][i] = x[i] + ((1.0 - sqrt(n + 1.0)) / n * scale);
      }
      // otherwise, first element of simplex has unchanged starting values
    } else {
      for (auto &vertex : init_simplex) {
        vertex = x;
      }
      for (size_t i = 1; i < init_simplex.size(); i++) {
        init_simplex[i][i - 1] += step;
      }
    }
    this->vals = init_simplex;
  }
  void replace(std::vector<scalar_t> &new_val, const size_t at) {
    this->vals[at] = new_val;
  }
  [[maybe_unused]] void replace(std::vector<scalar_t> &new_val, const size_t at,
                                const std::vector<scalar_t> &upper,
                                const std::vector<scalar_t> &lower,
                                const scalar_t inversion_eps = 0.00001) {
    for (size_t i = 0; i < new_val.size(); i++) {
      this->vals[at][i] = new_val[i] < lower[i]   ? lower[i] + inversion_eps
                          : new_val[i] > upper[i] ? upper[i] - inversion_eps
                                                  : new_val[i];
    }
  }
  [[nodiscard]] size_t size() const { return this->vals.size(); }
  std::vector<std::vector<scalar_t>> vals;
};

template <typename scalar_t = double>
inline void update_centroid(std::vector<scalar_t> &centroid,
                            const simplex<scalar_t> &x, const size_t except) {
  // reset centroid - fill with 0
  std::fill(centroid.begin(), centroid.end(), 0.0);
  size_t i = 0;
  for (; i < except; i++) {
    // TODO(JSzitas): SIMD Candidate
    for (size_t j = 0; j < centroid.size(); j++) {
      centroid[j] += x.vals[i][j];
    }
  }
  i = except + 1;
  for (; i < x.size(); i++) {
    for (size_t j = 0; j < centroid.size(); j++) {
      centroid[j] += x.vals[i][j];
    }
  }
  for (auto &val : centroid) val /= static_cast<scalar_t>(i - 1);
}
// bound version
template <typename scalar_t = double, const bool reflect = false,
          const bool bound = false>
inline void simplex_transform(const std::vector<scalar_t> &point,
                              const std::vector<scalar_t> &centroid,
                              std::vector<scalar_t> &result,
                              const scalar_t coef,
                              const std::vector<scalar_t> &upper,
                              const std::vector<scalar_t> &lower) {
  for (size_t i = 0; i < point.size(); i++) {
    scalar_t temp = 0.0;
    // TODO(JSzitas): SIMD Candidate
    if constexpr (reflect) {
      temp = centroid[i] + coef * (centroid[i] - point[i]);
    } else {
      temp = centroid[i] + coef * (point[i] - centroid[i]);
    }
    if constexpr (bound) {
      temp = std::clamp(temp, lower[i], upper[i]);
    }
    result[i] = temp;
  }
}

template <typename scalar_t = double>
inline void shrink(simplex<scalar_t> &current_simplex, const size_t best,
                   const scalar_t sigma) {
  // take a reference to the best vector
  std::vector<scalar_t> &best_val = current_simplex.vals[best];
  const size_t n = best_val.size();
  for (size_t i = 0; i < best; i++) {
    // update all items in current vector using the best vector -
    // hopefully the contiguous data here can help a bit with cache
    // locality
    for (size_t j = 0; j < n; j++) {
      // TODO(JSzitas): SIMD Candidate
      current_simplex.vals[i][j] =
          best_val[j] + sigma * (current_simplex.vals[i][j] - best_val[j]);
    }
  }
  // skip the best point - this uses separate loops, so we do not have to do
  // extra work (e.g. check i == best) which could lead to a branch
  // mis-prediction
  for (size_t i = best + 1; i < current_simplex.size(); i++) {
    for (size_t j = 0; j < n; j++) {
      // TODO(JSzitas): SIMD Candidate
      current_simplex.vals[i][j] =
          best_val[j] + sigma * (current_simplex.vals[i][j] - best_val[j]);
    }
  }
}

template <typename scalar_t = double>
static inline scalar_t std_err(const std::vector<scalar_t> &x) {
  // fewer than two samples have no spread; avoid /0 (size 1 -> /(i-1)=/0, and
  // size 0 -> /i=/0) which would otherwise return NaN.
  if (x.size() < 2) return 0;
  size_t i = 0;
  scalar_t mean_val = 0, result = 0;
  // TODO(JSzitas): SIMD Candidate
  for (; i < x.size(); i++) {
    mean_val += x[i];
  }
  mean_val /= static_cast<scalar_t>(i);
  i = 0;
  for (; i < x.size(); i++) {
    result += pow(x[i] - mean_val, 2);
  }
  result /= static_cast<scalar_t>(i - 1);
  return sqrt(result);
}

template <typename scalar_t = double>
struct solver_status {
  solver_status<scalar_t>(const scalar_t f_val, const size_t iter_used,
                          const size_t f_calls_used,
                          const size_t grad_evals_used = 0ul,
                          const size_t hess_evals_used = 0ul,
                          const bool converged = true)
      : f_value(f_val),
        iteration(iter_used),
        function_calls_used(f_calls_used),
        gradient_evals_used(grad_evals_used),
        hessian_evals_used(hess_evals_used),
        converged_(converged) {}
  // whether the solver reached a valid result (e.g. a root finder returns false
  // here when the initial interval did not bracket a root)
  [[nodiscard]] bool success() const { return converged_; }
  void print() const {
    if (!converged_) {
      std::cout << "Solver did NOT converge (result is not valid)" << std::endl;
    }
    std::cout << "Function calls used: " << this->function_calls_used
              << std::endl;
    std::cout << "Algorithm iterations used: " << this->iteration << std::endl;
    if (gradient_evals_used > 0) {
      std::cout << "Gradient evaluations used: " << this->gradient_evals_used
                << std::endl;
    }
    if (hessian_evals_used > 0) {
      std::cout << "Hessian evaluations used: " << this->hessian_evals_used
                << std::endl;
    }
    std::cout << "With final function value of " << this->f_value << std::endl;
  }
  std::tuple<size_t, size_t, scalar_t, size_t, size_t> get_summary() const {
    return std::make_tuple(this->function_calls_used, this->iteration,
                           this->f_value, this->gradient_evals_used,
                           this->hessian_evals_used);
  }
  void add(const solver_status<scalar_t> &additional_runs) {
    auto other = additional_runs.get_summary();
    this->function_calls_used += std::get<0>(other);
    this->iteration += std::get<1>(other);
    this->f_value = std::get<2>(other);
    this->gradient_evals_used += std::get<3>(other);
    this->hessian_evals_used += std::get<4>(other);
  }

 private:
  scalar_t f_value;
  size_t iteration, function_calls_used, gradient_evals_used,
      hessian_evals_used;
  bool converged_ = true;
};

template <typename Callable, typename scalar_t = double>
class NelderMead {
 private:
  Callable &f;
  const scalar_t step, alpha, gamma, rho, sigma;
  scalar_t eps;
  std::vector<scalar_t> point_values;
  const size_t max_iter, no_change_best_tol, restarts;

 public:
  // constructor
  explicit NelderMead<Callable, scalar_t>(
      Callable &f, const scalar_t step = -1, const scalar_t alpha = 1,
      const scalar_t gamma = 2, const scalar_t rho = 0.5,
      const scalar_t sigma = 0.5, const scalar_t eps = 1e-6,
      const size_t max_iter = 500, const size_t no_change_best_tol = 20,
      const size_t restarts = 0)
      : f(f),
        step(step),
        alpha(alpha),
        gamma(gamma),
        rho(rho),
        sigma(sigma),
        eps(eps),
        max_iter(max_iter),
        no_change_best_tol(no_change_best_tol),
        restarts(restarts) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> upper, lower;
    auto res = this->solve<true, false>(x, upper, lower);
    for (size_t i = 0; i < this->restarts; i++) {
      res.add(this->solve<true, false>(x, upper, lower));
    }
    return res;
  }
  // minimize with known bounds interface; bounds are passed as (lower, upper)
  // for consistency with the other solvers
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    auto res = this->solve<true, true>(x, upper, lower);
    for (size_t i = 0; i < this->restarts; i++) {
      res.add(this->solve<true, true>(x, upper, lower));
    }
    return res;
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> upper, lower;
    auto res = this->solve<false, false>(x, upper, lower);
    for (size_t i = 0; i < this->restarts; i++) {
      res.add(this->solve<false, false>(x, upper, lower));
    }
    return res;
  }
  // maximize with known bounds interface; bounds are passed as (lower, upper)
  // for consistency with the other solvers
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    auto res = this->solve<false, true>(x, upper, lower);
    for (size_t i = 0; i < this->restarts; i++) {
      res.add(this->solve<false, true>(x, upper, lower));
    }
    return res;
  }

 private:
  template <const bool minimize = true, const bool bound = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &upper,
                                const std::vector<scalar_t> &lower) {
    // set up simplex
    simplex<scalar_t> current_simplex(x, this->step);
    std::vector<scalar_t> scores(current_simplex.size());
    /* this basically ensures that for minimization we are seeking
     * minimum of function **f**, and for maximization we are seeking minimum of
     * **-f** */
    size_t function_calls_used = 0;
    auto f_lam = [&](decltype(x) &x_) {
      constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
      auto val = f_multiplier * f(x_);
      function_calls_used++;
      return val;
    };
    // score simplex values
    for (size_t i = 0; i < current_simplex.size(); i++) {
      scores[i] = f_lam(current_simplex.vals[i]);
    }
    // relative-plus-absolute convergence tolerance, computed once from the
    // magnitude of the initial simplex. (Previously this overwrote the member
    // eps, which corrupted the tolerance and compounded it across restarts.)
    const scalar_t conv_tol = this->eps * (std::abs(scores[0]) + 1.0);
    // best / worst / second-worst vertices and the stall counter
    size_t best = 0, worst = 0, second_worst = 0, prev_worst = 0,
           no_change_iter = 0;
    scalar_t last_best_val = std::numeric_limits<scalar_t>::infinity();

    size_t iter = 0;
    std::vector<scalar_t> centroid(x.size()), temp_reflect(x.size()),
        temp_expand(x.size()), temp_contract(x.size());
    scalar_t ref_score, exp_score, cont_score, fun_std_err;
    bool we_called_shrink = false;
    // simplex iteration
    while (true) {
      // find best (argmin), worst (argmax) and second-worst vertices with a
      // correct scan (the previous chained if/else left second_worst stale)
      best = 0;
      prev_worst =
          worst;  // used to decide whether the centroid needs a refresh
      worst = 0;
      fun_std_err = std_err(scores);
      for (size_t i = 1; i < scores.size(); i++) {
        if (scores[i] < scores[best]) best = i;
        if (scores[i] > scores[worst]) worst = i;
      }
      second_worst = (worst == 0) ? 1 : 0;
      for (size_t i = 0; i < scores.size(); i++) {
        if (i != worst && scores[i] > scores[second_worst]) second_worst = i;
      }
      // stall detection on the best *value*. Nelder-Mead only ever replaces the
      // worst vertex, so the best *index* can stay fixed while the best value
      // keeps improving - tracking the index caused premature termination.
      if (last_best_val - scores[best] > conv_tol) {
        no_change_iter = 0;  // meaningful improvement this iteration
      } else {
        no_change_iter++;
      }
      if (scores[best] < last_best_val) last_best_val = scores[best];
      // check whether we should stop - either by exceeding iterations or by
      // reaching tolerance
      if (iter >= this->max_iter || fun_std_err < conv_tol ||
          no_change_iter >= this->no_change_best_tol) {
        x = current_simplex.vals[best];
        return solver_status<scalar_t>(scores[best], iter, function_calls_used);
      }
      iter++;
      // compute centroid of all points except for the worst one
      if (prev_worst != worst || we_called_shrink) {
        update_centroid(centroid, current_simplex, worst);
        we_called_shrink = false;
      }
      // reflect worst point
      simplex_transform<scalar_t, true, bound>(current_simplex.vals[worst],
                                               centroid, temp_reflect,
                                               this->alpha, upper, lower);
      // score reflected point
      ref_score = f_lam(temp_reflect);
      // if reflected point is better than second worst, not better than best
      if (ref_score >= scores[best] && ref_score < scores[second_worst]) {
        current_simplex.replace(temp_reflect, worst);
        scores[worst] = ref_score;
        // otherwise if this is the best score so far, expand
      } else if (ref_score < scores[best]) {
        simplex_transform<scalar_t, false, bound>(
            temp_reflect, centroid, temp_expand, this->gamma, upper, lower);
        // obtain score for expanded point
        exp_score = f_lam(temp_expand);
        // if this is better than the expanded point score, replace worst point
        // with the expanded point, otherwise replace it with reflected point
        current_simplex.replace(
            exp_score < ref_score ? temp_expand : temp_reflect, worst);
        scores[worst] = exp_score < ref_score ? exp_score : ref_score;
        // otherwise we have a point  worse than the 'second worst'
      } else {
        // contract (outside if the reflected point beat the worst, inside
        // otherwise) - moves toward the centroid, hence reflect=false
        simplex_transform<scalar_t, false, bound>(
            ref_score < scores[worst]
                ? temp_reflect
                :
                // or point is the worst point so far - contract inside
                current_simplex.vals[worst],
            centroid, temp_contract, this->rho, upper, lower);
        cont_score = f_lam(temp_contract);
        // if this contraction is better than the reflected point or worst
        if (cont_score <
            (ref_score < scores[worst] ? ref_score : scores[worst])) {
          // replace worst point with contracted point
          current_simplex.replace(temp_contract, worst);
          scores[worst] = cont_score;
          // otherwise shrink
        } else {
          // if we had not violated the bounds before shrinking, shrinking
          // will not cause new violations - hence no bounds applied here
          shrink(current_simplex, best, this->sigma);
          // only in this case do we have to score again
          for (size_t i = 0; i < best; i++) {
            scores[i] = f_lam(current_simplex.vals[i]);
          }
          // we have not updated the best value - hence no need to 'rescore'
          for (size_t i = best + 1; i < current_simplex.size(); i++) {
            scores[i] = f_lam(current_simplex.vals[i]);
          }
          we_called_shrink = true;
        }
      }
    }
  }
};

template <typename RNG, typename scalar_t = double>
static inline std::vector<scalar_t> generate_sequence(
    const std::vector<scalar_t> &offset, RNG &generator) {
  const size_t samples = offset.size();
  std::vector<scalar_t> result(samples);
  for (size_t i = 0; i < samples; i++) {
    // centre the population around the supplied point with a spread that has a
    // floor of 1, so a zero starting coordinate does not freeze that dimension
    // (the old `(gen-0.5)*offset[i]` centred at 0 and collapsed when offset==0)
    result[i] =
        offset[i] + (generator() - 0.5) * 2.0 * (std::abs(offset[i]) + 1.0);
  }
  return result;
}

template <typename RNG, typename scalar_t = double>
static inline std::vector<std::vector<scalar_t>> init_agents(
    const std::vector<scalar_t> &init, RNG &generator, const size_t n_agents) {
  std::vector<std::vector<scalar_t>> agents(n_agents);
  // first element of simplex is unchanged starting values
  for (auto &agent : agents) {
    agent = generate_sequence(init, generator);
  }
  return agents;
}
// static inline
template <typename RNG>
size_t generate_index(const size_t max, RNG &generator) {
  // a slightly evil typecast
  return static_cast<size_t>(generator() * max);
}

template <typename RNG>
static inline std::array<size_t, 4> generate_indices(const size_t fixed,
                                                     const size_t max,
                                                     RNG &generator) {
  // pick 4 distinct indices, the first being the reference agent `fixed`.
  // Only 4 values are ever compared, so a linear distinctness scan is far
  // cheaper than a std::unordered_set (no hashing, no per-call allocation -
  // this was a measurable hotspot). Requires max >= 4 to terminate.
  std::array<size_t, 4> result;  // NOLINT
  result[0] = fixed;
  size_t samples = 1;
  while (samples < 4) {
    const size_t proposal = generate_index(max, generator);
    bool distinct = true;
    for (size_t k = 0; k < samples; k++) {
      if (result[k] == proposal) {
        distinct = false;
        break;
      }
    }
    if (distinct) result[samples++] = proposal;
  }
  return result;
}

template <typename RNG, typename scalar_t = double>
static inline void propose_new_agent(
    const std::array<size_t, 4> &ids, std::vector<scalar_t> &proposal,
    const std::vector<std::vector<scalar_t>> &agents,
    const scalar_t crossover_probability, const scalar_t diff_weight,
    RNG &generator) {
  // pick dimensionality to always change
  size_t dim = generate_index(proposal.size(), generator);
  for (size_t i = 0; i < proposal.size(); i++) {
    // check if we mutate
    if (generator() < crossover_probability || i == dim) {
      proposal[i] = agents[ids[1]][i] +
                    diff_weight * (agents[ids[2]][i] - agents[ids[3]][i]);
    } else {
      // no replacement
      proposal[i] = agents[ids[0]][i];
    }
  }
}

enum RecombinationStrategy { best, random };

template <typename Callable, typename RNG, typename scalar_t = double,
          RecombinationStrategy RecombinationType = random>
class DE {
 private:
  Callable &f;
  RNG &generator;
  const scalar_t crossover_prob, differential_weight, eps;
  const size_t pop_size, max_iter, best_value_no_change;

 public:
  // constructor
  [[maybe_unused]] DE<Callable, RNG, scalar_t, RecombinationType>(
      Callable &f, RNG &generator, const scalar_t crossover_prob = 0.9,
      const scalar_t differential_weight = 0.8, const scalar_t eps = 1e-3,
      const size_t pop_size = 50, const size_t max_iter = 1000,
      const size_t best_val_no_change = 50)
      : f(f),
        generator(generator),
        crossover_prob(crossover_prob),
        differential_weight(differential_weight),
        eps(eps),
        pop_size(pop_size),
        max_iter(max_iter),
        best_value_no_change(best_val_no_change) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    const size_t n = x.size();
    auto clamp_vec = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t j = 0; j < n; j++)
          v[j] =
              v[j] < lower[j] ? lower[j] : (v[j] > upper[j] ? upper[j] : v[j]);
      }
    };
    std::vector<std::vector<scalar_t>> agents;
    if constexpr (constrained) {
      // standard DE initialization: spread the population uniformly across the
      // feasible box rather than clustering it around the starting point
      agents.assign(this->pop_size, std::vector<scalar_t>(n));
      for (size_t i = 0; i < this->pop_size; i++)
        for (size_t j = 0; j < n; j++)
          agents[i][j] = lower[j] + this->generator() * (upper[j] - lower[j]);
    } else {
      agents = init_agents(x, this->generator, this->pop_size);
    }
    std::array<size_t, 4> new_indices = {0, 0, 0, 0};
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
    std::vector<scalar_t> proposal_temp(x.size());

    std::vector<scalar_t> scores(agents.size());
    // evaluate all randomly generated agents (projected into the box first)
    for (size_t i = 0; i < agents.size(); i++) {
      clamp_vec(agents[i]);
      scores[i] = f_multiplier * this->f(agents[i]);
    }
    size_t function_calls_used = agents.size();
    size_t iter = 0;
    size_t best_id = 0, val_no_change = 0;
    scalar_t best_val = std::numeric_limits<scalar_t>::infinity();
    while (true) {
      // find the best agent (argmin)
      best_id = 0;
      for (size_t i = 1; i < scores.size(); i++) {
        if (scores[i] < scores[best_id]) best_id = i;
      }
      // stall detection on the best *value*. Previously this keyed on whether
      // the best *index* changed, so an in-place improvement (better score at
      // the same agent) looked like a stall and stopped DE prematurely.
      if (best_val - scores[best_id] > this->eps) {
        val_no_change = 0;
      } else {
        val_no_change++;
      }
      if (scores[best_id] < best_val) best_val = scores[best_id];
      // if agents have stabilized, return
      if (iter >= this->max_iter ||
          val_no_change >= this->best_value_no_change ||
          std_err(scores) < this->eps) {
        x = agents[best_id];
        // report the true objective value (scores carry f_multiplier)
        return solver_status<scalar_t>(f_multiplier * scores[best_id], iter,
                                       function_calls_used);
      }
      // main loop - this can in principle be parallelized
      for (size_t i = 0; i < agents.size(); i++) {
        // generate agent indices - either using the best or the current agent
        if constexpr (RecombinationType == random) {
          new_indices = generate_indices(i, agents.size(), this->generator);
        }
        if constexpr (RecombinationType == best) {
          new_indices =
              generate_indices(best_id, agents.size(), this->generator);
        }
        // create new mutate proposal
        propose_new_agent(new_indices, proposal_temp, agents,
                          this->crossover_prob, this->differential_weight,
                          this->generator);
        // project the trial vector into the box and evaluate
        clamp_vec(proposal_temp);
        const scalar_t score = f_multiplier * f(proposal_temp);
        function_calls_used++;
        // if score is better than previous score, update agent
        if (score < scores[i]) {
          for (size_t j = 0; j < proposal_temp.size(); j++) {
            agents[i][j] = proposal_temp[j];
          }
          scores[i] = score;
        }
      }
      // increment iteration counter
      iter++;
    }
  }
};

template <typename scalar_t, typename RNG>
static inline scalar_t rnorm(RNG &generator) {
  // this is not a particularly good generator, but it is 'good enough' for
  // our purposes.
  constexpr scalar_t pi_ = 3.141593;
  return sqrt(-2 * log(generator())) * cos(2 * pi_ * generator());
}

template <typename scalar_t, typename RNG, typename RandomIt>
static inline void rnorm(RNG &generator, RandomIt first, RandomIt last) {
  // this is not a particularly good generator, but it is 'good enough' for
  // our purposes.
  for (; first != last; ++first) {
    (*first) = rnorm<scalar_t, RNG>(generator);
  }
}

enum PSOType { Vanilla, Accelerated };

template <typename Callable, typename RNG, typename scalar_t = double,
          PSOType Type = Vanilla>
class PSO {
 private:
  // user supplied
  RNG &generator;
  Callable &f;
  scalar_t init_inertia, inertia;
  const scalar_t cognitive_coef, social_coef;
  std::vector<scalar_t> lower_, upper_;
  // static, derived from above
  size_t n_dim;
  // internally created
  std::vector<std::vector<scalar_t>> particle_positions, particle_velocities,
      particle_best_positions;
  std::vector<scalar_t> particle_best_values, swarm_best_position;
  scalar_t swarm_best_value;
  // bookkeeping
  size_t val_no_change, f_evals;
  // static limits
  const size_t n_particles, max_iter, best_val_no_change;
  const scalar_t eps;

 public:
  [[maybe_unused]] PSO<Callable, RNG, scalar_t, Type>(
      Callable &f, RNG &generator, const scalar_t inertia = 0.8,
      const scalar_t cognitive_coef = 1.8, const scalar_t social_coef = 1.8,
      const size_t n_particles = 10, const size_t max_iter = 5000,
      const size_t best_val_no_change = 50, const scalar_t eps = 1e-3)
      : generator(generator),
        f(f),
        inertia(inertia),
        cognitive_coef(cognitive_coef),
        social_coef(social_coef),
        n_dim(0),
        val_no_change(0),
        f_evals(0),
        n_particles(n_particles),
        max_iter(max_iter),
        best_val_no_change(best_val_no_change),
        eps(eps) {
    this->particle_positions =
        std::vector<std::vector<scalar_t>>(this->n_particles);
    if constexpr (Type == Vanilla) {
      this->particle_velocities =
          std::vector<std::vector<scalar_t>>(this->n_particles);
    }
    if constexpr (Type == Accelerated) {
      // keep track of original inertia
      this->init_inertia = inertia;
    }
    this->particle_best_positions =
        std::vector<std::vector<scalar_t>>(this->n_particles);
  }
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower(x.size());
    std::vector<scalar_t> upper(x.size());
    for (size_t i = 0; i < x.size(); i++) {
      scalar_t temp = std::abs(x[i]);
      lower[i] = -temp;
      upper[i] = temp;
    }
    this->init_solver_state(lower, upper);
    return this->solve<true, false>(x);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower(x.size());
    std::vector<scalar_t> upper(x.size());
    for (size_t i = 0; i < x.size(); i++) {
      scalar_t temp = std::abs(x[i]);
      lower[i] = -temp;
      upper[i] = temp;
    }
    this->init_solver_state(lower, upper);
    return this->solve<false, false>(x);
  }
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    this->init_solver_state(lower, upper);
    return this->solve<true, true>(x);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    this->init_solver_state(lower, upper);
    return this->solve<false, true>(x);
  }

 private:
  template <const bool minimize = true, const bool constrained = true>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x) {
    size_t iter = 0;
    this->update_best_positions<minimize>();
    while (true) {
      // if particles have stabilized (no improvement in objective iteration or
      // no heterogeneity of particles) or we are over the limit, return
      if (iter >= this->max_iter || val_no_change >= best_val_no_change ||
          std_err(this->particle_best_values) < this->eps) {
        x = this->swarm_best_position;
        // best scores, iteration number and function calls used total
        return solver_status<scalar_t>(this->swarm_best_value, iter,
                                       this->f_evals);
      }
      if constexpr (Type == Vanilla) {
        // Vanilla velocity update
        this->update_velocities();
      }
      if constexpr (Type == Accelerated) {
        // update inertia - we might want to create a nicer way to do this
        // updating schedule... maybe a functor for it too?
        this->inertia = pow(this->init_inertia, iter);
        // for accelerated pso update_positions also updated velocities
      }
      this->update_positions();
      if constexpr (constrained) {
        this->threshold_positions();
      }
      this->update_best_positions<minimize>();
      // increment iteration counter
      iter++;
    }
  }
  // for repeated initializations we will initialize solver with new bounds
  void init_solver_state(const std::vector<scalar_t> &lower,
                         const std::vector<scalar_t> &upper) {
    this->n_dim = lower.size();
    this->upper_ = upper;
    this->lower_ = lower;
    this->swarm_best_value = std::numeric_limits<scalar_t>::max();
    this->f_evals = 0;
    this->val_no_change = 0;
    // create particles
    for (size_t i = 0; i < this->n_particles; i++) {
      this->particle_positions[i] = std::vector<scalar_t>(this->n_dim);
      if constexpr (Type == Vanilla) {
        this->particle_velocities[i] = std::vector<scalar_t>(this->n_dim);
      }
      this->particle_best_positions[i] = std::vector<scalar_t>(this->n_dim);
    }
    for (size_t i = 0; i < n_particles; i++) {
      for (size_t j = 0; j < this->n_dim; j++) {
        // update velocities and positions
        scalar_t temp = std::abs(upper[j] - lower[j]);
        this->particle_positions[i][j] =
            lower[j] + ((upper[j] - lower[j]) * this->generator());
        if constexpr (Type == Vanilla) {
          this->particle_velocities[i][j] =
              -temp + (this->generator() * 2.0 * temp);
        }
        // update particle best positions
        this->particle_best_positions[i][j] = this->particle_positions[i][j];
      }
    }
    this->particle_best_values = std::vector<scalar_t>(
        this->n_particles, std::numeric_limits<scalar_t>::max());
    // ensure swarm_best_position is always a valid, correctly-sized vector even
    // if no particle improves on the (now +inf) initial best on the first pass
    this->swarm_best_position = this->particle_positions[0];
  }
  void update_velocities() {
    // scalar_t r_p = 0, r_g = 0;
    for (size_t i = 0; i < this->n_particles; i++) {
      for (size_t j = 0; j < this->n_dim; j++) {
        // generate random movements
        const scalar_t r_p = generator(), r_g = generator();
        // update current velocity for current particle - inertia update
        this->particle_velocities[i][j] =
            (this->inertia * this->particle_velocities[i][j]) +
            // cognitive update (moving more if further away from 'best'
            // position of particle )
            this->cognitive_coef * r_p *
                (particle_best_positions[i][j] - particle_positions[i][j]) +
            // social update (moving more if further away from 'best' position
            // of swarm)
            this->social_coef * r_g *
                (this->swarm_best_position[j] - particle_positions[i][j]);
      }
    }
  }
  void update_positions() {
    if constexpr (Type == Vanilla) {
      for (size_t i = 0; i < this->n_particles; i++) {
        for (size_t j = 0; j < this->n_dim; j++) {
          // update positions using current velocity
          this->particle_positions[i][j] += this->particle_velocities[i][j];
        }
      }
    }
    if constexpr (Type == Accelerated) {
      // accelerated PSO (Yang): x <- (1-beta) x + beta g* + alpha * span * eps.
      // beta in (0,1) is the pull toward the global best (a convex combination,
      // not the old -0.8 x + 1.8 g* which overshot); alpha is the decaying
      // inertia schedule and the random step is scaled by the domain width.
      const scalar_t beta =
          this->social_coef / (this->cognitive_coef + this->social_coef);
      for (size_t i = 0; i < this->n_particles; i++) {
        for (size_t j = 0; j < this->n_dim; j++) {
          const scalar_t span = this->upper_[j] - this->lower_[j];
          this->particle_positions[i][j] =
              (1.0 - beta) * this->particle_positions[i][j] +
              beta * this->swarm_best_position[j] +
              this->inertia * span * rnorm<scalar_t>(this->generator);
        }
      }
    }
  }
  void threshold_positions() {
    for (size_t i = 0; i < this->n_particles; i++) {
      for (size_t j = 0; j < this->n_dim; j++) {
        // threshold velocities between lower and upper
        this->particle_positions[i][j] =
            this->particle_positions[i][j] < this->lower_[j]
                ? this->lower_[j]
                : this->particle_positions[i][j];
        this->particle_positions[i][j] =
            this->particle_positions[i][j] > this->upper_[j]
                ? this->upper_[j]
                : this->particle_positions[i][j];
      }
    }
  }
  template <const bool minimize = true>
  void update_best_positions() {
    size_t best_index = 0;
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
    bool update_happened = false;
    for (size_t i = 0; i < this->n_particles; i++) {
      scalar_t temp = f_multiplier * f(particle_positions[i]);
      if (temp < this->swarm_best_value) {
        this->swarm_best_value = temp;
        // save update of swarm best position for after the loop, so we do not
        // by chance do many copies here
        best_index = i;
        update_happened = true;
      }
      if (temp < this->particle_best_values[i]) {
        this->particle_best_values[i] = temp;
      }
    }
    this->f_evals += n_particles;
    if (update_happened) {
      this->swarm_best_position = this->particle_positions[best_index];
    }
    // either increment to indicate no change in the best objective value,
    // or reset to 0
    this->val_no_change = update_happened ? 0 : (this->val_no_change + 1);
  }
};

template <typename Callable, typename RNG, typename scalar_t = double>
class SANN {
 private:
  // user supplied
  RNG &generator;
  Callable &f;
  // bookkeeping
  size_t f_evals;
  // static limits
  const size_t max_iter, temperature_iter;
  const scalar_t temperature_max;

 public:
  [[maybe_unused]] SANN<Callable, RNG, scalar_t>(
      Callable &f, RNG &generator, const size_t max_iter = 5000,
      const size_t temperature_iter = 10, const scalar_t temperature_max = 10.0)
      : generator(generator),
        f(f),
        f_evals(0),
        max_iter(max_iter),
        temperature_iter(temperature_iter),
        temperature_max(temperature_max) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0,
                       e_minus_1 = 1.7182818;
    const size_t n_dim = x.size();
    auto clamp_pt = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n_dim; i++)
          v[i] =
              v[i] < lower[i] ? lower[i] : (v[i] > upper[i] ? upper[i] : v[i]);
      }
    };
    clamp_pt(x);
    scalar_t best_val = f_multiplier * this->f(x), current_energy = best_val,
             scale = 1.0 / temperature_max;
    this->f_evals++;
    std::vector<scalar_t> p = x, ptry = x;
    size_t iter = 0;
    while (true) {
      if (iter >= this->max_iter) {
        // report the true objective value (best_val carries f_multiplier)
        return solver_status<scalar_t>(f_multiplier * best_val, iter,
                                       this->f_evals);
      }
      // temperature annealing schedule - cooling
      const scalar_t t =
          temperature_max / std::log(static_cast<scalar_t>(iter) + e_minus_1);
      for (size_t j = 0; j < this->temperature_iter; j++) {
        const scalar_t current_scale = t * scale;
        // use random normal variates - this should allow user specified values
        for (size_t i = 0; i < n_dim; i++) {
          // generate new candidate function values
          ptry[i] = p[i] + current_scale * rnorm<scalar_t>(this->generator);
        }
        clamp_pt(ptry);  // keep the candidate inside the box
        const scalar_t current_val = f_multiplier * f(ptry);
        this->f_evals++;
        // Metropolis: accept if the candidate improves on the *current* state's
        // energy, or probabilistically when uphill. (Previously this compared
        // against the best-ever value, so an accepted uphill move never moved
        // the reference energy and the chain could not ratchet up.)
        const scalar_t difference = current_val - current_energy;
        if ((difference <= 0.0) || (this->generator() < exp(-difference / t))) {
          for (size_t k = 0; k < n_dim; k++) p[k] = ptry[k];
          current_energy = current_val;
          // track the best point seen so far separately
          if (current_val < best_val) {
            for (size_t k = 0; k < n_dim; k++) x[k] = p[k];
            best_val = current_val;
          }
        }
      }
      iter++;
    }
  }
};
enum GradientStepType {
  Linesearch,
  Fixed,
  Bigstep,
  Anneal,
  PAGE
  // Momentum
};

template <const size_t level>
constexpr size_t bigstep_offset() {
  if constexpr (level == 1) return 0;
  if constexpr (level == 2) return 2;
  if constexpr (level == 3) return 5;
  if constexpr (level == 4) return 12;
  if constexpr (level == 5) return 27;
  if constexpr (level == 6) return 58;
  if constexpr (level == 7) return 121;
  return 0;
}

template <const size_t level>
constexpr size_t bigstep_len() {
  if constexpr (level == 1) return 2;
  if constexpr (level == 2) return 3;
  if constexpr (level == 3) return 7;
  if constexpr (level == 4) return 15;
  if constexpr (level == 5) return 31;
  if constexpr (level == 6) return 63;
  if constexpr (level == 7) return 127;
  return 0;
}
template <typename Callable, typename scalar_t, const size_t accuracy = 1>
struct fin_diff {
  // `accuracy` selects the central stencil: 0 is two-point (2n evaluations per
  // gradient), 1 is four-point (4n); see finite_difference_gradient
  void operator()(Callable &f, std::vector<scalar_t> &x,
                  std::vector<scalar_t> &gradient) {
    nlsolver::finite_difference::finite_difference_gradient<Callable, scalar_t,
                                                            accuracy>(f, x,
                                                                      gradient);
  }
};
// Forward differences, n + 1 evaluations per gradient, with the step
// sqrt(eps) * max(1, |x_d|) at which the O(h f'') truncation error and the
// eps / h rounding error balance near 1e-8 relative: the right tool when an
// evaluation budget is too tight for the 2n of a central stencil and cannot
// carry a polish beyond about 1e-6 anyway.
template <typename Callable, typename scalar_t>
struct fin_diff_fwd {
  void operator()(Callable &f, std::vector<scalar_t> &x,
                  std::vector<scalar_t> &gradient) {
    const scalar_t f0 = f(x);
    const scalar_t root_eps =
        std::sqrt(std::numeric_limits<scalar_t>::epsilon());
    for (size_t d = 0; d < x.size(); d++) {
      const scalar_t h = root_eps * std::max<scalar_t>(1.0, std::abs(x[d]));
      const scalar_t saved = x[d];
      x[d] = saved + h;
      gradient[d] = (f(x) - f0) / h;
      x[d] = saved;
    }
  }
};
template <typename T>
struct is_fin_diff_fwd : std::false_type {};
template <typename Callable, typename scalar_t>
struct is_fin_diff_fwd<fin_diff_fwd<Callable, scalar_t>> : std::true_type {};
// recognises any fin_diff instantiation and exposes its stencil
template <typename T>
struct is_fin_diff : std::false_type {};
template <typename Callable, typename scalar_t, const size_t accuracy>
struct is_fin_diff<fin_diff<Callable, scalar_t, accuracy>> : std::true_type {
  static constexpr size_t stencil = accuracy;
};
template <typename Callable, typename scalar_t>
struct fin_diff_h {
  void operator()(Callable &f, std::vector<scalar_t> &x,
                  std::vector<scalar_t> &hessian) {
    nlsolver::finite_difference::finite_difference_hessian<Callable, scalar_t,
                                                           1>(f, x, hessian);
  }
};
template <typename Callable, typename scalar_t,
          const GradientStepType step = GradientStepType::Fixed,
          const size_t bigstep_level = 5,
          const bool grad_norm_lipschitz_scaling = true,
          typename Grad = fin_diff<Callable, scalar_t>>
class GradientDescent {
  Callable &f;
  Grad g;
  const size_t max_iter, minibatch, minibatch_prime;
  const scalar_t grad_eps, alpha;
  nlsolver::rng::xorshift<scalar_t> generator;
  constexpr static std::array<scalar_t, 248> fixed_steps = {
      2.9, 1.5,                                 // pattern length 2 => type 1
      1.5, 4.9,  1.5,                           // type 2
      1.5, 2.2,  1.5,  12.0, 1.5,  2.2,   1.5,  // type 3
      1.4, 2.0,  1.4,  4.5,  1.4,  2.0,   1.4, 29.7, 1.4,  2.0, 1.4, 4.5,   1.4,
      2.0, 1.4,  // type 4
      1.4, 2.0,  1.4,  3.9,  1.4,  2.0,   1.4, 8.2,  1.4,  2.0, 1.4, 3.9,   1.4,
      2.0, 1.4,  72.3, 1.4,  2.0,  1.4,   3.9, 1.4,  2.0,  1.4, 8.2, 1.4,   2.0,
      1.4, 3.9,  1.4,  2.0,  1.4,  // type 5
      1.4, 2.0,  1.4,  3.9,  1.4,  2.0,   1.4, 7.2,  1.4,  2.0, 1.4, 3.9,   1.4,
      2.0, 1.4,  14.2, 1.4,  2.0,  1.4,   3.9, 1.4,  2.0,  1.4, 7.2, 1.4,   2.0,
      1.4, 3.9,  1.4,  2.0,  1.4,  164.0, 1.4, 2.0,  1.4,  3.9, 1.4, 2.0,   1.4,
      7.2, 1.4,  2.0,  1.4,  3.9,  1.4,   2.0, 1.4,  14.2, 1.4, 2.0, 1.4,   3.9,
      1.4, 2.0,  1.4,  7.2,  1.4,  2.0,   1.4, 3.9,  1.4,  2.0, 1.4,  // type 6
      1.4, 2.0,  1.4,  3.9,  1.4,  2.0,   1.4, 7.2,  1.4,  2.0, 1.4, 3.9,   1.4,
      2.0, 1.4,  12.6, 1.4,  2.0,  1.4,   3.9, 1.4,  2.0,  1.4, 7.2, 1.4,   2.0,
      1.4, 3.9,  1.4,  2.0,  1.4,  23.5,  1.4, 2.0,  1.4,  3.9, 1.4, 2.0,   1.4,
      7.2, 1.4,  2.0,  1.4,  3.9,  1.4,   2.0, 1.4,  12.6, 1.4, 2.0, 1.4,   3.9,
      1.4, 2.0,  1.4,  7.2,  1.4,  2.0,   1.4, 3.9,  1.4,  2.0, 1.4, 370.0, 1.4,
      2.0, 1.4,  3.9,  1.4,  2.0,  1.4,   7.2, 1.4,  2.0,  1.4, 3.9, 1.4,   2.0,
      1.4, 12.6, 1.4,  2.0,  1.4,  3.9,   1.4, 2.0,  1.4,  7.2, 1.4, 2.0,   1.4,
      3.9, 1.4,  2.0,  1.4,  23.5, 1.4,   2.0, 1.4,  3.9,  1.4, 2.0, 1.4,   7.5,
      1.4, 2.0,  1.4,  3.9,  1.4,  2.0,   1.4, 12.6, 1.4,  2.0, 1.4, 3.9,   1.4,
      2.0, 1.4,  7.2,  1.4,  2.0,  1.4,   3.9, 1.4,  2.0,  1.4  // type 7
  };
  std::vector<scalar_t> search_direction, linesearch_temp, gradient_temp;

 public:
  explicit GradientDescent<Callable, scalar_t, step, bigstep_level,
                           grad_norm_lipschitz_scaling, Grad>(
      Callable &f, const scalar_t alpha = 1, const size_t max_iter = 500,
      const scalar_t grad_eps = 1e-12, const size_t minibatch_b = 128,
      const size_t minibatch_b_prime = 11,
      Grad g = fin_diff<Callable, scalar_t>())
      : f(f),
        g(g),
        max_iter(max_iter),
        minibatch(minibatch_b),
        minibatch_prime(minibatch_b_prime),
        grad_eps(grad_eps),
        alpha(alpha),
        generator(nlsolver::rng::xorshift<scalar_t>()) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces (projected gradient); bounds are (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    const size_t n_dim = x.size();
    int i_dim = static_cast<int>(n_dim);
    auto project = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n_dim; i++)
          v[i] =
              v[i] < lower[i] ? lower[i] : (v[i] > upper[i] ? upper[i] : v[i]);
      }
    };
    project(x);
    std::vector<scalar_t> gradient = std::vector<scalar_t>(n_dim, 0.0),
                          prev_gradient = std::vector<scalar_t>(n_dim, 0.0);
    if constexpr (step == GradientStepType::Linesearch) {
      // we need additional temporaries for linesearch
      this->search_direction = std::vector<scalar_t>(n_dim, 0.0);
      this->linesearch_temp = std::vector<scalar_t>(n_dim, 0.0);
      this->gradient_temp = std::vector<scalar_t>(n_dim, 0.0);
    }
    scalar_t alpha_ = this->alpha;
    size_t iter = 0, function_calls_used = 0, grad_evals_used = 0;
    constexpr scalar_t f_multiplier = minimize ? -1.0 : 1.0;
    scalar_t max_grad_norm = 0;
    // only necessary and interesting for PAGE (float division - this was
    // size_t/size_t and always truncated to 0, so the resample branch never
    // fired)
    const scalar_t p = static_cast<scalar_t>(minibatch) /
                       static_cast<scalar_t>(minibatch_prime + minibatch);
    const scalar_t ratio = static_cast<scalar_t>(minibatch) /
                           static_cast<scalar_t>(minibatch_prime);
    // construct lambda that takes f and enables function evaluation counting
    auto f_lam = [&](decltype(x) &coef) {
      function_calls_used++;
      return this->f(coef);
    };
    auto g_lam = [&](decltype(x) &coef, decltype(gradient) &grad) {
      grad_evals_used++;
      // simple optimization - finite difference gradient is actually stateless
      // so in that case this wrapper is entirely valid
      if constexpr (std::is_same<Grad, fin_diff<Callable, scalar_t>>::value) {
        fin_diff<decltype(f_lam), scalar_t>()(f_lam, coef, grad);
        return;
      }
      /* otherwise we cannot keep track of function evaluations that way
       * and our users will have to implement their own gradient function (that
       * might be smarter, actually) - we can however implement our own counter
       * for gradient evaluations
       */
      this->g(this->f, coef, grad);
      return;
    };
    // compute gradient
    g_lam(x, gradient);
    while (true) {
      const scalar_t grad_norm = math::norm(gradient.data(), i_dim);
      max_grad_norm = std::max(max_grad_norm, grad_norm);
      if (iter >= this->max_iter || grad_norm < grad_eps ||
          std::isinf(grad_norm)) {
        // evaluate at current parameters
        scalar_t current_val = f_lam(x);
        return solver_status<scalar_t>(current_val, iter, function_calls_used,
                                       grad_evals_used);
      }
      if constexpr (step == GradientStepType::Linesearch) {
        nlsolver::math::a_mult_scalar_to_b(gradient.data(), f_multiplier,
                                           this->search_direction.data(),
                                           i_dim);
        for (size_t i = 0; i < n_dim; i++) {
          // this->search_direction[i] = f_multiplier * gradient[i];
          this->gradient_temp[i] = gradient[i];
        }
        alpha_ = nlsolver::linesearch::more_thuente_search(
            f_lam, x, this->gradient_temp, this->search_direction,
            this->linesearch_temp, this->alpha, g_lam);
      }
      // fixed step, optionally normalized by the gradient magnitude
      if constexpr (step == GradientStepType::Fixed) {
        alpha_ = this->alpha;
        if constexpr (grad_norm_lipschitz_scaling) alpha_ /= max_grad_norm;
      }
      if constexpr (step == GradientStepType::Anneal) {
        // cooling schedule; scaled by the gradient magnitude so the raw step is
        // not an unscaled ~alpha (which diverged/oscillated on most problems)
        alpha_ = this->alpha / (1.0 + (static_cast<scalar_t>(iter) / max_iter));
        if constexpr (grad_norm_lipschitz_scaling) alpha_ /= max_grad_norm;
      }
      if constexpr (step == GradientStepType::PAGE) {
        // Lipschitz-scaled step keeps the variance-reduced gradient bounded
        // (without this PAGE diverged to ~1e+300)
        alpha_ = this->alpha;
        if constexpr (grad_norm_lipschitz_scaling) alpha_ /= max_grad_norm;
      }
      if constexpr (step == GradientStepType::Bigstep) {
        constexpr size_t offset = bigstep_offset<bigstep_level>();
        constexpr size_t step_len = bigstep_len<bigstep_level>();
        const size_t current_step = offset + iter % step_len;
        if constexpr (bigstep_level == 0) {
          alpha_ = ((current_step == 0) * (fixed_steps[current_step] - alpha)) +
                   ((current_step != 0) * fixed_steps[current_step]);
        }
        if constexpr (bigstep_level != 0) {
          alpha_ = fixed_steps[current_step];
        }
        if constexpr (grad_norm_lipschitz_scaling) {
          alpha_ /= max_grad_norm;
        }
      }
      // update parameters: x += (f_multiplier * alpha_) * gradient. Use a local
      // signed step so alpha_ is never sign-flipped/persisted across iterations
      // (which previously alternated the step direction for Fixed/PAGE).
      const scalar_t signed_step = f_multiplier * alpha_;
      nlsolver::math::a_mult_scalar_add_b(gradient.data(), signed_step,
                                          x.data(), i_dim);
      project(x);  // projected gradient: keep x inside the box
      if constexpr (step == GradientStepType::PAGE) {
        for (size_t i = 0; i < n_dim; i++) prev_gradient[i] = gradient[i];
      }
      // compute gradient
      g_lam(x, gradient);
      if constexpr (step == GradientStepType::PAGE) {
        if (generator() > p) {
          // only do a small update where new gradient is old gradient
          // + difference between gradients
          nlsolver::math::a_minus_b_mult_scalar_add_c(
              gradient.data(), prev_gradient.data(), ratio, gradient.data(),
              i_dim);
        }
      }
      iter++;
    }
  }
};

template <typename Callable, typename scalar_t,
          typename Grad = fin_diff<Callable, scalar_t>>
class ConjugatedGradientDescent {
  Callable &f;
  Grad g;
  const size_t max_iter;
  const scalar_t grad_eps, alpha;

 public:
  explicit ConjugatedGradientDescent<Callable, scalar_t, Grad>(
      Callable &f, Grad g = fin_diff<Callable, scalar_t>(),
      const size_t max_iter = 500, const scalar_t grad_eps = 5e-3,
      const scalar_t alpha =
          1.0)  // line-search initial step (armijo shrinks it)
      : f(f), g(g), max_iter(max_iter), grad_eps(grad_eps), alpha(alpha) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces (projected gradient); bounds are (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    const size_t n_dim = x.size();
    int i_dim = static_cast<int>(n_dim);
    auto project = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n_dim; i++)
          v[i] =
              v[i] < lower[i] ? lower[i] : (v[i] > upper[i] ? upper[i] : v[i]);
      }
    };
    project(x);
    std::vector<scalar_t> gradient = std::vector<scalar_t>(n_dim, 0.0);
    // we need additional temporaries for linesearch
    std::vector<scalar_t> search_direction = std::vector<scalar_t>(n_dim, 0.0);
    std::vector<scalar_t> linesearch_temp = std::vector<scalar_t>(n_dim, 0.0);
    scalar_t alpha_ = this->alpha;
    size_t iter = 0, function_calls_used = 0, grad_evals_used = 0;
    constexpr scalar_t f_multiplier = minimize ? -1.0 : 1.0;
    // construct lambda that takes f and enables function evaluation counting
    auto f_lam = [&](decltype(x) &coef) {
      function_calls_used++;
      return this->f(coef);
    };
    auto g_lam = [&](decltype(x) &coef, decltype(gradient) &grad) {
      grad_evals_used++;
      // simple optimization - finite difference gradient is actually stateless
      // so in that case this wrapper is entirely valid
      if constexpr (std::is_same<Grad, fin_diff<Callable, scalar_t>>::value) {
        fin_diff<decltype(f_lam), scalar_t>()(f_lam, coef, grad);
        return;
      }
      /* otherwise we cannot keep track of function evaluations that way
       * and our users will have to implement their own gradient function (that
       * might be smarter, actually) - we can however implement our own counter
       * for gradient evaluations
       */
      this->g(this->f, coef, grad);
      return;
    };
    g_lam(x, gradient);
    // set search direction for linesearch
    nlsolver::math::a_mult_scalar_to_b(gradient.data(), f_multiplier,
                                       search_direction.data(), i_dim);
    size_t since_restart = 0;
    while (true) {
      // compute gradient
      const scalar_t grad_norm = math::norm(gradient.data(), i_dim);
      if (iter >= this->max_iter || grad_norm < grad_eps ||
          std::isinf(grad_norm)) {
        // evaluate at current parameters
        scalar_t current_val = f_lam(x);
        return solver_status<scalar_t>(current_val, iter, function_calls_used,
                                       grad_evals_used);
      }
      alpha_ = nlsolver::linesearch::armijo_search(
          f_lam, x, gradient, search_direction, linesearch_temp, this->alpha);
      // update parameters; x[i] += search_direction[i] * alpha_;
      nlsolver::math::a_mult_scalar_add_b(search_direction.data(), alpha_,
                                          x.data(), i_dim);
      project(x);  // projected gradient: keep x inside the box
      // recompute gradient, compute new search direction using conjugation
      // first, compute gradient.dot(gradient) with existing gradient,
      // then compute new gradient and compute the same, then compute their
      // ratio
      scalar_t denominator = math::dot(gradient.data(), gradient.data(), i_dim);
      g_lam(x, gradient);
      // figure out the numerator from new gradient
      scalar_t numerator = math::dot(gradient.data(), gradient.data(), i_dim);
      // Fletcher-Reeves beta, guarded against a vanishing denominator
      scalar_t search_update =
          denominator > 1e-30 ? numerator / denominator : 0.0;
      // periodic restart to steepest descent - Fletcher-Reeves CG stalls
      // (loses conjugacy) without it
      if (++since_restart >= n_dim) {
        search_update = 0.0;
        since_restart = 0;
      }
      // update search direction
      nlsolver::math::a_mul_scalar(search_direction.data(), search_update,
                                   i_dim);
      nlsolver::math::a_mult_scalar_add_b(gradient.data(), f_multiplier,
                                          search_direction.data(), i_dim);
      iter++;
    }
  }
};
template <typename scalar_t>
void update_inverse_hessian(std::vector<scalar_t> &inv_hessian,
                            const std::vector<scalar_t> &step,
                            const std::vector<scalar_t> &grad_diff,
                            std::vector<scalar_t> &grad_diff_inv_hess,
                            const scalar_t rho) {
  // precompute temporaries needed in the hessian update
  const size_t n_dim = grad_diff.size();
  int i_dim = static_cast<int>(n_dim);
  for (size_t i = 0; i < n_dim; i++) {
    grad_diff_inv_hess[i] =
        math::dot(grad_diff.data(), inv_hessian.data() + (i * n_dim), i_dim);
  }
  // BFGS inverse update, rho = 1 / (y's):
  //   H+ = (I - rho s y') H (I - rho y s') + rho s s'
  //      = H - rho (s (Hy)' + (Hy) s') + rho (rho y'Hy + 1) s s'
  // The rank-one s s' term is ADDED; subtracting it (as this once did) breaks
  // positive definiteness even when y's > 0, so the quasi-Newton direction
  // stops being a descent direction every few iterations and the method
  // collapses to steepest descent.
  const scalar_t y_h_y =
      math::dot(grad_diff.data(), grad_diff_inv_hess.data(), i_dim);
  const scalar_t ss_coef = rho * (rho * y_h_y + 1.0);
  for (size_t j = 0; j < n_dim; j++) {
    for (size_t i = 0; i < n_dim; i++) {
      inv_hessian[j * n_dim + i] +=
          ss_coef * step[i] * step[j] - rho * (step[i] * grad_diff_inv_hess[j] +
                                               grad_diff_inv_hess[i] * step[j]);
    }
  }
}
// Line search of BFGS. More-Thuente enforces the strong Wolfe conditions and
// evaluates a gradient at every trial step, which with finite differences
// makes each trial cost 2n + 1 evaluations; Armijo backtracking needs only
// function values per trial and relies on the curvature-pair test (y's > 0)
// to keep the inverse Hessian positive definite. Armijo is the default: on
// the smooth suite problems it reaches 1e-10 in the same evaluations as
// L-BFGS-B, two to three times fewer than More-Thuente.
enum BFGSLineSearch { MoreThuente, Armijo };

template <typename Callable, typename scalar_t = double,
          typename Grad = fin_diff<Callable, scalar_t, 0>>
class BFGS {
 private:
  Callable &f;
  Grad g;
  // stopping
  const size_t max_iter;
  const scalar_t grad_eps, alpha;
  const BFGSLineSearch line_search;

 public:
  // constructor
  explicit BFGS<Callable, scalar_t, Grad>(
      Callable &f, Grad g = Grad(), const size_t max_iter = 200,
      const scalar_t grad_eps = 1e-7, const scalar_t alpha = 1,
      const BFGSLineSearch line_search = Armijo)
      : f(f),
        g(g),
        max_iter(max_iter),
        grad_eps(grad_eps),
        alpha(alpha),
        line_search(line_search) {}

  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces (projected); bounds are (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    // maximize f by minimizing -f (the library-wide convention); every
    // function/gradient value flows through f_lam/g_lam below, so the whole
    // method - line search, Hessian update - operates consistently on sign*f.
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const size_t n_dim = x.size();
    auto project = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n_dim; i++)
          v[i] =
              v[i] < lower[i] ? lower[i] : (v[i] > upper[i] ? upper[i] : v[i]);
      }
    };
    project(x);
    int i_dim = static_cast<int>(n_dim);
    std::vector<scalar_t> inverse_hessian =
        std::vector<scalar_t>(n_dim * n_dim);
    std::vector<scalar_t> search_direction = std::vector<scalar_t>(n_dim, 0.0),
                          gradient = std::vector<scalar_t>(n_dim, 0.0),
                          prev_gradient = std::vector<scalar_t>(n_dim, 0.0),
                          grad_update = std::vector<scalar_t>(n_dim, 0.0),
                          s = std::vector<scalar_t>(n_dim, 0.0),
                          linesearch_temp = std::vector<scalar_t>(n_dim, 0.0),
                          grad_diff_inv_hess(n_dim);
    // initialize to identity matrix
    for (size_t i = 0; i < n_dim; i++) inverse_hessian[i + (i * n_dim)] = 1.0;
    // an identity H carries no scale; the first curvature pair after a (re)set
    // rescales it to (s'y / y'y) I before the update (Nocedal & Wright, eq.
    // 6.20), otherwise the unexplored subspace keeps unit curvature and the
    // line search must shrink every step on a stiff problem
    bool fresh_hessian = true;
    size_t iter = 0, function_calls_used = 0, grad_evals_used = 0;
    auto f_lam = [&](decltype(x) &coef) {
      function_calls_used++;
      return sign * this->f(coef);
    };
    auto g_lam = [&](decltype(x) &coef, decltype(gradient) &grad) {
      grad_evals_used++;
      // simple optimization - finite difference gradient is actually stateless
      // so in that case this wrapper is entirely valid
      if constexpr (is_fin_diff<Grad>::value) {
        // f_lam already carries the sign, so the finite-difference gradient is
        // the gradient of sign*f, with the stencil the user selected
        fin_diff<decltype(f_lam), scalar_t, is_fin_diff<Grad>::stencil>()(
            f_lam, coef, grad);
        return;
      }
      /* otherwise we cannot keep track of function evaluations that way
       * and our users will have to implement their own gradient function (that
       * might be smarter, actually) - we can however implement our own counter
       * for gradient evaluations
       */
      this->g(this->f, coef, grad);
      if constexpr (!minimize) {
        for (auto &gi : grad) gi = -gi;  // gradient of -f
      }
      return;
    };
    g_lam(x, gradient);
    scalar_t current_grad_norm = nlsolver::math::norm(gradient.data(), i_dim);
    // objective at the iterate, carried along by the Armijo line search
    scalar_t fx = f_lam(x);
    std::vector<scalar_t> x_new(n_dim, 0.0);
    auto may_evaluate = []() { return true; };
    // Stop on a small gradient, or once a step moves x by less than
    // sqrt(eps) relative (a failed line search returns a zero step), which is
    // the finite-difference resolution below which no further progress can be
    // verified. A test on the change of the gradient norm between iterations
    // is not a convergence test: in a curved valley the norm changes slowly
    // while the point is still far from the optimum.
    const scalar_t min_rel_step =
        std::sqrt(std::numeric_limits<scalar_t>::epsilon());
    scalar_t step_inf = std::numeric_limits<scalar_t>::infinity();
    while (true) {
      scalar_t x_inf = 0;
      for (size_t j = 0; j < n_dim; j++)
        x_inf = std::max(x_inf, std::abs(x[j]));
      if (iter >= this->max_iter || current_grad_norm < grad_eps ||
          step_inf <= min_rel_step * (1.0 + x_inf) ||
          !std::isfinite(current_grad_norm)) {
        // report the true (un-negated) value; the Armijo path carries f(x),
        // the More-Thuente path evaluates it
        const scalar_t current_val =
            this->line_search == Armijo ? sign * fx : sign * f_lam(x);
        return solver_status<scalar_t>(current_val, iter, function_calls_used,
                                       grad_evals_used);
      }
      // update search direction vector using -inverse_hessian * gradient
      for (size_t j = 0; j < n_dim; j++) {
        search_direction[j] = -math::dot(inverse_hessian.data() + (j * n_dim),
                                         gradient.data(), i_dim);
      }
      // Reset the inverse Hessian only when the model direction is not a
      // descent direction. Resetting on an increase of the gradient norm
      // discards the curvature model on every step along a curved valley,
      // where the norm is not monotone, and degrades the method to steepest
      // descent.
      scalar_t phi = math::dot(gradient.data(), search_direction.data(), i_dim);
      if ((phi >= 0) || std::isnan(phi)) {
        fresh_hessian = true;
        std::fill(inverse_hessian.begin(), inverse_hessian.end(), 0.0);
        // reset hessian approximation and search_direction
        for (size_t i = 0; i < n_dim; i++) {
          inverse_hessian[i + (i * n_dim)] = 1.0;
          search_direction[i] = -gradient[i];
        }
      }
      prev_gradient = gradient;
      if (this->line_search == Armijo) {
        // unit step once curvature is known; the first step after a (re)set is
        // scaled by the gradient magnitude, as in L-BFGS-B
        scalar_t g_inf = 0;
        for (size_t j = 0; j < n_dim; j++) {
          g_inf = std::max(g_inf, std::abs(gradient[j]));
        }
        const scalar_t step =
            fresh_hessian ? 1.0 / std::max(static_cast<scalar_t>(1), g_inf)
                          : 1.0;
        scalar_t f_new = fx;
        const bool accepted =
            armijo_backtrack(f_lam, x, fx, gradient, search_direction, step,
                             project, may_evaluate, x_new, f_new);
        if (!accepted) {
          return solver_status<scalar_t>(sign * fx, iter, function_calls_used,
                                         grad_evals_used);
        }
        fx = f_new;
        // the projected displacement is the curvature-pair step
        for (size_t j = 0; j < n_dim; j++) s[j] = x_new[j] - x[j];
        x = x_new;
      } else {
        const scalar_t rate = nlsolver::linesearch::more_thuente_search(
            f_lam, x, gradient, search_direction, linesearch_temp, this->alpha,
            g_lam);
        // update parameters
        nlsolver::math::a_mult_scalar_to_b(search_direction.data(), rate,
                                           s.data(), i_dim);
        if constexpr (constrained) {
          // projected step: clamp x into the box, then make s the *actual*
          // displacement so the curvature pair (s, y) stays consistent
          for (size_t j = 0; j < n_dim; j++) {
            const scalar_t x_old = x[j];
            x[j] = x_old + s[j];
            x[j] = x[j] < lower[j] ? lower[j]
                                   : (x[j] > upper[j] ? upper[j] : x[j]);
            s[j] = x[j] - x_old;
          }
        } else {
          nlsolver::math::a_plus_b(x.data(), s.data(), i_dim);
        }
      }
      step_inf = 0;
      for (size_t j = 0; j < n_dim; j++) {
        step_inf = std::max(step_inf, std::abs(s[j]));
      }
      // we also need to compute the gradient at this new point
      // update it by reference
      g_lam(x, gradient);
      current_grad_norm = nlsolver::math::norm(gradient.data(), i_dim);
      // Update grad difference, rho and inverse hessian
      nlsolver::math::a_minus_b_to_c(gradient.data(), prev_gradient.data(),
                                     grad_update.data(), i_dim);
      const scalar_t sy =
          nlsolver::math::dot(grad_update.data(), s.data(), i_dim);
      // BFGS curvature condition: only apply the (Sherman-Morrison) inverse
      // Hessian update when y.s is sufficiently positive. Otherwise rho =
      // 1/(y.s) is +-inf or negative and corrupts the approximation, so we
      // skip the update and keep the previous inverse Hessian.
      const scalar_t y_norm = nlsolver::math::norm(grad_update.data(), i_dim);
      const scalar_t s_norm = nlsolver::math::norm(s.data(), i_dim);
      if (sy > std::numeric_limits<scalar_t>::epsilon() * y_norm * s_norm) {
        if (fresh_hessian) {
          const scalar_t scale = sy / (y_norm * y_norm);
          for (size_t i = 0; i < n_dim; i++) {
            inverse_hessian[i + (i * n_dim)] = scale;
          }
          fresh_hessian = false;
        }
        const scalar_t rho = 1 / sy;
        update_inverse_hessian(inverse_hessian, s, grad_update,
                               grad_diff_inv_hess, rho);
      }
      iter++;
    }
  }
};
template <typename Callable, typename scalar_t>
class [[maybe_unused]] Brent {
  Callable &f;
  const scalar_t tol, eps;
  const size_t max_iter;

 public:
  explicit Brent(Callable &f, const scalar_t tol = 1e-12,
                 const scalar_t eps = 1e-12, const size_t max_iter = 200)
      : f(f), tol(tol), eps(eps), max_iter(max_iter) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(scalar_t &x,
                                                    const scalar_t lower = -5,
                                                    const scalar_t upper = 5) {
    return this->solve<true>(x, upper, lower);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(scalar_t &x,
                                                    const scalar_t lower = -5,
                                                    const scalar_t upper = 5) {
    return this->solve<false>(x, upper, lower);
  }

 private:
  template <const bool minimize = true>
  solver_status<scalar_t> solve(scalar_t &x_, const scalar_t upper,
                                const scalar_t lower) {
    // lightly adapted from R's C level Brent_fmin
    // values
    size_t f_evals_used = 0;
    constexpr scalar_t f_mult = minimize ? 1.0 : -1.0;
    // functor
    auto f_lam = [&](auto x) {
      const auto val = f_mult * this->f(x);
      f_evals_used++;
      return val;
    };
    /*  c is the squared inverse of the golden ratio */
    const scalar_t c = (3. - std::sqrt(5.)) * .5;
    /* Local variables */
    scalar_t a, b, d, e, p, q, r, u, v, w, x, t2, fu, fv, fw, fx, xm, tol1,
        tol3;
    /*  eps is approximately the square root of the relative machine precision.
     */
    tol1 = eps + 1.; /* the smallest 1.000... > 1 */
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    d = 0.; /* -Wall */
    e = 0.;
    fx = f_lam(x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    /*  main loop starts here ----------------------------------- */
    size_t iter = 0;
    for (; iter < max_iter; ++iter) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
      /* check stopping criterion */
      if (std::fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (std::fabs(e) > tol1) { /* fit parabola */
        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2.;
        if (q > 0.) {
          p = -p;
        } else {
          q = -q;
        }
        r = e;
        e = d;
      }
      if (std::fabs(p) >= std::fabs(q * .5 * r) || p <= q * (a - x) ||
          p >= q * (b - x)) { /* a golden-section step */
        if (x < xm) {
          e = b - x;
        } else {
          e = a - x;
        }
        d = c * e;
      } else { /* a parabolic-interpolation step */
        d = p / q;
        u = x + d;
        /* f must not be evaluated too close to ax or bx */
        if (u - a < t2 || b - u < t2) {
          d = tol1;
          if (x >= xm) {
            d = -d;
          }
        }
      }
      /* f must not be evaluated too close to x */
      if (std::fabs(d) >= tol1) {
        u = x + d;
      } else if (d > 0.) {
        u = x + tol1;
      } else {
        u = x - tol1;
      }
      fu = f_lam(u);
      /*  update  a, b, v, w, and x */
      if (fu <= fx) {
        if (u < x)
          b = x;
        else
          a = x;
        v = w;
        w = x;
        x = u;
        fv = fw;
        fw = fx;
        fx = fu;
      } else {
        if (u < x) {
          a = u;
        } else {
          b = u;
        }
        if (fu <= fw || w == x) {
          v = w;
          fv = fw;
          w = u;
          fw = fu;
        } else if (fu <= fv || v == x || v == w) {
          v = u;
          fv = fu;
        }
      }
    }
    x_ = x;
    return solver_status<scalar_t>(fx, iter, f_evals_used);
  }
};
template <typename Callable, typename scalar_t,
          typename Grad = fin_diff<Callable, scalar_t>,
          typename Hess = fin_diff_h<Callable, scalar_t>>
class [[maybe_unused]] LevenbergMarquardt {
 private:
  Callable &f;
  Grad g;
  Hess h;
  scalar_t lambda;
  const scalar_t upward_mult, downward_mult;
  const size_t max_iter;
  const scalar_t f_delta;

 public:
  // constructor
  explicit LevenbergMarquardt<Callable, scalar_t, Grad, Hess>(
      Callable &f, const scalar_t lambda = 10, const scalar_t upward_mult = 10,
      const scalar_t downward_mult = 10, const size_t max_iter = 100,
      const scalar_t f_delta = 1e-12, Grad g = fin_diff<Callable, scalar_t>(),
      Hess h = fin_diff_h<Callable, scalar_t>())
      : f(f),
        g(g),
        h(h),
        lambda(lambda),
        upward_mult(upward_mult),
        downward_mult(downward_mult),
        max_iter(max_iter),
        f_delta(f_delta) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    return this->solve<true>(x);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    return this->solve<false>(x);
  }

 private:
  template <const bool minimize = true>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x) {
    // maximize f by minimizing -f: near a maximum of f the Hessian of -f is
    // positive-definite, so the LM normal equations stay well-posed.
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const size_t n = x.size();
    size_t iter = 0;
    // scalar_t current_f_value;
    size_t function_calls_used = 0, grad_evals_used = 0, hess_evals_used = 0;
    std::vector<scalar_t> gradient(n, 0.0);
    std::vector<scalar_t> hessian(n * n, 0.0);
    // construct lambda that takes f and enables function evaluation counting
    auto f_lam = [&](decltype(x) &coef) {
      function_calls_used++;
      return sign * this->f(coef);
    };
    auto g_lam = [&](decltype(x) &coef, decltype(gradient) &grad) {
      grad_evals_used++;
      // simple optimization - finite difference gradient is actually stateless
      // so in that case this wrapper is entirely valid
      if constexpr (std::is_same<Grad, fin_diff<Callable, scalar_t>>::value) {
        // f_lam already carries the sign
        fin_diff<decltype(f_lam), scalar_t>()(f_lam, coef, grad);
        return;
      }
      this->g(this->f, coef, grad);
      if constexpr (!minimize) {
        for (auto &gi : grad) gi = -gi;
      }
      return;
    };
    auto h_lam = [&](decltype(x) &coef, decltype(hessian) &hess_) {
      hess_evals_used++;
      // simple optimization - finite difference gradient is actually stateless
      // so in that case this wrapper is entirely valid
      if constexpr (std::is_same<Hess, fin_diff_h<Callable, scalar_t>>::value) {
        // f_lam already carries the sign
        fin_diff_h<decltype(f_lam), scalar_t>()(f_lam, coef, hess_);
        return;
      }
      this->h(this->f, coef, hess_);
      if constexpr (!minimize) {
        for (auto &hi : hess_) hi = -hi;
      }
      return;
    };
    g_lam(x, gradient);
    h_lam(x, hessian);
    scalar_t current_f_value = f_lam(x);
    std::vector<scalar_t> update(n), x_trial(n), damped_hessian(n * n);
    while (true) {
      if (iter >= max_iter || !std::isfinite(lambda) || lambda > 1e12) {
        return solver_status<scalar_t>(sign * current_f_value, iter,
                                       function_calls_used, grad_evals_used,
                                       hess_evals_used);
      }
      // damp a *copy* of the Hessian: H + lambda*I (never mutate the Hessian
      // itself, so a rejected step does not accumulate damping)
      damped_hessian = hessian;
      for (size_t i = 0; i < n; i++) damped_hessian[i * n + i] += lambda;
      // trial step: x_trial = x - (H + lambda I)^-1 g
      nlsolver::math::get_update_with_hessian(update, damped_hessian, gradient);
      for (size_t i = 0; i < n; i++) x_trial[i] = x[i] - update[i];
      const scalar_t trial_f = f_lam(x_trial);
      iter++;
      if (std::isfinite(trial_f) && trial_f < current_f_value) {
        // accept the step, decrease damping, refresh gradient/Hessian
        const scalar_t f_delta = std::abs(current_f_value - trial_f);
        x = x_trial;
        current_f_value = trial_f;
        lambda /= downward_mult;
        g_lam(x, gradient);
        h_lam(x, hessian);
        if (f_delta < this->f_delta) {
          return solver_status<scalar_t>(sign * current_f_value, iter,
                                         function_calls_used, grad_evals_used,
                                         hess_evals_used);
        }
      } else {
        // reject: keep x, increase damping (toward gradient descent) and retry
        lambda *= upward_mult;
      }
    }
  }
};
template <typename Callable, typename RNG, typename scalar_t = double>
class NelderMeadPSO {
 private:
  // user supplied
  RNG &generator;
  Callable &f;
  // initialized once
  const scalar_t alpha, gamma, rho, sigma, inertia, cognitive_coef, social_coef;
  scalar_t eps;
  // used during optimization
  std::vector<std::vector<scalar_t>> particle_positions, particle_velocities;
  std::vector<scalar_t> particle_current_values;
  const size_t max_iter, no_change_best_iter;
  size_t function_calls_used = 0;

 public:
  // constructor
  [[maybe_unused]] NelderMeadPSO(
      Callable &f, RNG &generator, const scalar_t alpha = 1,
      const scalar_t gamma = 2, const scalar_t rho = 0.5,
      const scalar_t sigma = 0.5, const scalar_t inertia = 0.8,
      const scalar_t cognitive_coef = 1.8, const scalar_t social_coef = 1.8,
      const scalar_t eps = 1e-6, const size_t max_iter = 1000,
      const size_t no_change_best_iter = 20)
      : generator(generator),
        f(f),
        alpha(alpha),
        gamma(gamma),
        rho(rho),
        sigma(sigma),
        inertia(inertia),
        cognitive_coef(cognitive_coef),
        social_coef(social_coef),
        eps(eps),
        max_iter(max_iter),
        no_change_best_iter(no_change_best_iter) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    const size_t n_dim = x.size();
    // compute implied upper and lower bounds
    std::vector<scalar_t> lower(n_dim);
    std::vector<scalar_t> upper(n_dim);
    for (size_t i = 0; i < n_dim; i++) {
      scalar_t temp = std::abs(2.5 * x[i]);
      lower[i] = -temp;
      upper[i] = temp;
    }
    return this->solve<true, false>(x, upper, lower);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    const size_t n_dim = x.size();
    // compute implied upper and lower bounds
    std::vector<scalar_t> lower(n_dim);
    std::vector<scalar_t> upper(n_dim);
    for (size_t i = 0; i < n_dim; i++) {
      scalar_t temp = std::abs(2.5 * x[i]);
      lower[i] = -temp;
      upper[i] = temp;
    }
    return this->solve<false, false>(x, upper, lower);
  }
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, upper, lower);
  }
  // maximize helper
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, upper, lower);
  }

 private:
  template <const bool minimize = true, const bool bound = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &upper,
                                const std::vector<scalar_t> &lower) {
    const size_t n_dim = x.size();
    if (n_dim < 2) {
      std::cout
          << "You are trying to optimize a one dimensional function "
          << "you should probably be using vanilla NelderMead (or vanilla PSO)"
          << " - our implementation does not support this in the "
             "NelderMead-PSO hybrid."
          << std::endl;
      // return some invalid solver state
      return solver_status<scalar_t>(999999, 0, 0);
    }
    // initialize solver state
    const size_t n_simplex_particles = n_dim + 1, n_pso_particles = 2 * n_dim,
                 n_particles = n_pso_particles + n_simplex_particles;
    init_solver_state<minimize>(x, upper, lower, n_simplex_particles,
                                n_pso_particles, n_dim);
    // temporaries for the simplex
    std::vector<scalar_t> centroid(n_dim), temp_reflect(n_dim),
        temp_expand(n_dim), temp_contract(n_dim);
    // create an index for the current particle order (best to worst)
    std::vector<size_t> current_order(n_particles);
    std::iota(current_order.begin(), current_order.end(), 0);
    size_t iter = 0;
    scalar_t best_val = this->particle_current_values[0];
    size_t best_val_no_change = 0;
    while (true) {
      // sort particles from best to worst
      std::sort(current_order.begin(), current_order.end(),
                [&](size_t left_id, size_t right_id) -> bool {
                  return this->particle_current_values[left_id] <
                         this->particle_current_values[right_id];
                });
      // stall detection on the best value (best_val must actually be updated -
      // previously it was fixed at the initial value, so this never fired)
      const scalar_t cur_best = this->particle_current_values[current_order[0]];
      if (best_val - cur_best > this->eps) {
        best_val_no_change = 0;  // meaningful improvement
      } else {
        best_val_no_change++;
      }
      if (cur_best < best_val) best_val = cur_best;
      // stopping criteria
      if (iter >= this->max_iter ||
          best_val_no_change >= this->no_change_best_iter ||
          // this should be applied only over the simplex particles
          this->simplex_std_err(current_order, n_simplex_particles) <
              this->eps) {
        x = this->particle_positions[current_order[0]];
        // best scores, iteration number and function calls used total
        return solver_status<scalar_t>(
            this->particle_current_values[current_order[0]], iter,
            this->function_calls_used);
      }
      // use top N+1 particles to form simplex and apply the simplex update
      apply_simplex<minimize, bound>(centroid, temp_reflect, temp_expand,
                                     temp_contract, current_order,
                                     n_simplex_particles, n_dim, upper, lower);
      // update the rest of the particles - i.e. the non-simplex ones
      // using the regular PSO update
      apply_pso<minimize, bound>(current_order, n_simplex_particles,
                                 n_particles, n_dim, upper, lower);
      iter++;
    }
  }
  template <const bool minimize = true>
  void init_solver_state(const std::vector<scalar_t> &x,
                         const std::vector<scalar_t> &upper,
                         const std::vector<scalar_t> &lower,
                         const size_t nm_particles, const size_t pso_particles,
                         const size_t n_dim) {
    const size_t n_particles = nm_particles + pso_particles;
    this->particle_positions = std::vector<std::vector<scalar_t>>(n_particles);
    this->particle_velocities = std::vector<std::vector<scalar_t>>(n_particles);
    this->particle_current_values = std::vector<scalar_t>(n_particles);
    this->function_calls_used = 0;
    size_t i = 0;
    for (; i < n_particles; i++) {
      this->particle_positions[i] = std::vector<scalar_t>(n_dim);
      this->particle_velocities[i] = std::vector<scalar_t>(n_dim, 0.0);
    }
    // create particles - first x.size() + 1 particles should be initialized
    // as in NM; this follows Gao and Han, see:
    // 'Proper initialization is crucial for the Nelder–Mead simplex search.'
    // (2019), Wessing, S.  Optimization Letters 13, p. 847–856
    // (also at https://link.springer.com/article/10.1007/s11590-018-1284-4)
    particle_positions[0] = x;
    for (i = 1; i < nm_particles; i++) {
      scalar_t x_inf_norm = max_abs_vec(x);
      // if smaller than 1, set to 1
      scalar_t a = x_inf_norm < 1.0 ? 1.0 : x_inf_norm;
      // if larger than 10, set to 10
      scalar_t scale = a < 10 ? a : 10;
      for (i = 1; i < nm_particles; i++) {
        particle_positions[i] = x;
        particle_positions[i][i - 1] = x[i - 1] + scale;
      }
      // update first simplex point
      auto n = static_cast<scalar_t>(x.size());
      for (i = 0; i < x.size(); i++) {
        particle_positions[0][i] = x[i] + ((1.0 - sqrt(n + 1.0)) / n * scale);
      }
    }
    // the rest according to PSO
    for (i = nm_particles; i < n_particles; i++) {
      for (size_t j = 0; j < n_dim; j++) {
        // update velocities and positions
        scalar_t temp = std::abs(upper[j] - lower[j]);
        this->particle_positions[i][j] =
            lower[j] + ((upper[j] - lower[j]) * generator());
        this->particle_velocities[i][j] = -temp + (generator() * 2.0 * temp);
      }
    }
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
    for (i = 0; i < n_particles; i++) {
      this->particle_current_values[i] =
          f_multiplier * f(particle_positions[i]);
      this->function_calls_used++;
    }
  }
  template <const bool minimize, const bool bound>
  void apply_simplex(std::vector<scalar_t> &centroid,
                     std::vector<scalar_t> &temp_reflect,
                     std::vector<scalar_t> &temp_expand,
                     std::vector<scalar_t> &temp_contract,
                     std::vector<size_t> &current_order,
                     const size_t nm_particles, const size_t n_dim,
                     const std::vector<scalar_t> &upper,
                     const std::vector<scalar_t> &lower) {
    const scalar_t best_score = this->particle_current_values[current_order[0]];
    const size_t worst_id = current_order[nm_particles - 1];
    const size_t second_worst_id = current_order[nm_particles - 2];
    // update centroid of all points except for the worst one
    this->update_centroid(centroid, current_order, nm_particles - 1);
    // reflect worst point
    simplex_transform<scalar_t, true, bound>(this->particle_positions[worst_id],
                                             centroid, temp_reflect,
                                             this->alpha, upper, lower);
    // set constant multiplier for minimization
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
    // score reflected point
    const scalar_t ref_score = f_multiplier * f(temp_reflect);
    this->function_calls_used++;
    // if reflected point is better than second worst, but not better than best
    if (ref_score >= best_score &&
        ref_score < this->particle_current_values[second_worst_id]) {
      this->particle_positions[worst_id] = temp_reflect;
      this->particle_current_values[worst_id] = ref_score;
      // otherwise if this is the best score so far, expand
    } else if (ref_score < best_score) {
      simplex_transform<scalar_t, false, bound>(
          temp_reflect, centroid, temp_expand, this->gamma, upper, lower);
      // obtain score for expanded point
      const scalar_t exp_score = f_multiplier * f(temp_expand);
      this->function_calls_used++;
      // if this is better than the expanded point score, replace the worst
      // point with the expanded point, otherwise replace it with
      // the reflected point
      std::vector<scalar_t> &replacement =
          exp_score < ref_score ? temp_expand : temp_reflect;
      this->particle_positions[worst_id] = replacement;
      this->particle_current_values[worst_id] =
          exp_score < ref_score ? exp_score : ref_score;
      // otherwise we have a point  worse than the 'second worst'
    } else {
      // contract outside - here we overwrite the 'temp_expand' and it
      // functionally becomes 'temp_contract'
      const scalar_t worst_score = this->particle_current_values[worst_id];
      simplex_transform<scalar_t, false, bound>(
          ref_score < worst_score
              ? temp_reflect
              :
              // or point is the worst point so far - contract inside
              this->particle_positions[worst_id],
          centroid, temp_contract, this->rho, upper, lower);
      const scalar_t cont_score = f_multiplier * f(temp_contract);
      this->function_calls_used++;
      // if this contraction is better than the reflected point or worst point
      if (cont_score < std::min(ref_score, worst_score)) {
        // replace worst point with contracted point
        this->particle_positions[worst_id] = temp_contract;
        this->particle_current_values[worst_id] = cont_score;
        // otherwise shrink
      } else {
        // if we had not violated the bounds before shrinking, shrinking
        // will not cause new violations - hence no bounds applied here
        shrink(current_order, this->sigma, nm_particles, n_dim);
        // only in this case do we have to score again
        for (size_t i = 1; i < nm_particles; i++) {
          this->particle_current_values[current_order[i]] =
              f_multiplier * f(this->particle_positions[current_order[i]]);
        }
        this->function_calls_used += nm_particles - 1;
        // re-sort, as PSO step is called after this
        std::sort(current_order.begin(), current_order.end(),
                  [&](size_t left_id, size_t right_id) -> bool {
                    return this->particle_current_values[left_id] <
                           this->particle_current_values[right_id];
                  });
      }
    }
  }
  template <const bool minimize, const bool bound>
  void apply_pso(const std::vector<size_t> &current_order,
                 const size_t n_simplex_particles,
                 const size_t n_total_particles, const size_t n_dim,
                 const std::vector<scalar_t> &upper,
                 const std::vector<scalar_t> &lower) {
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
    bool order_flip = false;
    size_t best_in_pair = current_order[n_simplex_particles];
    const std::vector<scalar_t> &best =
        this->particle_positions[current_order[0]];
    for (size_t i = n_simplex_particles; i < n_total_particles; i++) {
      const size_t id = current_order[i];
      if (order_flip) {
        best_in_pair = current_order[i + 1];
      }
      order_flip = static_cast<bool>((i - n_simplex_particles) % 2);
      // get references to current particle, current velocity and pairwise best
      // particle
      // these MUST be references (the `&` previously bound only to `particle`,
      // so velocity/pairwise_best were copies and velocity updates were lost)
      std::vector<scalar_t> &particle = particle_positions[id];
      std::vector<scalar_t> &velocity = this->particle_velocities[id];
      const std::vector<scalar_t> &pairwise_best =
          particle_positions[best_in_pair];
      for (size_t j = 0; j < n_dim; j++) {
        // generate random movements
        const scalar_t r_p = generator(), r_g = generator();
        // update current velocity for current particle - inertia update
        // TODO(JSzitas): SIMD Candidate
        const scalar_t temp =
            (this->inertia * velocity[j]) +
            // cognitive update - this should be based on better
            // particle of each 2 particle pairs
            this->cognitive_coef * r_p * (pairwise_best[j] - particle[j]) +
            // social update (moving more if further away from 'best' position)
            this->social_coef * r_g * (best[j] - particle[j]);
        velocity[j] = temp;
        particle[j] += temp;
        // clamp the resulting *position* (not the velocity) into the box
        if constexpr (bound) {
          particle[j] = std::clamp(particle[j], lower[j], upper[j]);
        }
      }
      // rerun function evaluation
      this->particle_current_values[id] = f_multiplier * f(particle);
      this->function_calls_used++;
    }
  }
  void update_centroid(std::vector<scalar_t> &centroid,
                       const std::vector<size_t> &current_order,
                       const size_t last_point) {
    // reset centroid - fill with 0
    std::fill(centroid.begin(), centroid.end(), 0.0);
    // iterate through 0 to last_point - 1 - last point taken to be the Nth
    // best point in an N+1 dimensional simplex
    size_t i = 0;
    for (; i < last_point; i++) {
      const std::vector<scalar_t> &particle =
          this->particle_positions[current_order[i]];
      // TODO(JSzitas): SIMD Candidate
      for (size_t j = 0; j < centroid.size(); j++) {
        centroid[j] += particle[j];
      }
    }
    for (auto &val : centroid) val /= static_cast<scalar_t>(i);
  }
  void shrink(const std::vector<size_t> &current_order, const scalar_t sigma_,
              const size_t nm_particles, const size_t n_dim) {
    // take a reference to the best vector
    const std::vector<scalar_t> &best =
        this->particle_positions[current_order[0]];
    for (size_t i = 1; i < nm_particles; i++) {
      // update all items in current vector using the best vector -
      // hopefully the contiguous data here can help a bit with cache
      // locality
      std::vector<scalar_t> &current =
          this->particle_positions[current_order[i]];
      // TODO(JSzitas): SIMD Candidate
      for (size_t j = 0; j < n_dim; j++) {
        current[j] = best[j] + sigma_ * (current[j] - best[j]);
      }
    }
  }
  scalar_t simplex_std_err(const std::vector<size_t> &current_order,
                           const size_t nm_particles) {
    size_t i = 0;
    scalar_t mean_val = 0, result = 0;
    for (; i < nm_particles; i++) {
      mean_val += this->particle_current_values[current_order[i]];
    }
    mean_val /= static_cast<scalar_t>(i);
    i = 0;
    for (; i < nm_particles; i++) {
      result +=
          pow(this->particle_current_values[current_order[i]] - mean_val, 2);
    }
    result /= (scalar_t)(i - 1);
    return sqrt(result);
  }
};
}  // namespace nlsolver

namespace nlsolver::rootfinder {
template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> bisection(
    Callable &f, scalar_t &x, const scalar_t lower = -100,
    const scalar_t upper = 100, const scalar_t eps = 1e-6,
    const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper;
  if (a > b) {
    std::swap(a, b);
  }
  if (f_lam(a) * f_lam(b) >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  scalar_t mid, val;
  size_t iter = 0;
  while (true) {
    // define midpoint
    mid = (a + b) / 2;
    val = f_lam(mid);
    if (std::abs(val) < eps || iter > max_iter) {
      //; solution found
      x = mid;
      return solver_status<scalar_t>(val, iter, f_evals_used);
    }
    if (val > 0) {
      b = mid;
    } else if (val < 0) {
      a = mid;
    }
    iter++;
  }
}
// inspired by this talk: https://www.youtube.com/watch?v=J48YTbdJNNc
// I hope Andrei won't hold a grudge; probably proportional to how bad the
// method is?
template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> alexandrescu_bisection(
    Callable &f, scalar_t &x, const scalar_t lower = -100,
    const scalar_t upper = 100, const scalar_t eps = 1e-6,
    const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper;
  if (a > b) {
    std::swap(a, b);
  }
  if (f_lam(a) * f_lam(b) >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  scalar_t mid, val;
  // leapfrogging idea
  mid = (a + b) / 2;
  val = f_lam(mid);
  if (std::abs(val) < eps) {
    // solution found
    x = mid;
    return solver_status<scalar_t>(val, 0, f_evals_used);
  }
  if (val > 0) {
    // attempt to leap-frog -
    auto temp = mid - (mid - a) / 4;
    // still bracketing
    if (f_lam(temp) > 0) {
      b = temp;
    } else {
      a = temp;
    }
  } else if (val < 0) {
    // opposite side leapfrog
    auto temp = mid + (b - mid) / 4;
    // still bracketing
    if (f_lam(temp) < 0) {
      a = temp;
    } else {
      b = temp;
    }
  }
  size_t iter = 1;
  // standard binary search
  while (true) {
    mid = (a + b) / 2;
    val = f_lam(mid);
    if (std::abs(val) < eps || iter > max_iter) {
      //; solution found
      x = mid;
      return solver_status<scalar_t>(val, iter, f_evals_used);
    }
    if (val > 0) {
      b = mid;
    } else if (val < 0) {
      a = mid;
    }
    iter++;
  }
}
template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> false_position(
    Callable &f, scalar_t &x, const scalar_t lower = -100,
    const scalar_t upper = 100, const scalar_t eps = 1e-6,
    const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t val_a = f_lam(lower), val_b = f_lam(upper);
  if (val_a * val_b >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  scalar_t a = lower, b = upper, mid, val;
  size_t iter = 0;
  while (true) {
    // define midpoint
    mid = a + ((b - a) * val_a) / (val_a - val_b);
    val = f_lam(mid);
    if (std::abs(val) < eps || iter > max_iter) {
      //; solution found
      x = mid;
      return solver_status<scalar_t>(val, iter, f_evals_used);
    }
    if (val < 0) {
      a = mid;
      val_a = val;
    } else if (val > 0) {
      b = mid;
      val_b = val;
    }
    iter++;
  }
}

template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> brent(Callable &f, scalar_t &x,
                                               const scalar_t lower,
                                               const scalar_t upper,
                                               const scalar_t tol = 1e-12,
                                               const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  // functor
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper;
  scalar_t val_a = f_lam(a), val_b = f_lam(b);
  if (val_a * val_b >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  scalar_t c = a, val_c = val_a, s, val_s, d;
  size_t iter = 0;
  bool flag = true;
  while (true) {
    // inverse quadratic interpolation
    if ((val_a != val_c) && (val_b != val_c)) {
      s = ((a * val_b * val_c) / ((val_a - val_b) * (val_a - val_c))) +
          ((b * val_a * val_c) / ((val_b - val_a) * (val_b - val_c))) +
          ((c * val_a * val_b) / ((val_c - val_a) * (val_c - val_b)));
    } else {
      // otherwise do the secant method
      s = b - val_b * ((b - a) / (val_b - val_a));
    }
    // bisection
    if (!(s > (((3 * a) + b) / 4) && s < b) ||
        (flag && (std::abs(s - b) >= (std::abs(b - c) / 2))) ||
        (!flag && std::abs(s - b) >= (std::abs(c - d) / 2)) ||
        (flag && std::abs(b - c) < tol) || (!flag && std::abs(c - d) < tol)) {
      s = (a + b) / 2;
      flag = true;
    } else {
      flag = false;
    }
    val_s = f_lam(s);
    d = c;
    c = b;
    val_c = val_b;
    if ((val_a * val_s) < 0) {
      b = s;
      val_b = val_s;
    } else {
      a = s;
      val_a = val_s;
    }
    if (std::abs(val_a) < std::abs(val_b)) {
      std::swap(a, b);
      std::swap(val_a, val_b);
    }
    if (std::abs(val_b) < tol || std::abs(val_s) < tol ||
        std::abs(b - a) < tol || iter >= max_iter) {
      x = b;
      return solver_status<scalar_t>(val_b, iter, f_evals_used);
    }
    iter++;
  }
}

template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> ridders(Callable &f, scalar_t &x,
                                                 const scalar_t lower,
                                                 const scalar_t upper,
                                                 const scalar_t tol = 1e-12,
                                                 const scalar_t eps = 1e-12,
                                                 const size_t max_iter = 5) {
  size_t f_evals_used = 0;
  // functor
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper;
  scalar_t val_a = f_lam(a), val_b = f_lam(b);
  if (val_a * val_b >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  size_t iter = 0;
  scalar_t mid, val_mid, new_mid;
  while (true) {
    // form an evaluate at midpoint
    // evaluate at midpoint
    mid = (a + b) / 2, val_mid = f_lam(mid);
    // new iterate
    new_mid =
        mid + (mid - a) * (std::copysign(1., val_a - val_b) * val_mid /
                           std::sqrt(std::pow(val_mid, 2) - (val_a * val_b)));
    // check if equal zero, if yes return
    scalar_t val_new_mid = f_lam(new_mid);
    // check tolerances
    if (std::min(std::abs(new_mid - a), std::abs(new_mid - b)) < tol ||
        std::abs(val_new_mid) < eps || iter >= max_iter) {
      x = new_mid;
      return solver_status<scalar_t>(val_new_mid, iter, f_evals_used);
    }
    // not matching signs
    if ((val_mid * val_new_mid) < 0) {
      a = mid;
      val_a = val_mid;
      b = new_mid;
      val_b = val_new_mid;
      // likewise, not matching signs
    } else if ((val_a * val_new_mid) < 0) {
      a = new_mid;
      val_a = val_new_mid;
    } else {
      b = new_mid;
      val_b = val_new_mid;
    }
    iter++;
  }
}

namespace nlsolver::internal::circulant {
template <typename T, const size_t size>
struct fixed_circulant {
 private:
  std::array<T, size> data = std::array<T, size>();
  size_t circle_index = 0;

 public:
  void load(std::array<T, size> &&x) { data = x; }
  T &operator[](const size_t i) {
    return this->data[(this->circle_index + i) % size];
  }
  void push_back(const T item) {
    this->data[this->circle_index] = item;
    this->circle_index = (this->circle_index + 1) % size;
  }
  auto last() { return this->data[this->circle_index]; }
};
};  // namespace nlsolver::internal::circulant
template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> tiruneh(
    Callable &f, scalar_t &x,
    const std::array<scalar_t, 3> x_k = {-100., 0., 100.},
    const scalar_t eps = 1e-6, const scalar_t tol = 1e-12,
    const size_t max_iter = 10) {
  size_t f_evals_used = 0;
  // functor
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  // initialize circulant
  nlsolver::internal::circulant::fixed_circulant<scalar_t, 3> k, f_k;
  for (auto &val : x_k) {
    k.push_back(val);
  }
  for (size_t i = 0; i < 3; i++) {
    f_k.push_back(f_lam(k[i]));
  }
  size_t iter = 0;
  scalar_t temp = 0.;
  while (true) {
    if (std::abs(f_k.last()) < tol || iter > max_iter ||
        std::abs(f_k[0] - f_k[1]) < eps) {
      x = k.last();
      return solver_status<scalar_t>(f_k.last(), iter, f_evals_used);
    }
    // iteration and replacement
    temp = k[2] - (f_k[2] * (f_k[0] - f_k[1])) /
                      (((f_k[0] - f_k[2]) / (k[0] - k[2]) * (f_k[0] - f_k[1])) -
                       f_k[0] * (((f_k[0] - f_k[2]) / (k[0] - k[2])) -
                                 ((f_k[1] - f_k[2]) / (k[1] - k[2]))));
    // update
    k.push_back(temp);
    f_k.push_back(f_lam(temp));
    iter++;
  }
}
template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> itp(
    Callable &f, scalar_t &x, const scalar_t lower, const scalar_t upper,
    const scalar_t kappa1 = 0.3, const scalar_t kappa2 = 2.1,
    const scalar_t n0 = 1, const scalar_t tol = 1e-12,
    const scalar_t eps = 1e-12, const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  // functor
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper;
  // standard bracketing condition
  scalar_t val_a = f_lam(a), val_b = f_lam(b);
  if (val_a * val_b >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  size_t iter = 0;
  const scalar_t two_eps = 2 * eps;
  const scalar_t n_max = std::log2((b - a) / two_eps) + n0;
  scalar_t mid, r, delta, b_min_a, interp, sigma, temp, val_itp = 100000.0;
  while (true) {
    b_min_a = b - a;
    if (b_min_a < two_eps || iter >= max_iter) {
      x = (a + b) / 2.;
      return solver_status<scalar_t>(val_itp, iter, f_evals_used);
    }
    mid = (a + b) / 2;
    r = eps * std::pow(2, n_max - 1) - (b_min_a / 2);
    delta = kappa1 * std::pow(b_min_a, kappa2);
    // interpolate
    interp = (val_b * a - val_a * b) / (val_b - val_a);
    // truncate
    temp = mid - interp;
    sigma = (temp) > 0;
    const bool project = temp <= r;
    if (delta <= std::abs(temp)) {
      interp += sigma * delta;
    } else {
      interp = mid;
    }
    // project
    if (project) {
      temp = interp;
    } else {
      temp = mid - sigma * r;
    }
    // update
    val_itp = f_lam(temp);
    if (val_itp > 0) {
      b = temp;
      val_b = val_itp;
    } else if (val_itp < 0) {
      a = temp;
      val_a = val_itp;
    } else {
      x = temp;
      return solver_status<scalar_t>(val_itp, iter, f_evals_used);
    }
    iter++;
  }
}

template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> chandrupatla(
    Callable &f, scalar_t &x, const scalar_t lower, const scalar_t upper,
    const scalar_t eps_m = 1e-10, const scalar_t eps_a = 2e-10,
    const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  // functor
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper, c = upper;
  // standard bracketing condition
  scalar_t val_a = f_lam(a), val_b = f_lam(b), val_c = val_b;
  if (val_a * val_b >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  size_t iter = 0;
  scalar_t t = 0.5, x_t, val_t, x_m, val_m;
  while (true) {
    // use t to linearly interpolate between a and b,
    // and evaluate this function as our newest estimate xt
    x_t = b + t * (a - b);
    val_t = f_lam(x_t);
    //  not matching signs
    if ((val_t * val_b) < 0) {
      c = a;
      a = b;
      val_c = val_a;
      val_a = val_b;
    } else {
      c = b;
      val_c = val_b;
    }
    b = x_t;
    val_b = val_t;
    // set xm so that f(xm) is the minimum magnitude of f(a) and f(b)
    const bool b_smaller_a = std::abs(val_b) < std::abs(val_a);
    x_m = b_smaller_a ? b : a;
    val_m = b_smaller_a ? val_b : val_a;
    if (std::abs(val_m) < eps_a || iter > max_iter) {
      x = x_m;
      return solver_status<scalar_t>(val_m, iter, f_evals_used);
    }
    // Figure out values xi and phi to determine which method we should use next
    const scalar_t tol = 2 * eps_m * std::abs(x_m) + eps_a;
    const scalar_t t_lim = tol / std::abs(a - c);
    if (t_lim > 0.5) {
      x = x_m;
      return solver_status<scalar_t>(val_m, iter, f_evals_used);
    }
    const scalar_t xi = (a - b) / (c - b);
    const scalar_t phi = (val_a - val_b) / (val_c - val_b);
    // inverse quadratic extrapolation bit
    if ((std::pow(phi, 2) < xi) && (std::pow(1 - phi, 2) < (1 - xi))) {
      t = val_a / (val_b - val_a) * val_c / (val_b - val_c) +
          (c - a) / (b - a) * val_a / (val_c - val_a) * val_b / (val_c - val_b);
    } else {
      t = 0.5;
    }
    // limit range
    t = std::min(1. - t_lim, std::max(t_lim, t));
    iter++;
  }
}
template <typename Callable, typename scalar_t>
[[maybe_unused]] solver_status<scalar_t> hybrid_composite(
    Callable &f, scalar_t &x, const scalar_t lower = -100,
    const scalar_t upper = 100, const scalar_t eps = 1e-6,
    const scalar_t tol = 1e-10, const scalar_t min_improvement = 1e-2,
    const size_t max_noleaps = 3, const size_t max_iter = 200) {
  size_t f_evals_used = 0;
  auto f_lam = [&](auto x) {
    auto val = f(x);
    f_evals_used++;
    return val;
  };
  scalar_t a = lower, b = upper;
  if (a > b) {
    std::swap(a, b);
  }
  if (f_lam(a) * f_lam(b) >= 0) {
    // initial interval does not bracket a root: signal failure (NaN value +
    // success()==false) instead of printing and returning a spurious 0
    return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                   0, f_evals_used, 0, 0, false);
  }
  scalar_t mid, val, last_val = 0.;
  size_t iter = 0, noleap = 0;
  while (true) {
    // leapfrogging idea - basically a somewhat faster binary search (sometimes)
    // followed by tiruneh for quicker refinement
    mid = (a + b) / 2;
    val = f_lam(mid);
    if (std::abs(val - last_val) < min_improvement) break;
    if (std::abs(val) < eps || iter >= max_iter) {
      // solution found
      x = mid;
      return solver_status<scalar_t>(val, iter, f_evals_used);
    }
    if (val > 0) {
      // attempt to leap-frog -
      auto temp = mid - (mid - a) / 4;
      // still bracketing - the opposite is actually quite advantageous for us
      if (f_lam(temp) > 0) {
        noleap++;
        b = temp;
      } else {
        a = temp;
      }
    } else if (val < 0) {
      // opposite side leapfrog
      auto temp = mid + (b - mid) / 4;
      // still bracketing
      if (f_lam(temp) < 0) {
        noleap++;
        a = temp;
      } else {
        b = temp;
      }
    }
    iter++;
    last_val = val;
    if (noleap > max_noleaps) break;
  }
  // if we got here it means we reached the maximum number of 'noleaps', or the
  // improvement was too small; time has come to use tiruneh
  // initialize circulant
  nlsolver::internal::circulant::fixed_circulant<scalar_t, 3> k, f_k;
  k.load(std::array<scalar_t, 3>{a, b, mid});
  f_k.load(std::array<scalar_t, 3>{f(a), f(b), val});
  f_evals_used += 2;
  scalar_t temp = 0.;
  while (true) {
    if (std::abs(f_k.last()) < tol || iter > max_iter ||
        std::abs(f_k[0] - f_k[1]) < eps) {
      x = k.last();
      return solver_status<scalar_t>(f_k.last(), iter, f_evals_used);
    }
    // iteration and replacement
    temp = k[2] - (f_k[2] * (f_k[0] - f_k[1])) /
                      (((f_k[0] - f_k[2]) / (k[0] - k[2]) * (f_k[0] - f_k[1])) -
                       f_k[0] * (((f_k[0] - f_k[2]) / (k[0] - k[2])) -
                                 ((f_k[1] - f_k[2]) / (k[1] - k[2]))));
    // update
    k.push_back(temp);
    f_k.push_back(f_lam(temp));
    iter++;
  }
}
};  // namespace nlsolver::rootfinder

namespace nlsolver {

// indices of the k smallest elements of v, ascending by value
template <typename scalar_t>
inline std::vector<size_t> index_partial_sort(const std::vector<scalar_t> &v,
                                              const size_t k) {
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i < v.size(); i++) idx[i] = i;
  const size_t kk = std::min(k, v.size());
  std::partial_sort(
      idx.begin(), idx.begin() + kk, idx.end(),
      [&v](const size_t a, const size_t b) { return v[a] < v[b]; });
  idx.resize(kk);
  return idx;
}

// Search distribution of a CMA-type method at exit: mean, global step sigma
// and covariance C (row-major n x n). On a smooth objective sigma^2 C converges
// to a multiple of the inverse Hessian, which is what Composite seeds L-BFGS-B
// with.
template <typename scalar_t>
struct search_distribution {
  std::vector<scalar_t> mean;
  scalar_t sigma;
  std::vector<scalar_t> C;
};

// (mu/mu_w, lambda)-CMA-ES, following Hansen, "The CMA Evolution Strategy: A
// Tutorial" (https://arxiv.org/abs/1604.00772). Operates entirely on flat
// std::vector storage; the covariance eigendecomposition is delegated to
// tinyqr::QRSolver (whose eigenvectors come back as the *rows* of the returned
// matrix, hence the transpose into the column-eigenvector matrix Bmat below).
template <typename Callable, typename RNG, typename scalar_t = double>
class [[maybe_unused]] CMAES {
 public:
  using distribution = search_distribution<scalar_t>;
  // default population size 4 + floor(3 ln n), at least 5 (Hansen's rule)
  static size_t default_lambda(const size_t n) {
    const scalar_t nn = static_cast<scalar_t>(n);
    const auto lambda = static_cast<size_t>(4 + std::floor(3.0 * std::log(nn)));
    return std::max<size_t>(lambda, 5);
  }
  // Hansen's default strategy parameters for dimension n and a population
  // multiplier: population size, parent number, recombination weights, and
  // the learning rates of the two evolution paths, the rank-one and rank-mu
  // covariance updates, the step-size damping and the expected norm of a
  // standard normal vector
  struct parameters {
    size_t lambda, mu;
    std::vector<scalar_t> w;
    scalar_t mu_eff, cc, cs, c1, cmu, damps, chi_n;
  };
  static parameters strategy_parameters(const size_t n,
                                        const size_t pop_mult = 1) {
    const scalar_t nn = static_cast<scalar_t>(n);
    parameters p;
    p.lambda = default_lambda(n) * pop_mult;
    p.mu = p.lambda / 2;
    p.w.resize(p.mu);
    for (size_t i = 0; i < p.mu; i++) {
      p.w[i] = std::log(static_cast<scalar_t>(p.mu) + 0.5) -
               std::log(static_cast<scalar_t>(i) + 1.0);
    }
    scalar_t w_sum = 0.0;
    for (const auto wi : p.w) w_sum += wi;
    for (auto &wi : p.w) wi /= w_sum;  // normalize so sum(w) == 1
    p.mu_eff = 0.0;
    for (const auto wi : p.w) p.mu_eff += wi * wi;
    p.mu_eff = 1.0 / p.mu_eff;
    p.cc = (4.0 + p.mu_eff / nn) / (nn + 4.0 + 2.0 * p.mu_eff / nn);
    p.cs = (p.mu_eff + 2.0) / (nn + p.mu_eff + 5.0);
    p.c1 = 2.0 / (std::pow(nn + 1.3, 2.0) + p.mu_eff);
    p.cmu = std::min(1.0 - p.c1, 2.0 * (p.mu_eff - 2.0 + 1.0 / p.mu_eff) /
                                     (std::pow(nn + 2.0, 2.0) + p.mu_eff));
    p.damps =
        1.0 + p.cs +
        2.0 * std::max(0.0, std::sqrt((p.mu_eff - 1.0) / (nn + 1.0)) - 1.0);
    p.chi_n = std::sqrt(nn) * (1.0 - 1.0 / (4.0 * nn) + 1.0 / (21.0 * nn * nn));
    return p;
  }
  // generations over which the covariance update forgets its past, 1 / (c1 +
  // cmu); a run shorter than a few of these has not learned the problem's shape
  static scalar_t covariance_horizon(const size_t n,
                                     const size_t pop_mult = 1) {
    const parameters p = strategy_parameters(n, pop_mult);
    return 1.0 / (p.c1 + p.cmu);
  }

 private:
  Callable &f;
  RNG &generator;
  scalar_t m_step;
  const size_t max_iter;
  const scalar_t condition, x_delta, f_delta;
  // population size multiplier; IPOP restarts double it on every restart
  const size_t pop_mult;
  distribution exit_distribution_;

 public:
  // constructor
  [[maybe_unused]] CMAES<Callable, RNG, scalar_t>(
      Callable &f, RNG &generator, const scalar_t m_step = 0.5,
      const size_t max_iter = 1000, const scalar_t condition = 1e14,
      const scalar_t x_delta = 1e-12, const scalar_t f_delta = 1e-12,
      const size_t pop_mult = 1)
      : f(f),
        generator(generator),
        m_step(m_step),
        max_iter(max_iter),
        condition(condition),
        x_delta(x_delta),
        f_delta(f_delta),
        pop_mult(pop_mult) {}
  // the search distribution when the last minimize / maximize call returned
  [[maybe_unused]] const distribution &exit_distribution() const {
    return exit_distribution_;
  }
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;  // empty bounds => unconstrained
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;  // empty bounds => unconstrained
    return this->solve<false, false>(x, lower, upper);
  }
  // minimize with box constraints; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  // maximize with box constraints; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    const size_t n = x.size();
    // we always minimize internally; for maximization we minimize -f
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const scalar_t nn = static_cast<scalar_t>(n);

    // ---- selection and strategy parameters ----
    const parameters par = strategy_parameters(n, this->pop_mult);
    const size_t lambda = par.lambda, mu = par.mu;
    const std::vector<scalar_t> &w = par.w;
    const scalar_t mu_eff = par.mu_eff, cc = par.cc, cs = par.cs, c1 = par.c1,
                   cmu = par.cmu, damps = par.damps, chiN = par.chi_n;

    // ---- dynamic state ----
    scalar_t sigma = m_step;
    std::vector<scalar_t> mean = x;
    // project the initial mean into the feasible box, if constrained
    if constexpr (constrained) {
      for (size_t i = 0; i < n; i++) {
        mean[i] = mean[i] < lower[i]   ? lower[i]
                  : mean[i] > upper[i] ? upper[i]
                                       : mean[i];
      }
    }
    std::vector<scalar_t> pc(n, 0.0), ps(n, 0.0);
    // covariance C, its eigenvectors Bmat (column j == eigenvector j) and the
    // square roots of its eigenvalues dvec; initialized to the identity.
    std::vector<scalar_t> C(n * n, 0.0), Bmat(n * n, 0.0), dvec(n, 1.0);
    for (size_t i = 0; i < n; i++) {
      C[i * n + i] = 1.0;
      Bmat[i * n + i] = 1.0;
    }
    tinyqr::QRSolver<scalar_t> qr(n);
    // refresh the eigensystem about every 1/(10 (c1+cmu)) generations so the
    // O(n^3) decomposition stays amortized (the extra n factor here drove the
    // ratio below 1, i.e. it decomposed every single generation)
    const size_t eigen_every =
        std::max<size_t>(1, static_cast<size_t>(1.0 / (10.0 * (c1 + cmu))));
    size_t last_eigen = 0;

    // per-generation sample storage (flattened, particle-major)
    std::vector<scalar_t> Z(n * lambda, 0.0), Y(n * lambda, 0.0),
        Xp(n * lambda, 0.0);
    std::vector<scalar_t> costs(lambda, 0.0);
    std::vector<scalar_t> y_w(n, 0.0), Bty(n, 0.0), tmp(n, 0.0),
        eval_buf(n, 0.0);

    // best sampled point so far - used only for the reported result
    std::vector<scalar_t> best_x = mean;
    scalar_t best_f = sign * this->f(mean);
    size_t f_evals = 1;

    size_t iter = 0;
    while (true) {
      // ---- sample: x_k = mean + sigma * B * (d .* z_k), z_k ~ N(0, I) ----
      rnorm<scalar_t>(this->generator, Z.begin(), Z.end());
      for (size_t k = 0; k < lambda; k++) {
        for (size_t i = 0; i < n; i++) {
          scalar_t yi = 0.0;
          for (size_t j = 0; j < n; j++) {
            yi += Bmat[i * n + j] * dvec[j] * Z[k * n + j];
          }
          scalar_t xi = mean[i] + sigma * yi;
          // project candidate into the feasible box and keep y == (x-m)/sigma
          // consistent so the mean / path / covariance updates stay coherent
          if constexpr (constrained) {
            xi = xi < lower[i] ? lower[i] : xi > upper[i] ? upper[i] : xi;
            yi = (xi - mean[i]) / sigma;
          }
          Y[k * n + i] = yi;
          Xp[k * n + i] = xi;
          eval_buf[i] = xi;
        }
        const scalar_t fx = this->f(eval_buf);
        costs[k] = sign * fx;
        f_evals++;
        if (sign * fx < best_f) {
          best_f = sign * fx;
          best_x = eval_buf;
        }
      }
      // ---- selection: the mu best (smallest internal cost) ----
      const std::vector<size_t> idx = index_partial_sort(costs, mu);
      // ---- recompute mean; accumulate weighted step y_w = sum_i w_i y_{i:l}
      // --
      const std::vector<scalar_t> old_mean = mean;
      for (size_t i = 0; i < n; i++) {
        scalar_t yw = 0.0;
        for (size_t r = 0; r < mu; r++) yw += w[r] * Y[idx[r] * n + i];
        y_w[i] = yw;
        mean[i] = old_mean[i] + sigma * yw;
      }
      // ---- step-size path: ps = (1-cs)ps + sqrt(cs(2-cs)mu_eff) C^-1/2 y_w
      // ---- C^-1/2 y_w = B diag(1/d) B^T y_w
      for (size_t j = 0; j < n; j++) {
        scalar_t s = 0.0;
        for (size_t i = 0; i < n; i++) s += Bmat[i * n + j] * y_w[i];
        Bty[j] = s / dvec[j];
      }
      for (size_t i = 0; i < n; i++) {
        scalar_t s = 0.0;
        for (size_t j = 0; j < n; j++) s += Bmat[i * n + j] * Bty[j];
        tmp[i] = s;
      }
      const scalar_t ps_coef = std::sqrt(cs * (2.0 - cs) * mu_eff);
      for (size_t i = 0; i < n; i++) {
        ps[i] = (1.0 - cs) * ps[i] + ps_coef * tmp[i];
      }
      scalar_t ps_norm = 0.0;
      for (size_t i = 0; i < n; i++) ps_norm += ps[i] * ps[i];
      ps_norm = std::sqrt(ps_norm);
      // ---- Heaviside step ----
      const scalar_t hsig_den = std::sqrt(
          1.0 - std::pow(1.0 - cs, 2.0 * static_cast<scalar_t>(iter + 1)));
      const bool hsig = (ps_norm / hsig_den / chiN) < (1.4 + 2.0 / (nn + 1.0));
      // ---- covariance path: pc = (1-cc)pc + hsig sqrt(cc(2-cc)mu_eff) y_w
      // ----
      const scalar_t pc_coef = std::sqrt(cc * (2.0 - cc) * mu_eff);
      for (size_t i = 0; i < n; i++) {
        pc[i] = (1.0 - cc) * pc[i] + (hsig ? pc_coef : 0.0) * y_w[i];
      }
      // ---- covariance update (rank-1 + rank-mu, with hsig variance loss) ----
      const scalar_t delta_hsig = (hsig ? 0.0 : cc * (2.0 - cc));
      const scalar_t c_decay = 1.0 - c1 - cmu + c1 * delta_hsig;
      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
          scalar_t rank_mu = 0.0;
          for (size_t r = 0; r < mu; r++) {
            rank_mu += w[r] * Y[idx[r] * n + i] * Y[idx[r] * n + j];
          }
          C[i * n + j] =
              c_decay * C[i * n + j] + c1 * pc[i] * pc[j] + cmu * rank_mu;
        }
      }
      // ---- step-size update ----
      sigma *= std::exp((cs / damps) * (ps_norm / chiN - 1.0));
      // ---- eigendecomposition of C on schedule ----
      if (iter - last_eigen >= eigen_every) {
        last_eigen = iter;
        // enforce exact symmetry before decomposing
        for (size_t i = 0; i < n; i++) {
          for (size_t j = i + 1; j < n; j++) {
            const scalar_t avg = 0.5 * (C[i * n + j] + C[j * n + i]);
            C[i * n + j] = avg;
            C[j * n + i] = avg;
          }
        }
        qr.solve(C, 100, 1e-12);
        const auto &ev = qr.eigenvalues();
        const auto &QQ = qr.eigenvectors();  // eigenvectors are ROWS of QQ
        for (size_t j = 0; j < n; j++) {
          const scalar_t lam = ev[j] > 1e-30 ? ev[j] : 1e-30;
          dvec[j] = std::sqrt(lam);
          for (size_t i = 0; i < n; i++) Bmat[i * n + j] = QQ[j * n + i];
        }
      }
      // ---- termination checks ----
      scalar_t dmax = dvec[0], dmin = dvec[0];
      for (size_t j = 1; j < n; j++) {
        dmax = std::max(dmax, dvec[j]);
        dmin = std::min(dmin, dvec[j]);
      }
      const scalar_t cond = (dmax * dmax) / (dmin * dmin);
      scalar_t cost_spread = 0.0;
      for (size_t r = 1; r < mu; r++) {
        cost_spread =
            std::max(cost_spread, std::abs(costs[idx[r]] - costs[idx[0]]));
      }
      const scalar_t x_step = sigma * dmax;
      iter++;
      if (iter >= this->max_iter || cond > this->condition ||
          cost_spread < this->f_delta || x_step < this->x_delta ||
          !std::isfinite(sigma)) {
        x = best_x;
        m_step = sigma;
        exit_distribution_ = distribution{mean, sigma, C};
        return solver_status<scalar_t>(sign * best_f, iter, f_evals);
      }
    }
  }
};

// (1+1)-CMA-ES: one candidate per evaluation, accepted when it does not
// worsen the incumbent. The step size follows the smoothed success rate and
// the covariance takes a rank-one update from each successful step, applied
// directly to its Cholesky factor A (C = A A'), so no decomposition is ever
// needed. [[CITATION]] Igel, Suttorp, Hansen, "A computational efficient
// covariance matrix update and a (1+1)-CMA for evolution strategies", GECCO
// 2006, and Suttorp, Hansen, Igel, "Efficient covariance matrix update for
// variable metric evolution strategies", Machine Learning 75 (2009). Their
// constants: target success rate 2/11, success smoothing 1/12, step damping
// 1 + n/2, covariance rate 2 / (n^2 + 6); the covariance update is skipped
// while the smoothed success rate exceeds 0.44, when steps succeed for
// reasons of scale rather than shape. The steady-state counterpart of CMAES
// for budgets of tens to hundreds of evaluations, where a population method
// completes only a handful of generations.
template <typename Callable, typename RNG, typename scalar_t = double>
class [[maybe_unused]] OnePlusOneCMAES {
 private:
  Callable &f;
  RNG &generator;
  const scalar_t m_step;
  const size_t max_evals;
  const scalar_t x_delta, damping;
  const bool pattern_moves;
  search_distribution<scalar_t> exit_distribution_;
  static constexpr scalar_t kTargetSuccess = 2.0 / 11.0,
                            kSuccessSmoothing = 1.0 / 12.0,
                            kSuccessThreshold = 0.44;
  // pattern move: after a success the displacement is doubled from the
  // accepted point while that keeps improving (the Nelder-Mead expansion
  // coefficient, applied as often as it pays)
  static constexpr scalar_t kExpand = 2.0;
  // budget_damping: the step transient may take this share of the run
  static constexpr scalar_t kTransientShare = 0.5;

 public:
  // Step-size damping d: the step changes by exp((p_s - p_t) / (d (1 - p_t)))
  // per evaluation, so a run of failures shrinks it by a decade in
  // ln(10) (1 - p_t) / p_t * d, about 10 d, evaluations. The reference
  // 1 + n / 2 (passed as 0) buys stationary efficiency at the price of that
  // O(n) transient. `pattern_moves` enables the doubling described above; the
  // expansions do not enter the success-rate statistics, but the covariance
  // learns the whole displacement.
  [[maybe_unused]] OnePlusOneCMAES<Callable, RNG, scalar_t>(
      Callable &f, RNG &generator, const scalar_t m_step = 0.5,
      const size_t max_evals = 1000, const scalar_t x_delta = 1e-12,
      const scalar_t damping = 0.0, const bool pattern_moves = false)
      : f(f),
        generator(generator),
        m_step(m_step),
        max_evals(max_evals),
        x_delta(x_delta),
        damping(damping),
        pattern_moves(pattern_moves) {}
  // The damping for a run of `max_evals` evaluations: the reference 1 + n / 2,
  // capped so that one decade of step adaptation costs at most
  // kTransientShare of the run, and never below 1. A run of tens of
  // evaluations cannot afford the reference transient; one of thousands is
  // not affected.
  static scalar_t budget_damping(const size_t n, const size_t max_evals) {
    const scalar_t decade =
        std::log(10.0) * (1.0 - kTargetSuccess) / kTargetSuccess;
    const scalar_t reference = 1.0 + static_cast<scalar_t>(n) / 2.0;
    const scalar_t affordable =
        kTransientShare * static_cast<scalar_t>(max_evals) / decade;
    return std::clamp(affordable, 1.0, reference);
  }
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }
  // the search distribution when the last minimize / maximize call returned
  [[maybe_unused]] const search_distribution<scalar_t> &exit_distribution()
      const {
    return exit_distribution_;
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const size_t n = x.size();
    const scalar_t nn = static_cast<scalar_t>(n);
    const scalar_t damping =
        this->damping > 0.0 ? this->damping : 1.0 + nn / 2.0;
    const scalar_t c_cov = 2.0 / (nn * nn + 6.0);
    auto project = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n; i++)
          v[i] = std::clamp(v[i], lower[i], upper[i]);
      }
    };
    project(x);
    scalar_t fx = sign * this->f(x);
    size_t evals = 1, iter = 0;
    scalar_t sigma = this->m_step, p_succ = kTargetSuccess;
    // Cholesky factor of the covariance, row-major, starts at the identity
    std::vector<scalar_t> A(n * n, 0.0), z(n), y(n), cand(n), disp(n);
    for (size_t i = 0; i < n; i++) A[i * n + i] = 1.0;
    while (evals < this->max_evals) {
      rnorm<scalar_t>(this->generator, z.begin(), z.end());
      for (size_t i = 0; i < n; i++) {
        scalar_t yi = 0.0;
        for (size_t j = 0; j < n; j++) yi += A[i * n + j] * z[j];
        y[i] = yi;
        cand[i] = x[i] + sigma * yi;
      }
      project(cand);
      const scalar_t fc = sign * this->f(cand);
      evals++;
      const bool success = fc <= fx;
      p_succ = (1.0 - kSuccessSmoothing) * p_succ +
               kSuccessSmoothing * (success ? 1.0 : 0.0);
      sigma *= std::exp((p_succ - kTargetSuccess) /
                        (damping * (1.0 - kTargetSuccess)));
      if (success) {
        // multiple of the sampled step contained in the accepted displacement
        scalar_t k = 1.0;
        if (this->pattern_moves) {
          for (size_t i = 0; i < n; i++) disp[i] = cand[i] - x[i];
        }
        x = cand;
        fx = fc;
        while (this->pattern_moves && evals < this->max_evals) {
          bool moved = false;
          for (size_t i = 0; i < n; i++) cand[i] = x[i] + disp[i];
          project(cand);
          for (size_t i = 0; i < n; i++) moved = moved || cand[i] != x[i];
          if (!moved) break;
          const scalar_t fe = sign * this->f(cand);
          evals++;
          if (fe > fx) break;
          for (size_t i = 0; i < n; i++) disp[i] = cand[i] - x[i] + disp[i];
          x = cand;
          fx = fe;
          k *= kExpand;
        }
        if (p_succ < kSuccessThreshold) {
          // A <- a A + b (A z) z', the rank-one Cholesky update, with z scaled
          // by the displacement multiple k (A (k z) = k y)
          scalar_t z2 = 0.0;
          for (const scalar_t zi : z) z2 += zi * zi;
          z2 *= k * k;
          const scalar_t a = std::sqrt(1.0 - c_cov);
          const scalar_t b =
              a / z2 * (std::sqrt(1.0 + c_cov / (1.0 - c_cov) * z2) - 1.0);
          for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
              A[i * n + j] = a * A[i * n + j] + b * (k * y[i]) * (k * z[j]);
            }
          }
        }
      }
      iter++;
      // largest coordinate standard deviation of the search distribution
      scalar_t max_var = 0.0;
      for (size_t i = 0; i < n; i++) {
        scalar_t row = 0.0;
        for (size_t j = 0; j < n; j++) row += A[i * n + j] * A[i * n + j];
        max_var = std::max(max_var, row);
      }
      if (sigma * std::sqrt(max_var) < this->x_delta || !std::isfinite(sigma)) {
        break;
      }
    }
    // C = A A'
    std::vector<scalar_t> C(n * n, 0.0);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        scalar_t cij = 0.0;
        for (size_t k = 0; k < n; k++) cij += A[i * n + k] * A[j * n + k];
        C[i * n + j] = cij;
      }
    }
    exit_distribution_ = search_distribution<scalar_t>{x, sigma, C};
    return solver_status<scalar_t>(sign * fx, iter, evals);
  }
};

// Pattern search, in its compass form: from the incumbent, poll +h and -h
// along each coordinate in turn and move to the first improvement; after a
// sweep without one the step halves. The sign that last succeeded on a
// coordinate is polled first. [[CITATION]] Hooke, Jeeves, "Direct search
// solution of numerical and statistical problems", J. ACM 8 (1961); Torczon,
// "On the convergence of pattern search algorithms", SIAM J. Optim. 7 (1997);
// Kolda, Lewis, Torczon, "Optimization by direct search", SIAM Review 45
// (2003), whose compass search this is. A poll changes one coordinate, so its
// progress per evaluation carries none of the perpendicular noise that makes
// an isotropic step lose efficiency with dimension, and it tries both signs, so
// it succeeds on any slope. On the humpday demos at 60 to 480 evaluations
// (race_demos.cpp) this alone outranks every stochastic method in the library,
// with or without the rotated disguise, and it is the ingredient Alloy cannot
// lose. Stops when the step falls below `x_delta` (converged) or the budget is
// spent; exit_step() reports the final step. Supports box constraints: a poll
// clipped back onto the incumbent is skipped unevaluated.
template <typename Callable, typename scalar_t = double>
class [[maybe_unused]] PatternSearch {
 private:
  Callable &f;
  const scalar_t m_step;
  const size_t max_evals;
  const scalar_t x_delta;
  scalar_t exit_step_ = 0.0;
  // the step halves after a sweep without improvement, the classical
  // contraction of pattern search
  static constexpr scalar_t kContract = 0.5;

 public:
  [[maybe_unused]] PatternSearch<Callable, scalar_t>(
      Callable &f, const scalar_t m_step = 1.0, const size_t max_evals = 1000,
      const scalar_t x_delta = 1e-12)
      : f(f), m_step(m_step), max_evals(max_evals), x_delta(x_delta) {}
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }
  // the step length when the search stopped
  [[maybe_unused]] scalar_t exit_step() const { return exit_step_; }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const size_t n = x.size();
    if (this->max_evals == 0) {
      // nothing was evaluated, so there is no valid result
      return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                     0, 0, 0, 0, false);
    }
    if constexpr (constrained) {
      for (size_t d = 0; d < n; d++)
        x[d] = std::clamp(x[d], lower[d], upper[d]);
    }
    scalar_t fx = sign * this->f(x);
    size_t evals = 1, sweeps = 0;
    scalar_t step = this->m_step;
    // the probe equals the incumbent except in the coordinate being polled
    std::vector<scalar_t> probe = x, first_sign(n, 1.0);
    while (evals < this->max_evals && step >= this->x_delta) {
      bool improved = false;
      for (size_t d = 0; d < n && evals < this->max_evals; d++) {
        for (int k = 0; k < 2 && evals < this->max_evals; k++) {
          const scalar_t dir = k == 0 ? first_sign[d] : -first_sign[d];
          probe[d] = x[d] + dir * step;
          if constexpr (constrained) {
            probe[d] = std::clamp(probe[d], lower[d], upper[d]);
          }
          if (probe[d] == x[d]) continue;
          const scalar_t fp = sign * this->f(probe);
          evals++;
          if (fp < fx) {
            x[d] = probe[d];
            fx = fp;
            first_sign[d] = dir;
            improved = true;
            break;
          }
          probe[d] = x[d];
        }
      }
      sweeps++;
      if (!improved) step *= kContract;
    }
    exit_step_ = step;
    return solver_status<scalar_t>(sign * fx, sweeps, evals, 0, 0,
                                   step < this->x_delta);
  }
};

// Adaptive Coordinate Descent (Loshchilov, Schoenauer, Sebag, 2011): coordinate
// descent whose coordinate system is continuously rotated by an Adaptive
// Encoding covariance update - the same covariance / eigendecomposition
// machinery as CMA-ES (here reusing tinyqr::QRSolver). Each iteration performs
// one Gauss-Seidel sweep of derivative-free (golden-section) line searches
// along the adapted principal axes, then ranks the resulting points and feeds
// them to the encoding update. Box constraints (optional) are handled by
// projecting every evaluated point into the feasible box.
template <typename Callable, typename scalar_t = double>
class [[maybe_unused]] AdaptiveCoordinateDescent {
 private:
  Callable &f;
  scalar_t m_step;
  const size_t max_iter;
  const scalar_t x_delta, f_delta;
  // Each coordinate line search stops once its bracket is below `line_tol`
  // times the sweep's initial step (or 1e-12 relative, whichever is larger).
  // Zero keeps the full 1e-12 resolution, about 50 evaluations per axis and
  // sweep, which is what lets the method reach 1e-10 on narrow valleys; 1e-2
  // cuts the cost per sweep by three to four times but on Rosenbrock the
  // sweeps then stall near 1e-7, since the step-size feedback tracks the
  // resolution of the line optima.
  const scalar_t line_tol;

 public:
  [[maybe_unused]] AdaptiveCoordinateDescent<Callable, scalar_t>(
      Callable &f, const scalar_t m_step = 1.0, const size_t max_iter = 1000,
      const scalar_t x_delta = 1e-12, const scalar_t f_delta = 1e-12,
      const scalar_t line_tol = 0.0)
      : f(f),
        m_step(m_step),
        max_iter(max_iter),
        x_delta(x_delta),
        f_delta(f_delta),
        line_tol(line_tol) {}
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    const size_t n = x.size();
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const scalar_t nn = static_cast<scalar_t>(n);
    // ---- selection / adaptation constants (as in CMA-ES) ----
    const size_t mu = std::max<size_t>(1, n / 2);
    std::vector<scalar_t> w(mu);
    for (size_t i = 0; i < mu; i++) {
      w[i] = std::log(static_cast<scalar_t>(mu) + 0.5) -
             std::log(static_cast<scalar_t>(i) + 1.0);
    }
    scalar_t ws = 0;
    for (const auto v : w) ws += v;
    for (auto &v : w) v /= ws;
    scalar_t mu_eff = 0;
    for (const auto v : w) mu_eff += v * v;
    mu_eff = 1.0 / mu_eff;
    const scalar_t cc = (4.0 + mu_eff / nn) / (nn + 4.0 + 2.0 * mu_eff / nn);
    const scalar_t c1 = 2.0 / (std::pow(nn + 1.3, 2.0) + mu_eff);
    const scalar_t cmu =
        std::min(1.0 - c1, 2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) /
                               (std::pow(nn + 2.0, 2.0) + mu_eff));
    // ---- state: encoding (B columns = axes, d = axis scales), covariance C,
    // evolution path pc, distribution mean ----
    std::vector<scalar_t> mean = x;
    auto project = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n; i++)
          v[i] =
              v[i] < lower[i] ? lower[i] : (v[i] > upper[i] ? upper[i] : v[i]);
      }
    };
    project(mean);
    std::vector<scalar_t> C(n * n, 0.0), B(n * n, 0.0), d(n, 1.0), pc(n, 0.0);
    for (size_t i = 0; i < n; i++) {
      C[i * n + i] = 1.0;
      B[i * n + i] = 1.0;
    }
    scalar_t sigma = m_step;
    const scalar_t sigma_cap = m_step * static_cast<scalar_t>(1e2);
    tinyqr::QRSolver<scalar_t> qr(n);
    std::vector<scalar_t> best_x = mean;
    scalar_t best_f = sign * this->f(mean);
    size_t f_evals = 1;
    std::vector<scalar_t> buf(n), Btv(n), tmp(n), dm(n), diff(n), work(n);

    auto eval_at = [&](std::vector<scalar_t> &p) -> scalar_t {
      project(p);
      f_evals++;
      const scalar_t fv = sign * this->f(p);
      if (fv < best_f) {
        best_f = fv;
        best_x = p;
      }
      return fv;
    };
    // C^{-1/2} v = B diag(1/d) B^T v
    auto Cinvhalf = [&](const std::vector<scalar_t> &v,
                        std::vector<scalar_t> &out) {
      for (size_t j = 0; j < n; j++) {
        scalar_t s = 0;
        for (size_t i = 0; i < n; i++) s += B[i * n + j] * v[i];
        Btv[j] = s / d[j];
      }
      for (size_t i = 0; i < n; i++) {
        scalar_t s = 0;
        for (size_t j = 0; j < n; j++) s += B[i * n + j] * Btv[j];
        out[i] = s;
      }
    };
    // golden-section line search along principal axis k from `base`; h0 is the
    // initial step. Returns t*, writes the located point to `out` and its
    // objective value (already evaluated during the search) to `f_out`.
    auto line_search = [&](const std::vector<scalar_t> &base, const size_t k,
                           const scalar_t h0, std::vector<scalar_t> &out,
                           scalar_t &f_out) -> scalar_t {
      auto phi = [&](const scalar_t t) -> scalar_t {
        for (size_t i = 0; i < n; i++) buf[i] = base[i] + t * B[i * n + k];
        return eval_at(buf);
      };
      const scalar_t gr = 0.6180339887498949;
      const scalar_t h = h0 > 1e-300 ? h0 : 1e-4;
      const scalar_t f0 = phi(0.0), fp = phi(h), fm = phi(-h);
      scalar_t a, b, c, fb, fc;
      if (f0 <= fp && f0 <= fm) {  // minimum already bracketed in [-h, h]
        a = -h;
        b = 0;
        c = h;
      } else {
        const scalar_t step = (fp <= fm) ? h : -h;
        // cap the bracketing expansion at a finite multiple of the initial
        // step, so a single line search cannot run off to huge magnitudes on
        // deceptive / flat-rayed objectives (e.g. unconstrained Beale)
        const scalar_t max_step = h * static_cast<scalar_t>(1e2);
        a = 0;
        b = step;
        fb = (fp <= fm) ? fp : fm;
        c = b + step / gr;
        fc = phi(c);
        while (fc < fb && std::abs(c) < max_step) {
          a = b;
          b = c;
          fb = fc;
          c = b + (b - a) / gr;
          fc = phi(c);
        }
      }
      scalar_t lo = std::min(a, c), hi = std::max(a, c);
      scalar_t x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo);
      scalar_t f1 = phi(x1), f2 = phi(x2);
      const scalar_t tol =
          std::max(line_tol * h, 1e-12 * (std::abs(lo) + std::abs(hi) + 1e-12));
      for (size_t it = 0; it < 50 && (hi - lo) > tol; it++) {
        if (f1 < f2) {
          hi = x2;
          x2 = x1;
          f2 = f1;
          x1 = hi - gr * (hi - lo);
          f1 = phi(x1);
        } else {
          lo = x1;
          x1 = x2;
          f1 = f2;
          x2 = lo + gr * (hi - lo);
          f2 = phi(x2);
        }
      }
      // the line optimum is the best of the two bracket points and the start
      // itself; a bracket resolved only to the tolerance above may straddle a
      // narrow valley with both of its points above f(base), and moving there
      // would climb
      scalar_t t = (f1 < f2) ? x1 : x2;
      f_out = (f1 < f2) ? f1 : f2;
      if (f0 <= f_out) {
        t = 0.0;
        f_out = f0;
      }
      for (size_t i = 0; i < n; i++) out[i] = base[i] + t * B[i * n + k];
      project(out);  // phi evaluated this same projected point
      return t;
    };

    std::vector<std::vector<scalar_t>> offspring(n, std::vector<scalar_t>(n));
    std::vector<scalar_t> off_f(n), step_len(n);
    size_t iter = 0;
    // Best-value stagnation. A sweep that fails to improve the best value by
    // f_delta (relative) leaves the next sweep only a shifted mean and a
    // shrinking initial bracket, and the golden-section searches of the failed
    // sweep already covered those brackets to 1e-12; so a few such sweeps in a
    // row mean convergence, or a stall no further sweep can break. A patience
    // of 30 + 4n sweeps spent 40 to 70 full sweeps (thousands of evaluations)
    // confirming a result reached in the first one.
    scalar_t last_improved_f = best_f;
    size_t stall = 0;
    constexpr size_t kStallSweeps = 3;
    while (true) {
      // ---- one coordinate-descent sweep over the n adapted principal axes
      // ----
      work = mean;
      for (size_t k = 0; k < n; k++) {
        const scalar_t h0 = sigma * d[k];
        const scalar_t t = line_search(work, k, h0, offspring[k], off_f[k]);
        step_len[k] = std::abs(t);
        work = offspring[k];  // Gauss-Seidel: descend immediately
      }
      // ---- adaptive encoding update from the mu best sweep points ----
      const std::vector<size_t> idx = index_partial_sort(off_f, mu);
      const std::vector<scalar_t> m_old = mean;
      for (size_t i = 0; i < n; i++) {
        scalar_t s = 0;
        for (size_t r = 0; r < mu; r++) s += w[r] * offspring[idx[r]][i];
        mean[i] = s;
      }
      // evolution path from the (per-vector normalized) mean shift
      for (size_t i = 0; i < n; i++) dm[i] = mean[i] - m_old[i];
      Cinvhalf(dm, tmp);
      scalar_t zn = 0;
      for (const auto v : tmp) zn += v * v;
      zn = std::sqrt(zn);
      const scalar_t a0 = zn > 1e-300 ? std::sqrt(nn) / zn : 0.0;
      const scalar_t pcoef = std::sqrt(cc * (2.0 - cc) * mu_eff);
      for (size_t i = 0; i < n; i++)
        pc[i] = (1.0 - cc) * pc[i] + pcoef * a0 * dm[i];
      // covariance: C = (1-c1-cmu) C + c1 pc pc^T + cmu sum w_i y_i y_i^T,
      // y_i the per-vector normalized displacement of selected offspring
      for (auto &cij : C) cij *= (1.0 - c1 - cmu);
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) C[i * n + j] += c1 * pc[i] * pc[j];
      for (size_t r = 0; r < mu; r++) {
        for (size_t i = 0; i < n; i++)
          diff[i] = offspring[idx[r]][i] - m_old[i];
        Cinvhalf(diff, tmp);
        scalar_t zz = 0;
        for (const auto v : tmp) zz += v * v;
        zz = std::sqrt(zz);
        const scalar_t ai = zz > 1e-300 ? std::sqrt(nn) / zz : 0.0;
        const scalar_t cw = cmu * w[r] * ai * ai;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++) C[i * n + j] += cw * diff[i] * diff[j];
      }
      // ---- eigendecompose C -> B (columns), d (sqrt eigenvalues) ----
      for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++) {
          const scalar_t av = 0.5 * (C[i * n + j] + C[j * n + i]);
          C[i * n + j] = C[j * n + i] = av;
        }
      qr.solve(C, 100, 1e-13);
      const auto &ev = qr.eigenvalues();
      const auto &QQ = qr.eigenvectors();  // eigenvectors are ROWS of QQ
      scalar_t dgeo = 0;
      for (size_t j = 0; j < n; j++) {
        d[j] = std::sqrt(std::max(ev[j], static_cast<scalar_t>(1e-30)));
        dgeo += std::log(d[j]);
      }
      dgeo = std::exp(dgeo / nn);  // geometric-mean axis scale
      for (size_t j = 0; j < n; j++) {
        d[j] /= dgeo;
        for (size_t i = 0; i < n; i++) B[i * n + j] = QQ[j * n + i];
      }
      for (auto &cij : C) cij /= (dgeo * dgeo);  // keep C unit-scaled
      // ---- step size: track the typical successful line-search step ----
      scalar_t smean = 0;
      for (size_t k = 0; k < n; k++) smean += step_len[k];
      smean /= nn;
      // bound sigma so the step-size cannot inflate without limit (the
      // smean feedback would otherwise let a descending ray run away)
      sigma =
          0.5 * sigma + 0.5 * std::max(smean, static_cast<scalar_t>(1e-300));
      sigma = std::min(sigma, sigma_cap);
      // ---- termination ----
      scalar_t dmax = d[0], dmin = d[0];
      for (size_t j = 1; j < n; j++) {
        dmax = std::max(dmax, d[j]);
        dmin = std::min(dmin, d[j]);
      }
      iter++;
      // stagnation: count iterations without a meaningful best-value gain
      if (last_improved_f - best_f > this->f_delta * (std::abs(best_f) + 1.0)) {
        last_improved_f = best_f;
        stall = 0;
      } else {
        stall++;
      }
      if (iter >= this->max_iter || sigma * dmax < this->x_delta ||
          !std::isfinite(sigma) || (dmax / dmin) > 1e14 ||
          stall >= kStallSweeps) {
        x = best_x;
        m_step = sigma;
        return solver_status<scalar_t>(sign * best_f, iter, f_evals);
      }
    }
  }
};

// Projected Armijo backtracking from x along dir, starting at `step`. Trial
// points are projected by `project` (the identity when unconstrained) and
// accepted once f(x_new) <= f(x) + c1 grad'(x_new - x). The step halves until
// acceptance or until the trial displacement falls below sqrt(eps) relative to
// x, the finite-difference resolution below which no decrease consistent with
// the gradient can be verified; the first trial is always taken. Returns true
// on acceptance, with x_new and f_new set; `may_evaluate` lets a caller stop
// on an evaluation budget.
template <typename Callable, typename Project, typename MayEvaluate,
          typename scalar_t>
static inline bool armijo_backtrack(
    Callable &f, const std::vector<scalar_t> &x, const scalar_t fx,
    const std::vector<scalar_t> &grad, const std::vector<scalar_t> &dir,
    scalar_t step, Project &project, MayEvaluate &may_evaluate,
    std::vector<scalar_t> &x_new, scalar_t &f_new) {
  constexpr scalar_t c1 = 1e-4;
  const scalar_t min_rel_step =
      std::sqrt(std::numeric_limits<scalar_t>::epsilon());
  const size_t n = x.size();
  scalar_t dir_inf = 0, x_inf = 0;
  for (size_t j = 0; j < n; j++) {
    dir_inf = std::max(dir_inf, std::abs(dir[j]));
    x_inf = std::max(x_inf, std::abs(x[j]));
  }
  do {
    for (size_t j = 0; j < n; j++) x_new[j] = x[j] + step * dir[j];
    project(x_new);
    const scalar_t f_trial = f(x_new);
    scalar_t pred = 0;  // predicted decrease along the (projected) step
    for (size_t j = 0; j < n; j++) pred += grad[j] * (x_new[j] - x[j]);
    if (f_trial <= fx + c1 * pred) {
      f_new = f_trial;
      return true;
    }
    step *= 0.5;
  } while (step * dir_inf > min_rel_step * (1.0 + x_inf) && may_evaluate());
  return false;
}

// L-BFGS-B: limited-memory BFGS with box constraints. This is the practical
// "projected" variant - the L-BFGS two-loop recursion gives the quasi-Newton
// search direction and a projected Armijo line search keeps iterates feasible
// (rather than the full Byrd et al. generalized-Cauchy-point / subspace
// minimization). Unconstrained (no bounds) it is plain L-BFGS. Only the most
// recent `history_size` curvature pairs are stored, so memory is O(history*n).
// The default gradient is a two-point central difference, 2n evaluations per
// gradient. An optional `initial_inverse_hessian` (row-major n x n) replaces
// the identity in the two-loop recursion; it is rescaled by s'y / y'H0y from
// the latest curvature pair, the usual gamma scaling with H0 in place of I.
// `max_evals` caps objective evaluations (a gradient or line search in
// progress may overrun it by at most one gradient plus one trial step).
template <typename Callable, typename scalar_t = double,
          typename Grad = fin_diff<Callable, scalar_t, 0>>
class [[maybe_unused]] LBFGSB {
 private:
  Callable &f;
  Grad g;
  const size_t max_iter, history_size;
  const scalar_t grad_eps;
  const std::vector<scalar_t> h0;  // empty means the identity
  const size_t max_evals;

 public:
  [[maybe_unused]] explicit LBFGSB<Callable, scalar_t, Grad>(
      Callable &f, Grad g = Grad(), const size_t max_iter = 200,
      const scalar_t grad_eps = 1e-7, const size_t history_size = 10,
      std::vector<scalar_t> initial_inverse_hessian = {},
      const size_t max_evals = std::numeric_limits<size_t>::max())
      : f(f),
        g(g),
        max_iter(max_iter),
        history_size(history_size),
        grad_eps(grad_eps),
        h0(std::move(initial_inverse_hessian)),
        max_evals(max_evals) {}
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<true, false>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;
    return this->solve<false, false>(x, lower, upper);
  }
  // box-constrained interfaces; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    const size_t n = x.size();
    const int i_dim = static_cast<int>(n);
    if (!h0.empty() && h0.size() != n * n) {
      throw std::invalid_argument(
          "LBFGSB: initial_inverse_hessian must be n x n for an n-dimensional "
          "x");
    }
    // maximize f by minimizing -f (the library-wide convention)
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    size_t function_calls_used = 0, grad_evals_used = 0;
    auto f_lam = [&](decltype(x) &coef) {
      function_calls_used++;
      return sign * this->f(coef);
    };
    auto g_lam = [&](decltype(x) &coef, std::vector<scalar_t> &grad) {
      grad_evals_used++;
      if constexpr (is_fin_diff<Grad>::value) {
        // f_lam already carries the sign, so this is grad of sign*f, with the
        // same stencil the user selected
        fin_diff<decltype(f_lam), scalar_t, is_fin_diff<Grad>::stencil>()(
            f_lam, coef, grad);
        return;
      }
      if constexpr (is_fin_diff_fwd<Grad>::value) {
        fin_diff_fwd<decltype(f_lam), scalar_t>()(f_lam, coef, grad);
        return;
      }
      this->g(this->f, coef, grad);
      if constexpr (!minimize) {
        for (auto &gi : grad) gi = -gi;
      }
    };
    auto project = [&](std::vector<scalar_t> &v) {
      if constexpr (constrained) {
        for (size_t i = 0; i < n; i++)
          v[i] =
              v[i] < lower[i] ? lower[i] : (v[i] > upper[i] ? upper[i] : v[i]);
      }
    };
    project(x);
    std::vector<scalar_t> grad(n, 0.0), grad_new(n, 0.0), dir(n, 0.0),
        q(n, 0.0), x_new(n, 0.0), h0q(n, 0.0);
    // objective at the current iterate; carried forward from the accepted line
    // search step so the point is never re-evaluated
    scalar_t fx = f_lam(x);
    g_lam(x, grad);
    // limited-memory curvature pairs (s = step, y = gradient change)
    std::vector<std::vector<scalar_t>> S, Y;
    std::vector<scalar_t> rho;
    size_t iter = 0;
    while (true) {
      // projected-gradient stationarity measure: || x - P(x - grad) ||_inf
      scalar_t pg_norm = 0;
      for (size_t i = 0; i < n; i++) {
        scalar_t xi = x[i] - grad[i];
        if constexpr (constrained)
          xi = xi < lower[i] ? lower[i] : (xi > upper[i] ? upper[i] : xi);
        pg_norm = std::max(pg_norm, std::abs(x[i] - xi));
      }
      // success() reports whether the projected gradient fell below grad_eps,
      // as opposed to an exit on the iteration or evaluation cap
      if (iter >= this->max_iter || pg_norm < this->grad_eps ||
          function_calls_used >= this->max_evals) {
        return solver_status<scalar_t>(sign * fx, iter, function_calls_used,
                                       grad_evals_used, 0,
                                       pg_norm < this->grad_eps);
      }
      // ---- two-loop recursion: dir = -H * grad ----
      q = grad;
      const size_t k = S.size();
      std::vector<scalar_t> alpha_h(k);
      for (size_t ii = k; ii-- > 0;) {
        alpha_h[ii] = rho[ii] * math::dot(S[ii].data(), q.data(), i_dim);
        for (size_t j = 0; j < n; j++) q[j] -= alpha_h[ii] * Y[ii][j];
      }
      // initial inverse Hessian gamma * H0, with H0 the identity or the seed;
      // gamma = s'y / y'H0y matches it to the latest curvature pair
      scalar_t gamma = 1.0;
      if (!h0.empty()) {
        for (size_t i = 0; i < n; i++) {
          h0q[i] = math::dot(h0.data() + i * n, q.data(), i_dim);
        }
        q = h0q;
        if (k > 0) {
          const scalar_t ys =
              math::dot(Y[k - 1].data(), S[k - 1].data(), i_dim);
          scalar_t yhy = 0.0;
          for (size_t i = 0; i < n; i++) {
            yhy += Y[k - 1][i] *
                   math::dot(h0.data() + i * n, Y[k - 1].data(), i_dim);
          }
          if (yhy > 0) gamma = ys / yhy;
        }
      } else if (k > 0) {
        const scalar_t ys = math::dot(Y[k - 1].data(), S[k - 1].data(), i_dim);
        const scalar_t yy = math::dot(Y[k - 1].data(), Y[k - 1].data(), i_dim);
        if (yy > 0) gamma = ys / yy;
      }
      for (size_t j = 0; j < n; j++) q[j] *= gamma;
      for (size_t ii = 0; ii < k; ii++) {
        const scalar_t beta =
            rho[ii] * math::dot(Y[ii].data(), q.data(), i_dim);
        for (size_t j = 0; j < n; j++) q[j] += S[ii][j] * (alpha_h[ii] - beta);
      }
      for (size_t j = 0; j < n; j++) dir[j] = -q[j];
      // fall back to steepest descent if the model direction is not downhill
      if (math::dot(dir.data(), grad.data(), i_dim) >= 0) {
        for (size_t j = 0; j < n; j++) dir[j] = -grad[j];
      }
      // ---- projected backtracking (Armijo) line search ----
      // first step (no curvature info yet) scaled by the gradient magnitude
      scalar_t g_inf = 0;
      for (size_t j = 0; j < n; j++) g_inf = std::max(g_inf, std::abs(grad[j]));
      const scalar_t step =
          (k == 0) ? 1.0 / std::max(static_cast<scalar_t>(1), g_inf) : 1.0;
      auto may_evaluate = [&]() {
        return function_calls_used < this->max_evals;
      };
      scalar_t f_new = fx;
      const bool accepted = armijo_backtrack(
          f_lam, x, fx, grad, dir, step, project, may_evaluate, x_new, f_new);
      if (!accepted) {
        // no decrease consistent with the gradient could be found: not a
        // convergence, the gradient is unreliable here (noise, a kink)
        return solver_status<scalar_t>(sign * fx, iter, function_calls_used,
                                       grad_evals_used, 0, false);
      }
      g_lam(x_new, grad_new);
      // ---- curvature pair update (skip if it would break
      // positive-definiteness)
      std::vector<scalar_t> s(n), y(n);
      for (size_t j = 0; j < n; j++) {
        s[j] = x_new[j] - x[j];
        y[j] = grad_new[j] - grad[j];
      }
      const scalar_t sy = math::dot(s.data(), y.data(), i_dim);
      const scalar_t y_norm = math::norm(y.data(), i_dim);
      const scalar_t s_norm = math::norm(s.data(), i_dim);
      if (sy > std::numeric_limits<scalar_t>::epsilon() * y_norm * s_norm) {
        S.push_back(s);
        Y.push_back(y);
        rho.push_back(1.0 / sy);
        if (S.size() > this->history_size) {
          S.erase(S.begin());
          Y.erase(Y.begin());
          rho.erase(rho.begin());
        }
      }
      x = x_new;
      fx = f_new;
      grad = grad_new;
      iter++;
    }
  }
};

// Alloy, a blended derivative-free method for small evaluation budgets. The
// program was generated by a language model asked to blend Nelder-Mead,
// Differential Evolution, CMA-ES, pattern search and simulated annealing in
// equal parts, then selected on disguised objectives and validated on
// twenty-nine held-out problems at budgets of 60 to 480 evaluations, where it
// had the best mean rank at every budget.
// [[CITATION]] Cotton, P., "The Inspiration Simplex", working paper, and the
// reference implementation humpday/optimizers/alloy.py
// (https://github.com/microprediction/humpday). Control flow and constants are
// kept verbatim so that the validated behaviour carries over.
//
// Mechanism. Draw a uniform population on the box and keep the n+1 best points
// as a simplex. Each step applies one of four moves, chosen uniformly at
// random: a Nelder-Mead reflect / expand / contract / shrink; a Differential
// Evolution mutation with binomial crossover; a Gaussian perturbation of the
// best vertex with a success-adapted diagonal covariance; or a Hooke-Jeeves
// coordinate probe with a shrinking step. Every candidate replaces the worst
// vertex only if it passes a Metropolis test under a geometrically cooled
// temperature. After a run of non-improving steps the temperature is reheated
// and the simplex is rebuilt around the best point seen so far.
//
// The search lives on the unit cube [0, 1]^n, mapped affinely onto the user's
// box at evaluation time, so every step constant below is a fraction of the box
// width. The bounded overloads clip each candidate to the cube, as the
// reference does. The unbounded overloads centre a box of half-width
// clamp(||x0||_inf, 1, 10) on the starting point (the scale rule of the
// Nelder-Mead simplex initialisation above), take the step scale from it, and
// do not clip, so the search may leave that box. There is no convergence test:
// the method spends the whole evaluation budget and reports the best point
// seen.
template <typename Callable, typename RNG, typename scalar_t = double>
class [[maybe_unused]] Alloy {
 private:
  Callable &f;
  RNG &generator;
  const size_t max_evals;

  // Reference constants (humpday alloy.py). The method was validated as a
  // whole; none of these was tuned on its own, so none should move alone.
  // move selection: four generators with equal probability
  static constexpr scalar_t kMoveShare = 0.25;
  // Nelder-Mead coefficients
  static constexpr scalar_t kReflect = 1.0, kExpand = 2.0, kContract = 0.5,
                            kShrink = 0.5;
  // Differential Evolution: weight F, crossover rate CR, and the probability of
  // the rand/1 mutant over the current-to-best/1 mutant
  static constexpr scalar_t kDiffWeight = 0.6, kCrossover = 0.9,
                            kRandMutantShare = 0.5;
  // Gaussian sampler: global scale sigma with its growth on success, decay on
  // failure and bounds; diagonal covariance start, memory, innovation and floor
  static constexpr scalar_t kSigma0 = 0.2, kSigmaGrow = 1.05,
                            kSigmaDecay = 0.97, kSigmaMax = 0.5,
                            kSigmaMin = 1e-3;
  static constexpr scalar_t kCov0 = 0.04, kCovMemory = 0.8,
                            kCovInnovation = 0.2, kCovFloor = 1e-8;
  // Hooke-Jeeves probe: step, halving on failure, reset once below the minimum
  static constexpr scalar_t kStep0 = 0.25, kStepShrink = 0.5, kStepMin = 1e-6;
  // annealing: T0 = max(kTemp0Min, initial simplex spread + kTemp0Pad), cooled
  // by kCooling per step down to kTempMin; no uphill moves below kTempAcceptEps
  static constexpr scalar_t kTemp0Min = 1e-6, kTemp0Pad = 1e-3,
                            kCooling = 0.995, kTempMin = 1e-9,
                            kTempAcceptEps = 1e-12;
  // restart after max(kStallMin, kStallPerDim * n) non-improving steps: reheat
  // to max(kReheatFrac * T0, kReheatMult * T) and jitter the best point by
  // uniform(-kRestartJitter, kRestartJitter) per coordinate
  static constexpr scalar_t kReheatFrac = 0.5, kReheatMult = 5.0,
                            kRestartJitter = 0.3;
  static constexpr size_t kStallMin = 12, kStallPerDim = 4;
  // population size max(n + 1, min(kPopBase + kPopPerDim * n,
  // max(kPopMin, budget / kPopBudgetShare))); a short simplex is completed by
  // stepping the best member kSimplexFill along one coordinate at a time
  static constexpr size_t kPopBase = 8, kPopPerDim = 2, kPopMin = 5,
                          kPopBudgetShare = 6;
  static constexpr scalar_t kSimplexFill = 0.1;
  // unbounded overloads: half-width of the box centred on the starting point
  static constexpr scalar_t kHalfWidthMin = 1.0, kHalfWidthMax = 10.0;

 public:
  [[maybe_unused]] Alloy<Callable, RNG, scalar_t>(Callable &f, RNG &generator,
                                                  const size_t max_evals = 5000)
      : f(f), generator(generator), max_evals(max_evals) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;  // empty bounds => unconstrained
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;  // empty bounds => unconstrained
    return this->solve<false, false>(x, lower, upper);
  }
  // minimize with box constraints; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  // maximize with box constraints; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    // we always minimize internally; for maximization we minimize -f
    constexpr scalar_t f_multiplier = minimize ? 1.0 : -1.0;
    constexpr scalar_t nan = std::numeric_limits<scalar_t>::quiet_NaN();
    constexpr scalar_t inf = std::numeric_limits<scalar_t>::infinity();
    const size_t n = x.size();
    size_t f_evals = 0;
    if (n == 0) {
      // nothing to search over; the reference spends one evaluation on the
      // empty point if the budget allows it
      scalar_t f_value = nan;
      if (this->max_evals > 0) {
        f_value = this->f(x);
        f_evals = 1;
      }
      return solver_status<scalar_t>(f_value, 0, f_evals);
    }
    // unit cube -> box: x = origin + u * width
    std::vector<scalar_t> origin(n), width(n);
    if constexpr (constrained) {
      for (size_t i = 0; i < n; i++) {
        origin[i] = lower[i];
        width[i] = upper[i] - lower[i];
      }
    } else {
      const scalar_t half_width =
          std::clamp<scalar_t>(max_abs_vec(x), kHalfWidthMin, kHalfWidthMax);
      for (size_t i = 0; i < n; i++) {
        origin[i] = x[i] - half_width;
        width[i] = 2.0 * half_width;
      }
    }
    // simplex storage: n + 1 vertices, vertex `id` occupies row id * n
    const size_t n_vertices = n + 1;
    std::vector<scalar_t> vertices(n_vertices * n), vertex_f(n_vertices);
    std::vector<size_t> order(n_vertices), shuffled(n_vertices);
    std::vector<scalar_t> eval_buf(n), best_u(n, 0.0), centroid(n),
        reflected(n), candidate(n), delta(n), current(n), probe(n);
    scalar_t best_f = inf;

    // clip to the cube (bounded overloads), map onto the box, evaluate, and
    // keep the best point seen. A NaN objective value counts as +infinity so
    // the vertex ordering below stays well defined.
    auto evaluate = [&](std::vector<scalar_t> &u) -> scalar_t {
      for (size_t i = 0; i < n; i++) {
        if constexpr (constrained) u[i] = std::clamp<scalar_t>(u[i], 0.0, 1.0);
        eval_buf[i] = origin[i] + u[i] * width[i];
      }
      scalar_t value = f_multiplier * this->f(eval_buf);
      f_evals++;
      if (std::isnan(value)) value = inf;
      if (value < best_f) {
        best_f = value;
        best_u = u;
      }
      return value;
    };
    auto budget_left = [&]() { return f_evals < this->max_evals; };
    auto store_vertex = [&](const size_t id, const std::vector<scalar_t> &u,
                            const scalar_t value) {
      std::copy(u.begin(), u.end(), vertices.begin() + id * n);
      vertex_f[id] = value;
    };

    // ---- initial population (Differential Evolution graft) ----
    const size_t pop_size = std::max(
        n_vertices,
        std::min(kPopBase + kPopPerDim * n,
                 std::max(kPopMin, this->max_evals / kPopBudgetShare)));
    std::vector<scalar_t> population(pop_size * n), population_f(pop_size);
    size_t pop_count = 0;
    for (; pop_count < pop_size && budget_left(); pop_count++) {
      for (size_t d = 0; d < n; d++) candidate[d] = this->generator();
      population_f[pop_count] = evaluate(candidate);
      std::copy(candidate.begin(), candidate.end(),
                population.begin() + pop_count * n);
    }
    if (pop_count == 0) {
      // no budget at all: nothing was evaluated, so there is no valid result
      return solver_status<scalar_t>(nan, 0, 0, 0, 0, false);
    }
    // the simplex is the n + 1 best members. Ties fall back to population
    // order, which is what the reference's stable sort gives; sorting by
    // (value, id) reproduces it without stable_sort's per-call allocation
    std::vector<size_t> pop_order(pop_count);
    std::iota(pop_order.begin(), pop_order.end(), 0);
    std::sort(pop_order.begin(), pop_order.end(),
              [&](const size_t a, const size_t b) {
                return population_f[a] < population_f[b] ||
                       (population_f[a] == population_f[b] && a < b);
              });
    size_t n_filled = std::min(n_vertices, pop_count);
    for (size_t k = 0; k < n_filled; k++) {
      const size_t row = pop_order[k] * n;
      std::copy(population.begin() + row, population.begin() + row + n,
                vertices.begin() + k * n);
      vertex_f[k] = population_f[pop_order[k]];
    }
    // a population smaller than n + 1 (tiny budget) is completed by stepping
    // the best member along one coordinate at a time
    for (; n_filled < n_vertices && budget_left(); n_filled++) {
      std::copy(vertices.begin(), vertices.begin() + n, candidate.begin());
      const size_t j = (n_filled - 1) % n;
      candidate[j] += kSimplexFill;
      const scalar_t value = evaluate(candidate);
      store_vertex(n_filled, candidate, value);
    }

    // ---- adaptation state ----
    scalar_t step = kStep0, sigma = kSigma0;
    std::vector<scalar_t> cov_diag(n, kCov0);
    // annealing temperature from the spread of the initial simplex
    scalar_t f_lo = vertex_f[0], f_hi = vertex_f[0];
    for (size_t k = 1; k < n_filled; k++) {
      f_lo = std::min(f_lo, vertex_f[k]);
      f_hi = std::max(f_hi, vertex_f[k]);
    }
    const scalar_t temp0 = std::max(kTemp0Min, (f_hi - f_lo) + kTemp0Pad);
    scalar_t temp = temp0;
    // Metropolis gate against the worst vertex: an improvement is always
    // taken, a worse point with probability exp(-(f_new - f_old) / T)
    auto accept = [&](const scalar_t f_new, const scalar_t f_old) -> bool {
      if (f_new <= f_old) return true;
      if (temp <= kTempAcceptEps) return false;
      return this->generator() < std::exp(-(f_new - f_old) / temp);
    };
    const size_t stall_limit = std::max(kStallMin, kStallPerDim * n);
    size_t stagnation = 0, iter = 0;

    // Vertex ids by value, best first, worst last; ties by id, which is the
    // order the reference's stable sort gives. `order` stays a permutation
    // across iterations and is nearly sorted at each pass (usually one vertex
    // moved), so an ordinary sort reproduces that order exactly, without
    // stable_sort's per-call heap allocation.
    std::iota(order.begin(), order.end(), 0);
    auto by_value_then_id = [&](const size_t a, const size_t b) {
      return vertex_f[a] < vertex_f[b] || (vertex_f[a] == vertex_f[b] && a < b);
    };
    // Whenever budget remains at the top of the loop the simplex is complete:
    // the fill loops above and the restart below only stop short of n + 1
    // vertices when the budget is exhausted.
    while (budget_left()) {
      std::sort(order.begin(), order.end(), by_value_then_id);
      const size_t best_id = order[0], worst_id = order[n];
      const size_t best_row = best_id * n, worst_row = worst_id * n;
      const scalar_t worst_f = vertex_f[worst_id];

      const scalar_t move = this->generator();
      bool improved = false;
      if (move < kMoveShare) {
        // ---- Nelder-Mead: reflect, then expand, contract or shrink ----
        // centroid of every vertex but the worst (only this move needs it)
        std::fill(centroid.begin(), centroid.end(), 0.0);
        for (size_t id = 0; id < n_vertices; id++) {
          if (id == worst_id) continue;
          for (size_t d = 0; d < n; d++) centroid[d] += vertices[id * n + d];
        }
        for (auto &c : centroid) c /= static_cast<scalar_t>(n);
        for (size_t d = 0; d < n; d++) {
          reflected[d] =
              centroid[d] + kReflect * (centroid[d] - vertices[worst_row + d]);
        }
        const scalar_t f_reflected = evaluate(reflected);
        const std::vector<scalar_t> *chosen = nullptr;
        scalar_t chosen_f = inf;
        if (f_reflected < vertex_f[best_id]) {
          chosen = &reflected;
          chosen_f = f_reflected;
          if (budget_left()) {
            for (size_t d = 0; d < n; d++) {
              candidate[d] = centroid[d] +
                             kExpand * (centroid[d] - vertices[worst_row + d]);
            }
            const scalar_t f_expanded = evaluate(candidate);
            if (f_expanded < f_reflected) {
              chosen = &candidate;
              chosen_f = f_expanded;
            }
          }
        } else if (f_reflected < worst_f) {
          chosen = &reflected;
          chosen_f = f_reflected;
        } else if (budget_left()) {
          for (size_t d = 0; d < n; d++) {
            candidate[d] = centroid[d] +
                           kContract * (vertices[worst_row + d] - centroid[d]);
          }
          const scalar_t f_contracted = evaluate(candidate);
          if (f_contracted < worst_f) {
            chosen = &candidate;
            chosen_f = f_contracted;
          } else {
            // shrink every other vertex towards the best and re-evaluate it
            for (size_t k = 1; k < n_vertices && budget_left(); k++) {
              const size_t row = order[k] * n;
              for (size_t d = 0; d < n; d++) {
                candidate[d] =
                    vertices[best_row + d] +
                    kShrink * (vertices[row + d] - vertices[best_row + d]);
              }
              const scalar_t value = evaluate(candidate);
              store_vertex(order[k], candidate, value);
            }
          }
        }
        if (chosen != nullptr && accept(chosen_f, worst_f)) {
          store_vertex(worst_id, *chosen, chosen_f);
          if (chosen_f < worst_f) improved = true;
        }
      } else if (move < 2 * kMoveShare) {
        // ---- Differential Evolution: rand/1 or current-to-best/1 ----
        // three vertices from a random permutation of the sorted positions
        std::iota(shuffled.begin(), shuffled.end(), 0);
        for (size_t i = n_vertices - 1; i >= 1; i--) {
          const size_t j = generate_index(i + 1, this->generator);
          std::swap(shuffled[i], shuffled[j]);
        }
        const size_t row_a = order[shuffled[0]] * n;
        const size_t row_b = order[shuffled[1 % n_vertices]] * n;
        const size_t row_c = order[shuffled[2 % n_vertices]] * n;
        // the target is the worst vertex
        if (this->generator() < kRandMutantShare) {
          for (size_t d = 0; d < n; d++) {
            candidate[d] =
                vertices[row_a + d] +
                kDiffWeight * (vertices[row_b + d] - vertices[row_c + d]);
          }
        } else {
          for (size_t d = 0; d < n; d++) {
            candidate[d] =
                vertices[worst_row + d] +
                kDiffWeight *
                    (vertices[best_row + d] - vertices[worst_row + d]) +
                kDiffWeight * (vertices[row_b + d] - vertices[row_c + d]);
          }
        }
        // binomial crossover; one coordinate always comes from the mutant
        const size_t forced = generate_index(n, this->generator);
        for (size_t d = 0; d < n; d++) {
          const bool from_mutant =
              this->generator() < kCrossover || d == forced;
          if (!from_mutant) candidate[d] = vertices[worst_row + d];
        }
        const scalar_t f_trial = evaluate(candidate);
        if (accept(f_trial, worst_f)) {
          store_vertex(worst_id, candidate, f_trial);
          if (f_trial < worst_f) improved = true;
        }
      } else if (move < 3 * kMoveShare) {
        // ---- Gaussian sample around the best vertex, diagonal covariance ----
        // the unclipped offset feeds the covariance update
        for (size_t d = 0; d < n; d++) {
          const scalar_t z = rnorm<scalar_t>(this->generator);
          delta[d] = sigma * (z * std::sqrt(cov_diag[d]));
          candidate[d] = vertices[best_row + d] + delta[d];
        }
        const scalar_t f_sample = evaluate(candidate);
        if (accept(f_sample, worst_f)) {
          if (f_sample < worst_f) {
            improved = true;
            // exponential moving average of the squared successful offsets
            for (size_t d = 0; d < n; d++) {
              cov_diag[d] = kCovMemory * cov_diag[d] +
                            kCovInnovation * (delta[d] * delta[d] + kCovFloor);
            }
          }
          store_vertex(worst_id, candidate, f_sample);
        }
      } else {
        // ---- Hooke-Jeeves probe from the best vertex ----
        // +step then -step along each coordinate, keeping every improvement
        std::copy(vertices.begin() + best_row, vertices.begin() + best_row + n,
                  current.begin());
        scalar_t current_f = vertex_f[best_id];
        const scalar_t base_f = current_f;
        for (size_t d = 0; d < n && budget_left(); d++) {
          probe = current;
          probe[d] += step;
          const scalar_t f_plus = evaluate(probe);
          if (f_plus < current_f) {
            current = probe;
            current_f = f_plus;
            continue;
          }
          if (!budget_left()) break;
          probe = current;
          probe[d] -= step;
          const scalar_t f_minus = evaluate(probe);
          if (f_minus < current_f) {
            current = probe;
            current_f = f_minus;
          }
        }
        if (current_f < base_f) {
          improved = true;
          if (accept(current_f, worst_f)) {
            store_vertex(worst_id, current, current_f);
          }
        } else {
          step *= kStepShrink;
          if (step < kStepMin) step = kStep0;
        }
      }

      // ---- global scale adaptation and cooling ----
      if (improved) {
        stagnation = 0;
        sigma = std::min(kSigmaMax, sigma * kSigmaGrow);
      } else {
        stagnation++;
        sigma = std::max(kSigmaMin, sigma * kSigmaDecay);
      }
      temp *= kCooling;
      if (temp < kTempMin) temp = kTempMin;

      // ---- restart when stuck: reheat, rebuild around the best point ----
      if (stagnation > stall_limit && budget_left()) {
        stagnation = 0;
        temp = std::max(temp0 * kReheatFrac, temp * kReheatMult);
        step = kStep0;
        sigma = kSigma0;
        // jitter around a snapshot of the best point: an improving jitter point
        // moves best_u, but the remaining vertices still scatter about the
        // point the restart began from, as in the reference
        current = best_u;
        store_vertex(0, current, best_f);
        for (n_filled = 1; n_filled < n_vertices && budget_left(); n_filled++) {
          for (size_t d = 0; d < n; d++) {
            const scalar_t u = this->generator();
            candidate[d] =
                current[d] + (-kRestartJitter + 2.0 * kRestartJitter * u);
          }
          const scalar_t value = evaluate(candidate);
          store_vertex(n_filled, candidate, value);
        }
      }
      iter++;
    }

    for (size_t d = 0; d < n; d++) x[d] = origin[d] + best_u[d] * width[d];
    // report the true objective value (best_f carries f_multiplier)
    return solver_status<scalar_t>(f_multiplier * best_f, iter, f_evals);
  }
};

// Composite: CMA-ES bonded to L-BFGS-B. CMA-ES finds the basin and adapts a
// covariance that, on a smooth objective, converges to a multiple of the
// inverse Hessian. Once its search distribution has collapsed, when the global
// step sigma * sqrt(max eigenvalue of C) falls below `switch_step` times the
// widest box side (times the initial step when unbounded), L-BFGS-B starts
// from the best sampled point with sigma^2 C as its initial inverse Hessian.
// The first quasi-Newton step is then close to a Newton step, and the linear
// endgame of CMA-ES, which costs it several hundred evaluations per two
// decades of precision, is replaced by superlinear convergence at 2n
// evaluations per gradient.
//
// The hand-off threshold trades nothing measurable: on the 24-problem suite
// over ten seeds, thresholds of 1e-3, 1e-2 and 3e-2 pass 208, 211 and 209 of
// 240 runs (CMA-ES alone 205), while on the smooth problems 1e-2 reaches 1e-10
// in half to two thirds fewer evaluations than 1e-3.
//
// Elite polish recovers from a wrong basin using samples already paid for.
// By the time CMA-ES collapses it has evaluated hundreds of points, many in
// other basins; the best `elites` of them that lie at least kEliteDistance
// times the box width from the incumbent and from each other each get a
// short, unseeded L-BFGS-B polish capped at kElitePolishEvals evaluations
// (unseeded because the CMA covariance describes the incumbent's basin, and
// seeding with it drops the success rate from 97 to 89 percent). On the same
// sweep eight elites pass 236 of 240 runs at about 760 evaluations per solved
// problem; three IPOP restarts pass 230 at about 2550, since every restart
// re-runs CMA-ES with a doubled population whether or not it was needed.
//
// Budget-aware hand-off. Within a tight budget CMA-ES does not collapse at all
// (on the humpday demos at 60 to 480 evaluations it did so in under 6 percent
// of runs), so without a cap it spends everything and the polish never runs.
// The global phase may therefore use at most 1 - `polish_share` of the budget
// remaining at the start of a cycle; the rest goes to the seeded polish and
// then the elites. The default third was fixed before any measurement: two
// thirds to find the basin, one third to descend it, and since an L-BFGS-B
// iteration costs about 2n + 2 evaluations a third of the budget buys at least
// two quasi-Newton iterations whenever the budget exceeds 12 (n + 1). Zero
// restores the uncapped behaviour.
//
// The polish follows what its budget can afford. L-BFGS-B with central
// differences costs 2n + 2 evaluations per iteration and with forward
// differences n + 2; a polish budget that buys fewer than
// `polish_min_iterations` central iterations uses forward differences, and one
// that buys fewer forward iterations than that uses a (1+1) evolution strategy
// with the one-fifth success rule instead, whose every evaluation is a
// candidate. Six iterations is the smallest count at which quasi-Newton
// behaviour shows; zero disables the rule and always uses central differences.
//
// A polish that stops without converging (its line search found no decrease
// consistent with the gradient, which is what happens on rugged or noisy
// objectives) hands whatever remains of its budget to the (1+1) strategy from
// the point it reached, so finite-difference gradients failing does not leave
// budget unspent.
//
// `trust_start` off treats x0 as uninformative: one CMA-ES population is drawn
// uniformly over the box (about x0 when unbounded), evaluated together with
// x0, and the first cycle starts from the best of them. On by default, since a
// user-supplied start is usually worth more than a random sample.
//
// With `max_restarts` > 0 the cycle additionally repeats from a fresh random
// mean, uniform in the box or Gaussian about the starting point when
// unbounded, with the population doubled on every restart (the IPOP rule).
// Restarts are off by default for the cost reason above. A polish whose line
// search fails on a non-smooth objective simply ends. The best point over all
// phases is reported; `iteration` in the status counts CMA-ES cycles.
// Two behaviours of Alloy, abstracted rather than copied. Alloy's moves change
// their scale in a single evaluation (a Nelder-Mead contraction halves the
// simplex, a failed Hooke-Jeeves probe halves its step) where the (1+1)-CMA-ES
// step-size rule needs about 10 (1 + n / 2) evaluations per decade, the whole
// global phase at budgets of tens. With `budget_damping` (default) the
// steady-state phase uses OnePlusOneCMAES::budget_damping, the reference
// damping capped so that a decade of adaptation costs at most half the phase;
// budgets of thousands are unaffected. On the humpday demos this took the
// default from 41 to 44 percent of instances won against Alloy, with the gain
// at every budget; a fixed damping of one (the other extreme) won at 60 and
// 120 evaluations and lost at 240 and 480, which is the transient against
// stationary trade-off the rule resolves. Alloy also exploits a successful
// direction at once (the expansion and the pattern move); `pattern_moves`
// doubles a successful displacement while that keeps improving, in the
// steady-state phase and the strategy polish. Measured worse (7.3 against 7.0
// mean rank in an 18-solver race, and no gain on the suite) because on a bowl
// the doubling overshoots and on a rugged objective it lands in another bump,
// so it is off by default. The untrusted start takes its initial step from
// the spread of the design's best n + 1 points, the simplex a direct-search
// method would build; on the suite that lifts the untrusted start to 98.3
// percent solved at 8 percent more cost, on the demos it remains behind the
// trusted start.
// The global phase. Auto runs the population CMA-ES when its share of the
// budget buys at least kCovarianceHorizons times the covariance update's time
// constant, 1 / (c1 + cmu) generations, so that it can learn the problem's
// shape, and pattern search otherwise; the other values force one method.
// Below that horizon no covariance can be learnt, an evolution strategy is an
// isotropic random walk whose efficiency falls with dimension, and the
// coordinate polls of pattern search do not: on the humpday demos at 60 to
// 480 evaluations a bare compass search outranked every stochastic method
// (mean rank 4.6 of 13 against Alloy's 5.7 and 6.7 for the steady-state
// phase), rotated disguise included. The pattern phase hands over when its
// step falls below the switch threshold, the same collapse criterion as the
// distribution's, with no covariance to seed the polish; the strategy polish
// starts at the final step. SteadyState forces the (1+1)-CMA-ES (mean rank
// 5.1 against 5.5 for the population method in a 13-solver race; on the suite
// at 10000 evaluations 91.7 percent solved against 97.5, which the horizon
// rule separates). AlloyGlobal runs Alloy as the global phase: the best race
// result before the pattern phase and the second worst cost per solved problem
// on the suite, since Alloy never stops early; a measurement option, not a
// design.
enum class CompositeGlobal {
  Auto,
  Population,
  SteadyState,
  Pattern,
  AlloyGlobal
};
// The polish. Gradient (default): central differences when the budget affords
// polish_min_iterations of them, forward differences when it affords that many
// of those, the (1+1) strategy otherwise. Ladder inserts a quadratic model
// with a diagonal Hessian, fitted to the nearest archived points and driven by
// a trust region, between the two difference stencils; ModelFirst uses the
// model whenever the archive supports the fit. Both measured worse than
// Gradient on every budget and problem family of the humpday race except the
// rugged games, and ModelFirst loses 15 points of success on the suite, so
// they remain as options only.
enum class CompositePolish { Gradient, Ladder, ModelFirst };

template <typename Callable, typename RNG, typename scalar_t = double>
class [[maybe_unused]] Composite {
 public:
  // where the last run spent its evaluations, and which methods it used
  struct report {
    size_t global_evals = 0, polish_evals = 0, elite_evals = 0;
    bool steady_state = false;  // the (1+1)-CMA-ES ran the global phase
    bool pattern = false;       // pattern search ran the global phase
    char polish_method = '-';   // c central, f forward, m model, s strategy,
                                // p pattern search
    bool fallback = false;      // derivative-free polish took over after a
                                // gradient polish gave up
  };
  [[maybe_unused]] const report &exit_report() const { return report_; }

 private:
  Callable &f;
  RNG &generator;
  const scalar_t m_step;
  const size_t max_evals, elites, max_restarts;
  const scalar_t switch_step, grad_eps, polish_share;
  const bool trust_start;
  const size_t polish_min_iterations;
  const CompositeGlobal global_choice;
  const CompositePolish polish_choice;
  const bool budget_damping, pattern_moves, pattern_polish;
  report report_;
  // Auto runs the population method only when it affords this many covariance
  // time constants; below that the steady-state (1+1)-CMA-ES runs
  static constexpr scalar_t kCovarianceHorizons = 10.0;
  // quadratic model: coefficients 2n + 1, fitted to kModelPointsPerCoef times
  // as many nearest archived points; trust region grows by kTrustGrow above a
  // ratio of kTrustGood and shrinks by kTrustShrink below kTrustPoor
  static constexpr size_t kModelPointsPerCoef = 2;
  static constexpr scalar_t kTrustGrow = 2.0, kTrustShrink = 0.5,
                            kTrustGood = 0.75, kTrustPoor = 0.25;
  // per-polish caps: iterations and stored curvature pairs of L-BFGS-B
  static constexpr size_t kPolishIterations = 200, kPolishHistory = 10;
  // (1+1) evolution strategy: the step grows by kEsGrow on success and shrinks
  // by kEsGrow^(-1/4) on failure, which holds the success rate near one fifth
  static constexpr scalar_t kEsGrow = 1.5;
  // pattern move: a successful displacement is doubled while that improves
  static constexpr scalar_t kExpand = 2.0;
  // elite polish: minimum distance from the incumbent and from each other as a
  // fraction of the box width, and the evaluation cap of each elite's polish.
  // On the 24-problem sweep a cap of 150 solved every run a cap of 300 did,
  // at 15 percent less total cost.
  static constexpr scalar_t kEliteDistance = 0.05;
  static constexpr size_t kElitePolishEvals = 150;
  // CMA-ES stops on an ill-conditioned covariance or a collapsed cost spread
  // as it does standalone; the hand-off itself is driven by the step size
  static constexpr scalar_t kCondition = 1e14, kCostSpread = 1e-12;

 public:
  [[maybe_unused]] Composite<Callable, RNG, scalar_t>(
      Callable &f, RNG &generator, const scalar_t m_step = 0.5,
      const size_t max_evals = 10000, const size_t elites = 8,
      const size_t max_restarts = 0, const scalar_t switch_step = 1e-2,
      const scalar_t grad_eps = 1e-7, const scalar_t polish_share = 1.0 / 3.0,
      const bool trust_start = true, const size_t polish_min_iterations = 6,
      const CompositeGlobal global_choice = CompositeGlobal::Auto,
      const CompositePolish polish_choice = CompositePolish::Gradient,
      const bool budget_damping = true, const bool pattern_moves = false,
      const bool pattern_polish = true)
      : f(f),
        generator(generator),
        m_step(m_step),
        max_evals(max_evals),
        elites(elites),
        max_restarts(max_restarts),
        switch_step(switch_step),
        grad_eps(grad_eps),
        polish_share(polish_share),
        trust_start(trust_start),
        polish_min_iterations(polish_min_iterations),
        global_choice(global_choice),
        polish_choice(polish_choice),
        budget_damping(budget_damping),
        pattern_moves(pattern_moves),
        pattern_polish(pattern_polish) {}
  // minimize interface
  [[maybe_unused]] solver_status<scalar_t> minimize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;  // empty bounds => unconstrained
    return this->solve<true, false>(x, lower, upper);
  }
  // maximize interface
  [[maybe_unused]] solver_status<scalar_t> maximize(std::vector<scalar_t> &x) {
    std::vector<scalar_t> lower, upper;  // empty bounds => unconstrained
    return this->solve<false, false>(x, lower, upper);
  }
  // minimize with box constraints; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> minimize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<true, true>(x, lower, upper);
  }
  // maximize with box constraints; bounds are passed as (lower, upper)
  [[maybe_unused]] solver_status<scalar_t> maximize(
      std::vector<scalar_t> &x, const std::vector<scalar_t> &lower,
      const std::vector<scalar_t> &upper) {
    return this->solve<false, true>(x, lower, upper);
  }

 private:
  template <const bool minimize = true, const bool constrained = false>
  solver_status<scalar_t> solve(std::vector<scalar_t> &x,
                                const std::vector<scalar_t> &lower,
                                const std::vector<scalar_t> &upper) {
    // both phases minimize this signed, counted objective
    constexpr scalar_t sign = minimize ? 1.0 : -1.0;
    const size_t n = x.size();
    size_t evals = 0;
    // every point evaluated during a CMA-ES phase is archived (flat, row per
    // point) so the elite polish can revisit other basins without new samples
    bool archiving = false;
    std::vector<scalar_t> archive_x, archive_f;
    auto f_inner = [&](std::vector<scalar_t> &v) -> scalar_t {
      evals++;
      const scalar_t value = sign * this->f(v);
      if (archiving) {
        archive_x.insert(archive_x.end(), v.begin(), v.end());
        archive_f.push_back(value);
      }
      return value;
    };
    using Inner = decltype(f_inner);
    using Global = CMAES<Inner, RNG, scalar_t>;
    auto value_of = [](const solver_status<scalar_t> &st) {
      return std::get<2>(st.get_summary());
    };
    auto distance = [](const std::vector<scalar_t> &a,
                       const std::vector<scalar_t> &b) {
      scalar_t d2 = 0.0;
      for (size_t i = 0; i < a.size(); i++) d2 += (a[i] - b[i]) * (a[i] - b[i]);
      return std::sqrt(d2);
    };
    // the hand-off step: a fraction of the widest box side, or of the initial
    // step when there is no box
    scalar_t scale = this->m_step;
    if constexpr (constrained) {
      scale = 0.0;
      for (size_t i = 0; i < n; i++)
        scale = std::max(scale, upper[i] - lower[i]);
    }
    const scalar_t x_delta = this->switch_step * scale;
    // per-coordinate length unit for the model's distances: the box side, or
    // the hand-off scale when unbounded
    std::vector<scalar_t> coord_scale(n, scale);
    if constexpr (constrained) {
      for (size_t i = 0; i < n; i++) coord_scale[i] = upper[i] - lower[i];
    }
    report_ = report{};

    // L-BFGS-B from `xl` with the given gradient functor, seed (empty for the
    // identity) and evaluation cap; returns the value reached and whether the
    // projected gradient converged
    auto run_local = [&](auto grad, std::vector<scalar_t> &xl,
                         const std::vector<scalar_t> &seed, const size_t cap) {
      LBFGSB<Inner, scalar_t, decltype(grad)> local(
          f_inner, grad, kPolishIterations, this->grad_eps, kPolishHistory,
          seed, cap);
      solver_status<scalar_t> sl(0.0, 0, 0);
      if constexpr (constrained) {
        sl = local.minimize(xl, lower, upper);
      } else {
        sl = local.minimize(xl);
      }
      return std::make_pair(value_of(sl), sl.success());
    };
    // (1+1) evolution strategy from `xl` (value `fx`) with initial step
    // `sigma`, within `cap` evaluations; returns the value reached
    auto one_plus_one = [&](std::vector<scalar_t> &xl, scalar_t fx,
                            const size_t cap, scalar_t sigma) {
      const scalar_t shrink = std::pow(kEsGrow, -0.25);
      // converged once the step is below sqrt(eps) relative to the iterate,
      // the finite-difference resolution: at a point that is already optimal
      // the failures shrink sigma geometrically, so this ends the strategy
      // within a few hundred evaluations instead of the whole cap
      const scalar_t root_eps =
          std::sqrt(std::numeric_limits<scalar_t>::epsilon());
      std::vector<scalar_t> cand(n), disp(n);
      size_t used = 0;
      while (used < cap && evals < this->max_evals) {
        scalar_t x_inf = 0.0;
        for (size_t d = 0; d < n; d++) x_inf = std::max(x_inf, std::abs(xl[d]));
        if (sigma < root_eps * (1.0 + x_inf)) break;
        for (size_t d = 0; d < n; d++) {
          cand[d] = xl[d] + sigma * rnorm<scalar_t>(this->generator);
          if constexpr (constrained) {
            cand[d] = std::clamp(cand[d], lower[d], upper[d]);
          }
        }
        const scalar_t fc = f_inner(cand);
        used++;
        if (fc >= fx) {
          sigma *= shrink;
          continue;
        }
        for (size_t d = 0; d < n; d++) disp[d] = cand[d] - xl[d];
        fx = fc;
        xl = cand;
        sigma *= kEsGrow;
        // pattern move: double the displacement while it keeps improving;
        // these evaluations do not enter the one-fifth rule
        while (this->pattern_moves && used < cap && evals < this->max_evals) {
          bool moved = false;
          for (size_t d = 0; d < n; d++) {
            cand[d] = xl[d] + disp[d];
            if constexpr (constrained) {
              cand[d] = std::clamp(cand[d], lower[d], upper[d]);
            }
            moved = moved || cand[d] != xl[d];
          }
          if (!moved) break;
          const scalar_t fe = f_inner(cand);
          used++;
          if (fe >= fx) break;
          for (size_t d = 0; d < n; d++) disp[d] = cand[d] - xl[d] + disp[d];
          fx = fe;
          xl = cand;
        }
      }
      return fx;
    };
    // Quadratic model with a diagonal Hessian, m(x + d) = c + g'd + sum_i
    // h_i d_i^2 / 2, its 2n + 1 coefficients fitted by least squares to the
    // nearest archived points (distances in coord_scale units), driven by a
    // trust region: one evaluation per step, at the model minimiser within the
    // region and the box. The gradient and curvature come from samples already
    // paid for, so a step costs one evaluation instead of the 2n of a
    // finite-difference gradient. Returns the value reached; false in .second
    // when the archive is too thin to fit.
    const size_t model_coefs = 2 * n + 1;
    const size_t model_points = kModelPointsPerCoef * model_coefs;
    auto model_polish = [&](std::vector<scalar_t> &xl, scalar_t fx,
                            const size_t cap) -> std::pair<scalar_t, bool> {
      if (archive_f.size() < model_points) return {fx, false};
      const scalar_t root_eps =
          std::sqrt(std::numeric_limits<scalar_t>::epsilon());
      std::vector<std::pair<scalar_t, size_t>> nearest;
      std::vector<scalar_t> X(model_points * model_coefs), yv(model_points),
          cand(n), step(n);
      archiving = true;  // new points join the archive for the next fit
      scalar_t radius = 0.0;
      bool have_radius = false;
      for (size_t used = 0; used < cap && evals < this->max_evals; used++) {
        // the k nearest archived points, in scaled distance
        nearest.clear();
        for (size_t i = 0; i < archive_f.size(); i++) {
          scalar_t d2 = 0.0;
          for (size_t j = 0; j < n; j++) {
            const scalar_t d = (archive_x[i * n + j] - xl[j]) / coord_scale[j];
            d2 += d * d;
          }
          nearest.emplace_back(d2, i);
        }
        std::partial_sort(nearest.begin(), nearest.begin() + model_points,
                          nearest.end());
        if (!have_radius) {
          radius = std::sqrt(nearest[model_points - 1].first);
          have_radius = true;
        }
        if (radius < root_eps) break;
        // design matrix, column-major as tinyqr::lm expects: 1, d_j, d_j^2 / 2
        for (size_t r = 0; r < model_points; r++) {
          const size_t id = nearest[r].second;
          X[0 * model_points + r] = 1.0;
          for (size_t j = 0; j < n; j++) {
            const scalar_t d = (archive_x[id * n + j] - xl[j]) / coord_scale[j];
            X[(1 + j) * model_points + r] = d;
            X[(1 + n + j) * model_points + r] = 0.5 * d * d;
          }
          yv[r] = archive_f[id];
        }
        const std::vector<scalar_t> beta = tinyqr::lm(X, yv);
        bool finite = true;
        for (const scalar_t b : beta) finite = finite && std::isfinite(b);
        if (!finite) break;
        // model minimiser per coordinate (in scaled units), then clipped to the
        // trust region: a Newton step where the curvature is positive, a step
        // to the region's edge downhill where it is not
        scalar_t step_norm2 = 0.0;
        for (size_t j = 0; j < n; j++) {
          const scalar_t g = beta[1 + j], h = beta[1 + n + j];
          step[j] = h > 0.0 ? -g / h : (g > 0.0 ? -radius : radius);
          step_norm2 += step[j] * step[j];
        }
        const scalar_t step_norm = std::sqrt(step_norm2);
        const scalar_t shrink = step_norm > radius ? radius / step_norm : 1.0;
        scalar_t predicted = 0.0;
        for (size_t j = 0; j < n; j++) {
          step[j] *= shrink;
          predicted -=
              beta[1 + j] * step[j] + 0.5 * beta[1 + n + j] * step[j] * step[j];
          cand[j] = xl[j] + step[j] * coord_scale[j];
          if constexpr (constrained) {
            cand[j] = std::clamp(cand[j], lower[j], upper[j]);
          }
        }
        if (cand == xl) break;
        const scalar_t fc = f_inner(cand);
        const scalar_t actual = fx - fc;
        const scalar_t ratio =
            predicted > 0.0 ? actual / predicted : (actual > 0.0 ? 1.0 : -1.0);
        if (fc < fx) {
          fx = fc;
          xl = cand;
        }
        if (ratio > kTrustGood) {
          radius *= kTrustGrow;
        } else if (ratio < kTrustPoor) {
          radius *= kTrustShrink;
        }
      }
      archiving = false;
      return {fx, true};
    };
    // the polish of `xl` (value `fx`) within `cap` evaluations: the method the
    // cap can afford, see the class comment; `seed` is the L-BFGS-B initial
    // inverse Hessian (empty for the identity), `sigma` the strategy's step
    // The derivative-free polish: pattern search from the point reached, with
    // the given step, stopping at sqrt(eps) relative step as the strategy
    // does; or the (1+1) strategy when `pattern_polish` is off.
    auto derivative_free = [&](std::vector<scalar_t> &xl, const scalar_t fx,
                               const size_t cap, const scalar_t step) {
      if (!this->pattern_polish) {
        report_.polish_method = 's';
        return one_plus_one(xl, fx, cap, step);
      }
      report_.polish_method = 'p';
      if (cap == 0 || evals >= this->max_evals) return fx;
      scalar_t x_inf = 0.0;
      for (size_t d = 0; d < n; d++) x_inf = std::max(x_inf, std::abs(xl[d]));
      const scalar_t root_eps =
          std::sqrt(std::numeric_limits<scalar_t>::epsilon());
      PatternSearch<Inner, scalar_t> local(
          f_inner, step, std::min(cap, this->max_evals - evals),
          root_eps * (1.0 + x_inf));
      solver_status<scalar_t> st(fx, 0, 0);
      if constexpr (constrained) {
        st = local.minimize(xl, lower, upper);
      } else {
        st = local.minimize(xl);
      }
      return std::min(fx, value_of(st));
    };
    auto polish = [&](std::vector<scalar_t> &xl, const scalar_t fx,
                      const std::vector<scalar_t> &seed, const size_t cap,
                      const scalar_t sigma) -> scalar_t {
      const size_t central_iteration = 2 * n + 2, forward_iteration = n + 2;
      const size_t need = this->polish_min_iterations;
      const size_t before = evals;
      std::pair<scalar_t, bool> local(fx, false);
      const bool model_first =
          this->polish_choice == CompositePolish::ModelFirst &&
          archive_f.size() >= model_points;
      if (model_first) {
        local = model_polish(xl, fx, cap);
        report_.polish_method = 'm';
      } else if (need == 0 || cap >= need * central_iteration) {
        local = run_local(fin_diff<Inner, scalar_t, 0>{}, xl, seed, cap);
        report_.polish_method = 'c';
      } else if (this->polish_choice != CompositePolish::Gradient &&
                 archive_f.size() >= model_points) {
        local = model_polish(xl, fx, cap);
        report_.polish_method = 'm';
      } else if (cap >= need * forward_iteration) {
        local = run_local(fin_diff_fwd<Inner, scalar_t>{}, xl, seed, cap);
        report_.polish_method = 'f';
      } else {
        return derivative_free(xl, fx, cap, sigma);
      }
      const size_t spent = evals - before;
      if (local.second || spent >= cap) return local.first;
      // the gradient or model method gave up with budget left: continue
      // derivative-free
      report_.fallback = true;
      return derivative_free(xl, local.first, cap - spent, sigma);
    };

    const std::vector<scalar_t> x0 = x;
    std::vector<scalar_t> start = x, best_x = x, h0(n * n);
    scalar_t best_f = std::numeric_limits<scalar_t>::infinity();
    // the global phase's initial step; the untrusted start replaces it by the
    // spread of the design's best n + 1 points
    scalar_t start_step = this->m_step;
    size_t cycles = 0, pop_mult = 1;
    if (this->max_evals == 0) {
      // nothing was evaluated, so there is no valid result
      return solver_status<scalar_t>(std::numeric_limits<scalar_t>::quiet_NaN(),
                                     0, 0, 0, 0, false);
    }
    if (this->max_evals < Global::default_lambda(n) + 1) {
      // too small a budget for one CMA-ES generation: polish the start alone
      const scalar_t f0 = f_inner(x);
      const scalar_t fl =
          polish(x, f0, {}, this->max_evals - evals, kEliteDistance * scale);
      return solver_status<scalar_t>(sign * std::min(f0, fl), 1, evals);
    }
    if (!this->trust_start) {
      // ---- x0 carries no information: draw one population uniformly, keep
      // the best of it and x0 as the first mean ----
      archiving = true;
      best_f = f_inner(x);
      best_x = x;
      std::vector<scalar_t> cand(n);
      const size_t lambda0 = Global::default_lambda(n);
      for (size_t k = 0; k < lambda0 && evals < this->max_evals; k++) {
        for (size_t i = 0; i < n; i++) {
          if constexpr (constrained) {
            cand[i] = lower[i] + this->generator() * (upper[i] - lower[i]);
          } else {
            cand[i] = x0[i] + this->m_step * rnorm<scalar_t>(this->generator);
          }
        }
        const scalar_t fc = f_inner(cand);
        if (fc < best_f) {
          best_f = fc;
          best_x = cand;
        }
      }
      archiving = false;
      start = best_x;
      // ---- the initial step is the RMS coordinate distance from the best
      // point to the other members of the design's best n + 1 (the simplex a
      // direct-search method would build from them): a design whose best
      // points cluster starts a focused search, a spread one a wide search ----
      std::vector<size_t> by_value(archive_f.size());
      std::iota(by_value.begin(), by_value.end(), 0);
      std::sort(by_value.begin(), by_value.end(),
                [&](const size_t a, const size_t b) {
                  return archive_f[a] < archive_f[b];
                });
      const size_t simplex = std::min(n + 1, by_value.size());
      scalar_t spread2 = 0.0;
      for (size_t r = 1; r < simplex; r++) {
        const size_t id = by_value[r];
        for (size_t i = 0; i < n; i++) {
          const scalar_t d = archive_x[id * n + i] - best_x[i];
          spread2 += d * d;
        }
      }
      if (simplex > 1 && spread2 > 0.0) {
        start_step =
            std::sqrt(spread2 / static_cast<scalar_t>((simplex - 1) * n));
      }
    }
    while (true) {
      // ---- global phase: CMA-ES until its step size collapses, or until it
      // has spent its share of the remaining budget ----
      const size_t lambda = Global::default_lambda(n) * pop_mult;
      if (evals + lambda + 1 > this->max_evals) break;  // not one generation
      const auto remaining = static_cast<scalar_t>(this->max_evals - evals);
      const auto global_evals = static_cast<size_t>(
          std::floor((1.0 - this->polish_share) * remaining));
      // CMA-ES spends one evaluation on the mean plus lambda per generation
      const size_t generations =
          global_evals > lambda ? (global_evals - 1) / lambda : 1;
      const bool small_budget =
          static_cast<scalar_t>(generations) <
          kCovarianceHorizons * Global::covariance_horizon(n, pop_mult);
      const bool pattern =
          this->global_choice == CompositeGlobal::Pattern ||
          (this->global_choice == CompositeGlobal::Auto && small_budget);
      const bool steady = this->global_choice == CompositeGlobal::SteadyState;
      std::vector<scalar_t> xg = start;
      solver_status<scalar_t> sg(0.0, 0, 0);
      search_distribution<scalar_t> dist;
      bool have_dist = true;
      scalar_t pattern_step = 0.0;
      const size_t global_before = evals;
      archiving = true;
      if (this->global_choice == CompositeGlobal::AlloyGlobal) {
        Alloy<Inner, RNG, scalar_t> global(f_inner, this->generator,
                                           global_evals);
        if constexpr (constrained) {
          sg = global.minimize(xg, lower, upper);
        } else {
          sg = global.minimize(xg);
        }
        have_dist = false;
      } else if (pattern) {
        PatternSearch<Inner, scalar_t> global(f_inner, start_step, global_evals,
                                              x_delta);
        if constexpr (constrained) {
          sg = global.minimize(xg, lower, upper);
        } else {
          sg = global.minimize(xg);
        }
        have_dist = false;
        pattern_step = global.exit_step();
        report_.pattern = true;
      } else if (steady) {
        using Steady = OnePlusOneCMAES<Inner, RNG, scalar_t>;
        Steady global(
            f_inner, this->generator, start_step, global_evals, x_delta,
            this->budget_damping ? Steady::budget_damping(n, global_evals)
                                 : 0.0,
            this->pattern_moves);
        if constexpr (constrained) {
          sg = global.minimize(xg, lower, upper);
        } else {
          sg = global.minimize(xg);
        }
        dist = global.exit_distribution();
        report_.steady_state = true;
      } else {
        Global global(f_inner, this->generator, start_step, generations,
                      kCondition, x_delta, kCostSpread, pop_mult);
        if constexpr (constrained) {
          sg = global.minimize(xg, lower, upper);
        } else {
          sg = global.minimize(xg);
        }
        dist = global.exit_distribution();
      }
      archiving = false;
      report_.global_evals += evals - global_before;
      if (value_of(sg) < best_f) {
        best_f = value_of(sg);
        best_x = xg;
      }
      // ---- local phase: polish from the best sampled point, seeded with
      // sigma^2 C when a distribution exists; the strategy's step is the
      // distribution's RMS coordinate spread ----
      if (evals < this->max_evals) {
        // the strategy polish starts at the scale the global phase reached:
        // the pattern step, or the distribution's RMS coordinate spread
        scalar_t spread = pattern ? pattern_step : kEliteDistance * scale;
        std::vector<scalar_t> seed;
        if (have_dist) {
          const scalar_t sigma_sq = dist.sigma * dist.sigma;
          scalar_t trace = 0.0;
          for (size_t k = 0; k < n * n; k++) h0[k] = sigma_sq * dist.C[k];
          for (size_t d = 0; d < n; d++) trace += h0[d * n + d];
          spread = std::sqrt(trace / static_cast<scalar_t>(n));
          seed = h0;
        }
        std::vector<scalar_t> xl = xg;
        const size_t polish_before = evals;
        const scalar_t fl =
            polish(xl, value_of(sg), seed, this->max_evals - evals, spread);
        report_.polish_evals += evals - polish_before;
        if (fl < best_f) {
          best_f = fl;
          best_x = xl;
        }
      }
      // ---- elite polish (first cycle): other basins the samples visited ----
      const size_t elite_before = evals;
      if (cycles == 0 && this->elites > 0 && evals < this->max_evals) {
        std::vector<size_t> by_value(archive_f.size());
        std::iota(by_value.begin(), by_value.end(), 0);
        std::sort(by_value.begin(), by_value.end(),
                  [&](const size_t a, const size_t b) {
                    return archive_f[a] < archive_f[b];
                  });
        const scalar_t d_min = kEliteDistance * scale;
        // the incumbent and every elite already taken; a candidate must be at
        // least d_min from all of them
        std::vector<std::vector<scalar_t>> taken;
        taken.push_back(best_x);
        std::vector<scalar_t> xe(n);
        for (const size_t id : by_value) {
          if (taken.size() > this->elites || evals >= this->max_evals) break;
          std::copy(archive_x.begin() + id * n,
                    archive_x.begin() + (id + 1) * n, xe.begin());
          bool distinct = true;
          for (const auto &t : taken) {
            if (distance(xe, t) < d_min) {
              distinct = false;
              break;
            }
          }
          if (!distinct) continue;
          taken.push_back(xe);
          const scalar_t fe =
              polish(xe, archive_f[id], {},
                     std::min(kElitePolishEvals, this->max_evals - evals),
                     kEliteDistance * scale);
          if (fe < best_f) {
            best_f = fe;
            best_x = xe;
          }
        }
      }
      report_.elite_evals += evals - elite_before;
      cycles++;
      if (cycles > this->max_restarts || evals >= this->max_evals) break;
      // ---- IPOP restart: fresh mean, doubled population ----
      pop_mult *= 2;
      for (size_t i = 0; i < n; i++) {
        if constexpr (constrained) {
          start[i] = lower[i] + this->generator() * (upper[i] - lower[i]);
        } else {
          start[i] = x0[i] + this->m_step * rnorm<scalar_t>(this->generator);
        }
      }
    }
    x = best_x;
    // report the true objective value (best_f carries the sign)
    return solver_status<scalar_t>(sign * best_f, cycles, evals);
  }
};
};  // namespace nlsolver

#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#endif  // NLSOLVER_H_
