// BATS and TBATS exponential smoothing state space models.
//
// Single header, C++17. Depends on the standard library and on nlsolver.h
// (optimisation), which bundles tinyqr.h (linear algebra); both must be on
// the include path. DESIGN.md next to this file sets out the model, the
// state layout and the verification contract.
//
// [[CITATION]] De Livera, Hyndman and Snyder (2011), "Forecasting time series
// with complex seasonal patterns using exponential smoothing", Journal of the
// American Statistical Association 106(496), 1513-1527.
#ifndef TBATS_TBATS_H_
#define TBATS_TBATS_H_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "nlsolver.h"

namespace tbats {

using std::size_t;

// ===========================================================================
// 1. Dense matrix
// ===========================================================================

// Column-major dense matrix with its dimensions fixed at construction.
class Matrix {
 public:
  Matrix(const size_t rows, const size_t cols, const double fill = 0.0)
      : rows_(rows), cols_(cols), data_(rows * cols, fill) {}
  double &operator()(const size_t i, const size_t j) {
    return data_[j * rows_ + i];
  }
  double operator()(const size_t i, const size_t j) const {
    return data_[j * rows_ + i];
  }
  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  double *data() { return data_.data(); }
  const double *data() const { return data_.data(); }
  const std::vector<double> &storage() const { return data_; }
  std::vector<double> &storage() { return data_; }
  double *column(const size_t j) { return data_.data() + j * rows_; }
  const double *column(const size_t j) const { return data_.data() + j * rows_; }

 private:
  size_t rows_, cols_;
  std::vector<double> data_;
};

// ===========================================================================
// 2. Box-Cox transformation
// ===========================================================================

namespace box_cox {

inline double sign(const double x) { return (x > 0) - (x < 0); }

// y^(lambda) = (sign(y) |y|^lambda - 1) / lambda, log y for lambda = 0.
// A negative y with a negative lambda has no transform (NaN).
inline double transform(const double y, const double lambda) {
  if (lambda < 0 && y < 0) return std::numeric_limits<double>::quiet_NaN();
  if (lambda == 0) return std::log(y);
  return (sign(y) * std::pow(std::abs(y), lambda) - 1.0) / lambda;
}

// Inverse of transform(). For lambda < 0 the transform is bounded above by
// -1/lambda, so larger arguments have no inverse (NaN).
inline double inverse(const double z, const double lambda) {
  if (lambda < 0 && z > -1.0 / lambda) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (lambda == 0) return std::exp(z);
  const double zz = z * lambda + 1.0;
  return sign(zz) * std::pow(std::abs(zz), 1.0 / lambda);
}

// Back-transformed mean with the second order correction for the variance
// `variance` of z on the transformed scale:
//   E[y] ~ f(mu) (1 + variance (1 - lambda) / (2 f(mu)^(2 lambda))),
// f the inverse transform. For lambda = 0 this is the usual exp(mu)(1 + v/2).
inline double inverse_bias_adjusted(const double z, const double lambda,
                                    const double variance) {
  const double out = inverse(z, lambda);
  return out * (1.0 + 0.5 * variance * (1.0 - lambda) /
                          std::pow(out, 2.0 * lambda));
}

inline std::vector<double> transform(const std::vector<double> &y,
                                     const double lambda) {
  std::vector<double> out(y.size());
  for (size_t i = 0; i < y.size(); ++i) out[i] = transform(y[i], lambda);
  return out;
}

inline std::vector<double> inverse(const std::vector<double> &z,
                                   const double lambda) {
  std::vector<double> out(z.size());
  for (size_t i = 0; i < z.size(); ++i) out[i] = inverse(z[i], lambda);
  return out;
}

// Guerrero's criterion for a Box-Cox parameter: split the last floor(n /
// period) * period observations into consecutive blocks of `period`, and
// take the coefficient of variation across blocks of sd / mean^(1 - lambda).
// The lambda that makes this ratio constant across blocks stabilises the
// variance. Sample standard deviations (n - 1) throughout, as in R. A block
// with a negative mean has no real power and is left out, as R's na.rm does;
// the criterion is only meaningful for a positive series.
// [[CITATION]] Guerrero (1993), "Time-series analysis supported by power
// transformations", Journal of Forecasting 12(1), 37-48.
inline double guerrero_criterion(const std::vector<double> &y,
                                 const size_t period, const double lambda) {
  const size_t n = y.size();
  const size_t blocks = n / period;
  const size_t start = n - blocks * period;
  std::vector<double> ratio;
  ratio.reserve(blocks);
  for (size_t b = 0; b < blocks; ++b) {
    const double *block = y.data() + start + b * period;
    double mean = 0;
    for (size_t j = 0; j < period; ++j) mean += block[j];
    mean /= static_cast<double>(period);
    double ss = 0;
    for (size_t j = 0; j < period; ++j) {
      const double d = block[j] - mean;
      ss += d * d;
    }
    const double sd = std::sqrt(ss / static_cast<double>(period - 1));
    const double r = sd / std::pow(mean, 1.0 - lambda);
    if (!std::isnan(r)) ratio.push_back(r);
  }
  if (ratio.size() < 2) return std::numeric_limits<double>::quiet_NaN();
  double mean = 0;
  for (const double r : ratio) mean += r;
  mean /= static_cast<double>(ratio.size());
  double ss = 0;
  for (const double r : ratio) {
    const double d = r - mean;
    ss += d * d;
  }
  const double sd = std::sqrt(ss / static_cast<double>(ratio.size() - 1));
  return sd / mean;
}

// Tolerances of R's optimize(): the relative tolerance is the fourth root of
// the machine epsilon, and Brent's method adds the square root of the machine
// epsilon times |x|. The frozen Guerrero values were produced with these.
constexpr double kOptimizeTol = 1.220703125e-4;
constexpr double kOptimizeEps = 1.4901161193847656e-8;
constexpr size_t kOptimizeMaxIter = 1000;

// The Box-Cox parameter minimising Guerrero's criterion on [lower, upper].
// `period` is the seasonal period the blocks follow, at least 2. The lower
// bound is raised to zero when the series is not strictly positive, and a
// series shorter than two periods yields 1 (no transformation).
inline double guerrero_lambda(const std::vector<double> &y, size_t period,
                              double lower, const double upper) {
  period = std::max<size_t>(period, 2);
  for (const double v : y) {
    if (v <= 0) {
      lower = std::max(lower, 0.0);
      break;
    }
  }
  if (y.size() <= 2 * period) return 1.0;
  auto criterion = [&](const double lambda) {
    return guerrero_criterion(y, period, lambda);
  };
  nlsolver::Brent<decltype(criterion), double> brent(
      criterion, kOptimizeTol, kOptimizeEps, kOptimizeMaxIter);
  double lambda = 0.5 * (lower + upper);
  brent.minimize(lambda, lower, upper);
  return lambda;
}

}  // namespace box_cox

// ===========================================================================
// 3. Specification, state layout, parameters
// ===========================================================================

// How a seasonal component is represented: one state per position in the
// period (BATS), or two states per harmonic (TBATS).
enum class SeasonalType { kDummy, kTrigonometric };

struct SeasonalPeriod {
  double period;     // m_i; an integer for dummies, any value above 1 otherwise
  size_t harmonics;  // k_i; ignored for dummies
};

// Everything fixed before the parameters are optimised.
struct ModelSpec {
  bool box_cox = false;
  double box_cox_lower = 0.0;
  double box_cox_upper = 1.0;
  bool trend = false;
  bool damping = false;
  SeasonalType seasonal_type = SeasonalType::kTrigonometric;
  std::vector<SeasonalPeriod> seasonal;
  size_t p = 0;  // autoregressive order
  size_t q = 0;  // moving average order

  size_t seasonal_length(const size_t i) const {
    return seasonal_type == SeasonalType::kDummy
               ? static_cast<size_t>(seasonal[i].period)
               : 2 * seasonal[i].harmonics;
  }

  void validate() const {
    if (damping && !trend) {
      throw std::invalid_argument("tbats: damping requires a trend");
    }
    if (box_cox && !(box_cox_lower < box_cox_upper)) {
      throw std::invalid_argument("tbats: Box-Cox bounds must be ordered");
    }
    for (const auto &s : seasonal) {
      if (!(s.period > 1.0)) {
        throw std::invalid_argument("tbats: seasonal periods must exceed 1");
      }
      if (seasonal_type == SeasonalType::kDummy &&
          s.period != std::floor(s.period)) {
        throw std::invalid_argument(
            "tbats: dummy seasonal periods must be integers");
      }
      if (seasonal_type == SeasonalType::kTrigonometric && s.harmonics == 0) {
        throw std::invalid_argument(
            "tbats: a trigonometric seasonal period needs harmonics");
      }
    }
  }
};

// Position of every block in the state vector
//   x = (level, [trend], seasonal block 1, ..., seasonal block M,
//        d_t .. d_{t-p+1}, e_t .. e_{t-q+1}).
// Offsets are assigned by appending, here and nowhere else.
class StateLayout {
 public:
  explicit StateLayout(const ModelSpec &spec)
      : has_trend_(spec.trend),
        seasonal_offset_(offsets(spec)),
        seasonal_length_(lengths(spec)),
        ar_offset_(first_after_seasonal(spec)),
        p_(spec.p),
        ma_offset_(ar_offset_ + spec.p),
        q_(spec.q),
        dim_(ma_offset_ + spec.q) {}

  static constexpr size_t kLevel = 0;
  static constexpr size_t kTrend = 1;  // only when has_trend()

  bool has_trend() const { return has_trend_; }
  size_t seasonal_count() const { return seasonal_offset_.size(); }
  size_t seasonal_offset(const size_t i) const { return seasonal_offset_[i]; }
  size_t seasonal_length(const size_t i) const { return seasonal_length_[i]; }
  size_t ar_offset() const { return ar_offset_; }
  size_t p() const { return p_; }
  size_t ma_offset() const { return ma_offset_; }
  size_t q() const { return q_; }
  size_t dim() const { return dim_; }

 private:
  static size_t first_seasonal(const ModelSpec &spec) {
    return 1 + (spec.trend ? 1 : 0);
  }
  static std::vector<size_t> offsets(const ModelSpec &spec) {
    std::vector<size_t> out(spec.seasonal.size());
    size_t next = first_seasonal(spec);
    for (size_t i = 0; i < spec.seasonal.size(); ++i) {
      out[i] = next;
      next += spec.seasonal_length(i);
    }
    return out;
  }
  static std::vector<size_t> lengths(const ModelSpec &spec) {
    std::vector<size_t> out(spec.seasonal.size());
    for (size_t i = 0; i < spec.seasonal.size(); ++i) {
      out[i] = spec.seasonal_length(i);
    }
    return out;
  }
  static size_t first_after_seasonal(const ModelSpec &spec) {
    size_t next = first_seasonal(spec);
    for (size_t i = 0; i < spec.seasonal.size(); ++i) {
      next += spec.seasonal_length(i);
    }
    return next;
  }

  bool has_trend_;
  std::vector<size_t> seasonal_offset_, seasonal_length_;
  size_t ar_offset_, p_, ma_offset_, q_, dim_;
};

// The estimated quantities. A part the specification does not have keeps its
// neutral value and is never read: no damping is phi = 1, no trend leaves beta
// unused, no Box-Cox leaves lambda unused.
struct Parameters {
  double lambda = 1.0;  // Box-Cox parameter
  double alpha = 0.0;   // level smoothing
  double beta = 0.0;    // trend smoothing
  double phi = 1.0;     // damping
  std::vector<double> gamma_one;  // seasonal smoothing, one per period
  std::vector<double> gamma_two;  // second seasonal smoothing (harmonics only)
  std::vector<double> ar;
  std::vector<double> ma;
};

// Number of optimised parameters.
inline size_t parameter_count(const ModelSpec &spec) {
  const size_t gammas = spec.seasonal.size() *
                        (spec.seasonal_type == SeasonalType::kDummy ? 1 : 2);
  return (spec.box_cox ? 1 : 0) + 1 + (spec.damping ? 1 : 0) +
         (spec.trend ? 1 : 0) + gammas + spec.p + spec.q;
}

// Packed order: [lambda] alpha [phi] [beta] gamma_one... [gamma_two...] ar...
// ma..., the order of the forecast package, so that its parameter vectors
// compare directly.
inline std::vector<double> pack(const ModelSpec &spec, const Parameters &par) {
  std::vector<double> v;
  v.reserve(parameter_count(spec));
  if (spec.box_cox) v.push_back(par.lambda);
  v.push_back(par.alpha);
  if (spec.damping) v.push_back(par.phi);
  if (spec.trend) v.push_back(par.beta);
  v.insert(v.end(), par.gamma_one.begin(), par.gamma_one.end());
  if (spec.seasonal_type == SeasonalType::kTrigonometric) {
    v.insert(v.end(), par.gamma_two.begin(), par.gamma_two.end());
  }
  v.insert(v.end(), par.ar.begin(), par.ar.end());
  v.insert(v.end(), par.ma.begin(), par.ma.end());
  return v;
}

// Unpack into an existing Parameters whose vectors already have the right
// sizes, so that no memory is allocated.
inline void unpack_into(const ModelSpec &spec, const std::vector<double> &v,
                        Parameters &par) {
  if (v.size() != parameter_count(spec)) {
    throw std::invalid_argument("tbats: parameter vector has the wrong length");
  }
  size_t i = 0;
  if (spec.box_cox) par.lambda = v[i++];
  par.alpha = v[i++];
  if (spec.damping) par.phi = v[i++];
  if (spec.trend) par.beta = v[i++];
  const size_t m = spec.seasonal.size();
  par.gamma_one.assign(v.begin() + i, v.begin() + i + m);
  i += m;
  if (spec.seasonal_type == SeasonalType::kTrigonometric) {
    par.gamma_two.assign(v.begin() + i, v.begin() + i + m);
    i += m;
  }
  par.ar.assign(v.begin() + i, v.begin() + i + spec.p);
  i += spec.p;
  par.ma.assign(v.begin() + i, v.begin() + i + spec.q);
}

inline Parameters unpack(const ModelSpec &spec, const std::vector<double> &v) {
  Parameters par;
  unpack_into(spec, v, par);
  return par;
}

// Units in which the optimiser moves each packed parameter: the forecast
// package's parscale values, under which the frozen fits were produced.
// Harmonic seasonal smoothing is three orders of magnitude finer than dummy
// seasonal smoothing because its states feed back through rotations.
constexpr double kScaleLambda = 0.001;
constexpr double kScaleAlphaTrigonometric = 0.01;
constexpr double kScaleAlphaDummy = 0.1;
constexpr double kScalePhi = 0.01;
constexpr double kScaleBeta = 0.01;
constexpr double kScaleGammaTrigonometric = 1e-5;
constexpr double kScaleGammaDummy = 0.01;
constexpr double kScaleArma = 0.1;

inline std::vector<double> parameter_scales(const ModelSpec &spec) {
  const bool trig = spec.seasonal_type == SeasonalType::kTrigonometric;
  std::vector<double> s;
  s.reserve(parameter_count(spec));
  if (spec.box_cox) s.push_back(kScaleLambda);
  s.push_back(trig ? kScaleAlphaTrigonometric : kScaleAlphaDummy);
  if (spec.damping) s.push_back(kScalePhi);
  if (spec.trend) s.push_back(kScaleBeta);
  const size_t gammas = spec.seasonal.size() * (trig ? 2 : 1);
  s.insert(s.end(), gammas, trig ? kScaleGammaTrigonometric : kScaleGammaDummy);
  s.insert(s.end(), spec.p + spec.q, kScaleArma);
  return s;
}

// ===========================================================================
// 4. System and kernel
// ===========================================================================

// The innovations state space form
//   y_t = w' x_{t-1} + e_t,   x_t = F x_{t-1} + g e_t
// applied through its structure rather than through dense matrices. With the
// ARMA disturbance d_t = c' x_{t-1} + e_t (c holds the AR and MA coefficients
// in their blocks), the step is
//   level_t    = level + phi trend + alpha d_t
//   trend_t    = phi trend + beta d_t
//   harmonic j = rotation by 2 pi j / m of (s_j, s*_j) + (gamma1, gamma2) d_t
//   dummies    = cyclic shift, the entry leaving the back re-entering at the
//                front plus gamma d_t
//   AR block   = shift, d_t entering;  MA block = shift, e_t entering
// and the prediction w' x = level + phi trend + seasonal weights + c' x,
// where the seasonal weight of a harmonic block is the sum of its s_j and of
// a dummy block its last entry (the value one period back).
class System {
 public:
  System(const ModelSpec &spec, const Parameters &par)
      : spec_(validated(spec)),
        par_(par),
        layout_(spec),
        cos_(angles(spec, std::cos)),
        sin_(angles(spec, std::sin)) {
    check(par);
  }

  // Replace the parameters, keeping the specification, layout and angle
  // tables. The vectors of `par_` keep their sizes, so nothing is allocated.
  void set_parameters(const Parameters &par) {
    check(par);
    par_.lambda = par.lambda;
    par_.alpha = par.alpha;
    par_.beta = par.beta;
    par_.phi = par.phi;
    par_.gamma_one.assign(par.gamma_one.begin(), par.gamma_one.end());
    par_.gamma_two.assign(par.gamma_two.begin(), par.gamma_two.end());
    par_.ar.assign(par.ar.begin(), par.ar.end());
    par_.ma.assign(par.ma.begin(), par.ma.end());
  }

  const ModelSpec &spec() const { return spec_; }
  const Parameters &parameters() const { return par_; }
  const StateLayout &layout() const { return layout_; }
  size_t dim() const { return layout_.dim(); }

  // c' x: the ARMA part of the one-step prediction
  double arma_part(const double *x) const {
    double out = 0;
    const double *ar_state = x + layout_.ar_offset();
    for (size_t j = 0; j < layout_.p(); ++j) out += par_.ar[j] * ar_state[j];
    const double *ma_state = x + layout_.ma_offset();
    for (size_t j = 0; j < layout_.q(); ++j) out += par_.ma[j] * ma_state[j];
    return out;
  }

  // w' x: the one-step prediction on the model scale
  double predict(const double *x) const {
    double out = x[StateLayout::kLevel];
    if (layout_.has_trend()) out += par_.phi * x[StateLayout::kTrend];
    for (size_t i = 0; i < layout_.seasonal_count(); ++i) {
      const double *block = x + layout_.seasonal_offset(i);
      if (spec_.seasonal_type == SeasonalType::kTrigonometric) {
        for (size_t j = 0; j < spec_.seasonal[i].harmonics; ++j) {
          out += block[j];
        }
      } else {
        out += block[layout_.seasonal_length(i) - 1];
      }
    }
    return out + arma_part(x);
  }

  // out = F x + g e; `out` must not alias `x`
  void advance(const double *x, const double e, double *out) const {
    const double d = arma_part(x) + e;
    const double level = x[StateLayout::kLevel];
    const double trend = layout_.has_trend() ? x[StateLayout::kTrend] : 0.0;
    out[StateLayout::kLevel] = level + par_.phi * trend + par_.alpha * d;
    if (layout_.has_trend()) {
      out[StateLayout::kTrend] = par_.phi * trend + par_.beta * d;
    }
    size_t angle = 0;  // running index into cos_ / sin_
    for (size_t i = 0; i < layout_.seasonal_count(); ++i) {
      const size_t offset = layout_.seasonal_offset(i);
      const double *block = x + offset;
      double *next = out + offset;
      if (spec_.seasonal_type == SeasonalType::kTrigonometric) {
        const size_t k = spec_.seasonal[i].harmonics;
        const double g1 = par_.gamma_one[i] * d;
        const double g2 = par_.gamma_two[i] * d;
        for (size_t j = 0; j < k; ++j, ++angle) {
          const double s = block[j];
          const double s_star = block[k + j];
          next[j] = s * cos_[angle] + s_star * sin_[angle] + g1;
          next[k + j] = -s * sin_[angle] + s_star * cos_[angle] + g2;
        }
      } else {
        const size_t m = layout_.seasonal_length(i);
        next[0] = block[m - 1] + par_.gamma_one[i] * d;
        for (size_t r = 1; r < m; ++r) next[r] = block[r - 1];
      }
    }
    double *ar_next = out + layout_.ar_offset();
    const double *ar_state = x + layout_.ar_offset();
    if (layout_.p() > 0) {
      ar_next[0] = d;
      for (size_t j = 1; j < layout_.p(); ++j) ar_next[j] = ar_state[j - 1];
    }
    double *ma_next = out + layout_.ma_offset();
    const double *ma_state = x + layout_.ma_offset();
    if (layout_.q() > 0) {
      ma_next[0] = e;
      for (size_t j = 1; j < layout_.q(); ++j) ma_next[j] = ma_state[j - 1];
    }
  }

  // Dense w, g, F and D = F - g w', derived from the kernel: column i of F is
  // the step applied to the i-th unit vector with no innovation, g is the
  // step applied to the zero state with a unit innovation. The `_into`
  // forms write into buffers of size dim() (`unit` must hold zeros and is
  // left so); the others allocate.
  void observation_into(double *w, double *unit) const {
    for (size_t i = 0; i < dim(); ++i) {
      unit[i] = 1.0;
      w[i] = predict(unit);
      unit[i] = 0.0;
    }
  }
  void gain_into(double *g, const double *zero) const {
    advance(zero, 1.0, g);
  }
  void transition_into(Matrix &F, double *unit) const {
    for (size_t i = 0; i < dim(); ++i) {
      unit[i] = 1.0;
      advance(unit, 0.0, F.column(i));
      unit[i] = 0.0;
    }
  }
  void D_into(Matrix &D, double *w, double *g, double *unit) const {
    transition_into(D, unit);
    observation_into(w, unit);
    gain_into(g, unit);
    for (size_t j = 0; j < dim(); ++j) {
      for (size_t i = 0; i < dim(); ++i) D(i, j) -= g[i] * w[j];
    }
  }
  std::vector<double> observation() const {
    std::vector<double> w(dim()), unit(dim(), 0.0);
    observation_into(w.data(), unit.data());
    return w;
  }
  std::vector<double> gain() const {
    std::vector<double> g(dim()), zero(dim(), 0.0);
    gain_into(g.data(), zero.data());
    return g;
  }
  Matrix transition() const {
    Matrix F(dim(), dim());
    std::vector<double> unit(dim(), 0.0);
    transition_into(F, unit.data());
    return F;
  }
  Matrix D() const {
    Matrix out(dim(), dim());
    std::vector<double> w(dim()), g(dim()), unit(dim(), 0.0);
    D_into(out, w.data(), g.data(), unit.data());
    return out;
  }

 private:
  static const ModelSpec &validated(const ModelSpec &spec) {
    spec.validate();
    return spec;
  }
  // cos or sin of every harmonic angle 2 pi j / m_i, blocks concatenated
  static std::vector<double> angles(const ModelSpec &spec,
                                    double (*fn)(double)) {
    std::vector<double> out;
    if (spec.seasonal_type != SeasonalType::kTrigonometric) return out;
    constexpr double kTwoPi = 6.283185307179586476925286766559;
    for (const auto &s : spec.seasonal) {
      for (size_t j = 1; j <= s.harmonics; ++j) {
        out.push_back(fn(kTwoPi * static_cast<double>(j) / s.period));
      }
    }
    return out;
  }
  void check(const Parameters &par) const {
    const size_t m = spec_.seasonal.size();
    const bool trig = spec_.seasonal_type == SeasonalType::kTrigonometric;
    if (par.gamma_one.size() != m || (trig && par.gamma_two.size() != m) ||
        par.ar.size() != spec_.p || par.ma.size() != spec_.q) {
      throw std::invalid_argument(
          "tbats: parameters do not match the specification");
    }
  }

  ModelSpec spec_;
  Parameters par_;
  StateLayout layout_;
  std::vector<double> cos_, sin_;
};

// ---- kernel consumers -----------------------------------------------------

// Errors e_t = y_t - w' x_{t-1}, predictions w' x_{t-1} and the state history
// (column t holds x_t) of running the filter over the model-scale series `y`
// from the seed state `x0`.
struct FilterResult {
  FilterResult(const size_t n, const size_t dim)
      : errors(n), predictions(n), states(dim, n) {}
  std::vector<double> errors, predictions;
  Matrix states;
  std::vector<double> final_state() const {
    const size_t last = states.cols() - 1;
    return std::vector<double>(states.column(last),
                               states.column(last) + states.rows());
  }
};

inline FilterResult filter(const System &system, const std::vector<double> &y,
                           const std::vector<double> &x0) {
  if (x0.size() != system.dim()) {
    throw std::invalid_argument("tbats: seed state has the wrong dimension");
  }
  if (y.empty()) throw std::invalid_argument("tbats: empty series");
  FilterResult out(y.size(), system.dim());
  const double *previous = x0.data();
  for (size_t t = 0; t < y.size(); ++t) {
    out.predictions[t] = system.predict(previous);
    out.errors[t] = y[t] - out.predictions[t];
    system.advance(previous, out.errors[t], out.states.column(t));
    previous = out.states.column(t);
  }
  return out;
}

// Sum of squared one-step errors over the model-scale series `y` from `x0`,
// keeping no history: the likelihood's inner loop.
inline double sum_squared_errors(const System &system, const double *y,
                                 const size_t n, const double *x0,
                                 std::vector<double> &state,
                                 std::vector<double> &next) {
  state.assign(x0, x0 + system.dim());
  next.resize(system.dim());
  double sse = 0;
  for (size_t t = 0; t < n; ++t) {
    const double e = y[t] - system.predict(state.data());
    sse += e * e;
    system.advance(state.data(), e, next.data());
    state.swap(next);
  }
  return sse;
}

// Forecast means w' F^(h-1) x and variances sigma^2 (1 + sum_{j<h} c_j^2),
// c_j = w' F^(j-1) g, on the model scale, h steps from the state `x_last`.
struct ForecastModelScale {
  explicit ForecastModelScale(const size_t h) : mean(h), variance(h) {}
  std::vector<double> mean, variance;
};

inline ForecastModelScale forecast_model_scale(const System &system,
                                               const std::vector<double> &x_last,
                                               const double sigma2,
                                               const size_t h) {
  if (x_last.size() != system.dim()) {
    throw std::invalid_argument("tbats: state has the wrong dimension");
  }
  ForecastModelScale out(h);
  std::vector<double> x = x_last, next(system.dim());
  std::vector<double> v = system.gain();  // F^(j-1) g
  double multiplier = 1.0;                 // 1 + sum of c_j^2 so far
  for (size_t i = 0; i < h; ++i) {
    out.mean[i] = system.predict(x.data());
    out.variance[i] = sigma2 * multiplier;
    system.advance(x.data(), 0.0, next.data());
    x.swap(next);
    const double c = system.predict(v.data());
    multiplier += c * c;
    system.advance(v.data(), 0.0, next.data());
    v.swap(next);
  }
  return out;
}

// Sample paths on the model scale: one column per path of `innovations`
// (rows are steps), all starting from `x_last`.
inline Matrix simulate_model_scale(const System &system,
                                   const std::vector<double> &x_last,
                                   const Matrix &innovations) {
  if (x_last.size() != system.dim()) {
    throw std::invalid_argument("tbats: state has the wrong dimension");
  }
  const size_t h = innovations.rows(), paths = innovations.cols();
  Matrix out(h, paths);
  std::vector<double> x(system.dim()), next(system.dim());
  for (size_t path = 0; path < paths; ++path) {
    x = x_last;
    for (size_t t = 0; t < h; ++t) {
      const double e = innovations(t, path);
      out(t, path) = system.predict(x.data()) + e;
      system.advance(x.data(), e, next.data());
      x.swap(next);
    }
  }
  return out;
}

// ===========================================================================
// 5. Seed states
// ===========================================================================

namespace seed {

inline size_t gcd(size_t a, size_t b) {
  while (b != 0) {
    a %= b;
    std::swap(a, b);
  }
  return a;
}

// How many trailing columns of each dummy seasonal block leave the seed
// regression. The m dummy states of a block are only identified up to a
// constant (the level absorbs it), so one column goes; a period that divides
// a later period is entirely redundant with that block, so all of its columns
// go; two periods sharing a factor h describe the same h-periodic effect
// twice, so the later block loses h columns. The dropped states are restored
// afterwards so that every block sums to zero. The loop order follows the
// forecast package, where a later shared factor overrides an earlier one.
inline std::vector<size_t> dummy_columns_dropped(
    const std::vector<size_t> &periods) {
  const size_t count = periods.size();
  std::vector<size_t> dropped(count, 1);
  std::vector<bool> redundant(count, false);
  for (size_t i = count; i-- > 1;) {
    for (size_t j = 0; j < i; ++j) {
      if (periods[i] % periods[j] == 0) redundant[j] = true;
    }
  }
  for (size_t i = 0; i < count; ++i) {
    if (redundant[i]) dropped[i] = periods[i];
  }
  for (size_t s = count; s-- > 1;) {
    for (size_t j = s; j-- > 0;) {
      const size_t h = gcd(periods[s], periods[j]);
      if (h != 1 && !redundant[s] && !redundant[j]) dropped[s] = h;
    }
  }
  return dropped;
}

// Columns of the seed design that enter the regression: level, trend, every
// harmonic state, the leading (m - dropped) states of each dummy block, and
// no ARMA state (ARMA seeds are fixed at zero).
inline std::vector<size_t> kept_columns(const System &system,
                                        const std::vector<size_t> &dropped) {
  const StateLayout &layout = system.layout();
  std::vector<size_t> kept;
  kept.push_back(StateLayout::kLevel);
  if (layout.has_trend()) kept.push_back(StateLayout::kTrend);
  for (size_t i = 0; i < layout.seasonal_count(); ++i) {
    const size_t length = layout.seasonal_length(i);
    const size_t keep = dropped.empty() ? length : length - dropped[i];
    for (size_t j = 0; j < keep; ++j) {
      kept.push_back(layout.seasonal_offset(i) + j);
    }
  }
  return kept;
}

// Seed state x_0 by least squares. With x_0 = 0 the filter yields errors
// ytilde_t, and for any x_0 the errors are ytilde_t - wtilde_t' x_0 with
// wtilde_t' = w' D^(t-1), D = F - g w' (the recursion x_t = D x_{t-1} + g y_t
// is linear in x_0). The regression of ytilde on wtilde without intercept
// minimises the sum of squared errors over x_0. `y` is on the model scale.
inline std::vector<double> seed_states(const System &system,
                                       const std::vector<double> &y) {
  const size_t n = y.size(), dim = system.dim();
  const StateLayout &layout = system.layout();
  const std::vector<double> zero(dim, 0.0), w = system.observation();
  const FilterResult from_zero = filter(system, y, zero);
  const Matrix D = system.D();
  // design: row t = w' D^(t-1), column-major n by dim
  Matrix design(n, dim);
  for (size_t j = 0; j < dim; ++j) design(0, j) = w[j];
  for (size_t t = 1; t < n; ++t) {
    for (size_t j = 0; j < dim; ++j) {
      double acc = 0;
      for (size_t k = 0; k < dim; ++k) acc += design(t - 1, k) * D(k, j);
      design(t, j) = acc;
    }
  }
  std::vector<size_t> dropped;
  if (system.spec().seasonal_type == SeasonalType::kDummy) {
    std::vector<size_t> periods;
    for (const auto &s : system.spec().seasonal) {
      periods.push_back(static_cast<size_t>(s.period));
    }
    dropped = dummy_columns_dropped(periods);
  }
  const std::vector<size_t> kept = kept_columns(system, dropped);
  if (kept.size() >= n) {
    throw std::runtime_error("tbats: series too short to seed the states");
  }
  std::vector<double> X(n * kept.size());
  for (size_t c = 0; c < kept.size(); ++c) {
    std::copy(design.column(kept[c]), design.column(kept[c]) + n,
              X.begin() + c * n);
  }
  const std::vector<double> coef = tinyqr::lm(X, from_zero.errors);
  for (const double v : coef) {
    if (!std::isfinite(v)) {
      throw std::runtime_error("tbats: seed state regression is singular");
    }
  }
  std::vector<double> x0(dim, 0.0);
  size_t next = 0;
  x0[StateLayout::kLevel] = coef[next++];
  if (layout.has_trend()) x0[StateLayout::kTrend] = coef[next++];
  for (size_t i = 0; i < layout.seasonal_count(); ++i) {
    const size_t offset = layout.seasonal_offset(i);
    const size_t length = layout.seasonal_length(i);
    if (dropped.empty()) {
      for (size_t j = 0; j < length; ++j) x0[offset + j] = coef[next++];
    } else {
      // centre the block: the estimated states less their mean, the dropped
      // states at minus that mean, so the m states sum to zero
      const size_t estimated = length - dropped[i];
      double total = 0;
      for (size_t j = 0; j < estimated; ++j) total += coef[next + j];
      const double mean = total / static_cast<double>(length);
      for (size_t j = 0; j < estimated; ++j) x0[offset + j] = coef[next + j] - mean;
      for (size_t j = estimated; j < length; ++j) x0[offset + j] = -mean;
      next += estimated;
    }
  }
  return x0;
}

}  // namespace seed

// ===========================================================================
// 6. Admissibility, likelihood, optimisation
// ===========================================================================

namespace admissibility {

// Roots of the AR and MA polynomials must lie beyond 1 + kRootMargin and
// every eigenvalue of D within 1 + kRootMargin: a one percent margin on
// stationarity, invertibility and forecastability, as in the forecast
// package.
constexpr double kRootMargin = 0.01;
// Damping below this makes the trend die out within a few steps and the
// trend parameters unidentifiable.
constexpr double kPhiLower = 0.8;
// Trailing ARMA coefficients below this do not count toward the order of the
// polynomial whose roots are checked.
constexpr double kCoefficientEps = 1e-8;

// Buffers for the root checks of one specification, sized once.
struct RootScratch {
  explicit RootScratch(const ModelSpec &spec) {
    const size_t order = std::max(spec.p, spec.q);
    companion.reserve(order * order);
    coefficients.reserve(order);
    eigen.reserve(order);
  }
  std::vector<double> companion, coefficients;
  std::vector<std::complex<double>> eigen;
};

// True when every root of 1 + c_1 z + ... + c_p z^p has modulus at least
// 1 + kRootMargin. The reciprocals of the roots are the eigenvalues of the
// companion matrix of z^p + c_1 z^(p-1) + ... + c_p, so the condition is
// that no eigenvalue exceeds 1 / (1 + kRootMargin).
inline bool roots_outside(const std::vector<double> &c, RootScratch &scratch) {
  size_t order = 0;
  for (size_t i = 0; i < c.size(); ++i) {
    if (std::abs(c[i]) > kCoefficientEps) order = i + 1;
  }
  if (order == 0) return true;
  std::vector<double> &companion = scratch.companion;
  companion.assign(order * order, 0.0);  // column-major
  for (size_t j = 0; j < order; ++j) companion[j * order] = -c[j];
  for (size_t i = 1; i < order; ++i) companion[(i - 1) * order + i] = 1.0;
  try {
    tinyqr::eigenvalues(companion.data(), order, scratch.eigen);
  } catch (const std::runtime_error &) {
    return false;  // no convergence: nothing sensible lies here
  }
  const double bound = 1.0 / (1.0 + kRootMargin);
  for (const auto &v : scratch.eigen) {
    if (std::abs(v) > bound) return false;
  }
  return true;
}

inline bool all_finite(const std::vector<double> &v) {
  for (const double x : v) {
    if (!std::isfinite(x)) return false;
  }
  return true;
}

// The checks that need no filtering: finiteness, bounds and ARMA roots.
inline bool parameters_admissible(const ModelSpec &spec, const Parameters &par,
                                  RootScratch &scratch) {
  if (!std::isfinite(par.lambda) || !std::isfinite(par.alpha) ||
      !std::isfinite(par.beta) || !std::isfinite(par.phi) ||
      !all_finite(par.gamma_one) || !all_finite(par.gamma_two) ||
      !all_finite(par.ar) || !all_finite(par.ma)) {
    return false;
  }
  if (spec.box_cox &&
      (par.lambda <= spec.box_cox_lower || par.lambda >= spec.box_cox_upper)) {
    return false;
  }
  if (spec.trend && (par.phi < kPhiLower || par.phi > 1.0)) return false;
  std::vector<double> &ar = scratch.coefficients;
  ar.assign(par.ar.size(), 0.0);
  for (size_t i = 0; i < ar.size(); ++i) ar[i] = -par.ar[i];
  return roots_outside(ar, scratch) && roots_outside(par.ma, scratch);
}

inline bool parameters_admissible(const ModelSpec &spec, const Parameters &par) {
  RootScratch scratch(spec);
  return parameters_admissible(spec, par, scratch);
}

// Buffers for the dense matrices of one system, sized once.
struct SystemScratch {
  explicit SystemScratch(const size_t dim)
      : D(dim, dim), w(dim), g(dim), unit(dim, 0.0) {
    eigen.reserve(dim);
  }
  Matrix D;
  std::vector<double> w, g, unit;
  std::vector<std::complex<double>> eigen;
};

// Forecastability: the recursion x_t = D x_{t-1} + g y_t must not amplify.
inline bool system_admissible(const System &system, SystemScratch &scratch) {
  system.D_into(scratch.D, scratch.w.data(), scratch.g.data(),
                scratch.unit.data());
  try {
    tinyqr::eigenvalues(scratch.D.data(), system.dim(), scratch.eigen);
  } catch (const std::runtime_error &) {
    return false;  // no convergence: nothing sensible lies here
  }
  for (const auto &v : scratch.eigen) {
    if (!(std::abs(v) < 1.0 + kRootMargin)) return false;
  }
  return true;
}

inline bool system_admissible(const System &system) {
  SystemScratch scratch(system.dim());
  return system_admissible(system, scratch);
}

}  // namespace admissibility

// Value returned for an inadmissible or non-finite likelihood: large and
// finite, so the simplex statistics of Nelder-Mead stay defined (R's optim
// substitutes the same value for non-finite objectives).
constexpr double kInadmissibleValue = 1e35;

// When the forecastability check (the eigenvalues of D, cubic in the state
// dimension) runs during optimisation. A point that is worse than the best
// value seen so far can never become the returned optimum, so with
// kOnImprovement only improving points are verified before they are
// accepted as the new best and the returned optimum is admissible exactly
// as with kEveryEvaluation; the bounds and ARMA root checks always run.
// The two settings do not follow the same path, since with kOnImprovement a
// non-improving inadmissible point can sit in the simplex with its raw
// value, and on the frozen series that moved two searches to a worse local
// optimum. kAuto therefore checks every evaluation up to
// kForecastabilityAutoDim states and on improvement above.
enum class ForecastabilityCheck { kAuto, kEveryEvaluation, kOnImprovement };

// At 64 states the check costs about half a millisecond per evaluation, a
// second per fit; at 192 states (daily and weekly dummy periods together)
// about fifteen milliseconds, half a minute per fit, where checking only
// improving points is worth a different local optimum.
constexpr size_t kForecastabilityAutoDim = 64;

inline ForecastabilityCheck resolve_forecastability(const ForecastabilityCheck check,
                                                    const size_t dim) {
  if (check != ForecastabilityCheck::kAuto) return check;
  return dim > kForecastabilityAutoDim ? ForecastabilityCheck::kOnImprovement
                                       : ForecastabilityCheck::kEveryEvaluation;
}

// -2 log L up to constants, as a function of the packed parameters:
//   n log(sum e_t^2) - 2 (lambda - 1) sum log y_t,
// the second term only with a Box-Cox transformation (the Jacobian). The
// optimiser works on the parameters divided by parameter_scales().
//
// With a Box-Cox transformation the seed state is held on the original
// scale and transformed with the trial lambda at every evaluation, which is
// how the forecast package keeps one seed estimate valid across lambdas.
class Likelihood {
 public:
  Likelihood(const ModelSpec &spec, const std::vector<double> &y,
             const std::vector<double> &seed,
             const ForecastabilityCheck check = ForecastabilityCheck::kEveryEvaluation)
      : spec_(spec),
        check_(resolve_forecastability(check, StateLayout(spec).dim())),
        best_(std::numeric_limits<double>::infinity()),
        y_(y),
        seed_(seed),
        scales_(parameter_scales(spec)),
        log_sum_(spec.box_cox ? log_sum(y) : 0.0),
        system_(spec, unpack(spec, std::vector<double>(scales_.size(), 0.0))),
        par_(system_.parameters()),
        unscaled_(scales_.size()),
        roots_(spec),
        dense_(system_.dim()),
        y_model_(y.size()),
        x0_(seed.size()),
        state_(seed.size()),
        next_(seed.size()) {}

  size_t count() const { return scales_.size(); }
  const System &system() const { return system_; }
  std::vector<double> scale(const Parameters &par) const {
    std::vector<double> v = pack(spec_, par);
    for (size_t i = 0; i < v.size(); ++i) v[i] /= scales_[i];
    return v;
  }
  Parameters unscale(const std::vector<double> &scaled) const {
    std::vector<double> v = scaled;
    for (size_t i = 0; i < v.size(); ++i) v[i] *= scales_[i];
    return unpack(spec_, v);
  }
  // seed state on the model scale for these parameters
  const std::vector<double> &seed_for(const Parameters &par) {
    if (spec_.box_cox) {
      for (size_t i = 0; i < seed_.size(); ++i) {
        x0_[i] = box_cox::transform(seed_[i], par.lambda);
      }
    } else {
      x0_ = seed_;
    }
    return x0_;
  }

  // nlsolver's objective interface; allocates nothing once constructed
  double operator()(std::vector<double> &scaled) {
    unscaled_.assign(scaled.begin(), scaled.end());
    for (size_t i = 0; i < unscaled_.size(); ++i) unscaled_[i] *= scales_[i];
    unpack_into(spec_, unscaled_, par_);
    return evaluate(par_, false);
  }

  // The value at `par`, with the full admissibility check.
  double evaluate(const Parameters &par) { return evaluate(par, true); }

  double evaluate(const Parameters &par, const bool full_check) {
    if (!admissibility::parameters_admissible(spec_, par, roots_)) {
      return kInadmissibleValue;
    }
    system_.set_parameters(par);
    const System &system = system_;
    const double *y = y_.data();
    if (spec_.box_cox) {
      for (size_t t = 0; t < y_.size(); ++t) {
        y_model_[t] = box_cox::transform(y_[t], par.lambda);
      }
      y = y_model_.data();
    }
    const std::vector<double> &x0 = seed_for(par);
    const double sse = sum_squared_errors(system, y, y_.size(), x0.data(),
                                          state_, next_);
    double value = static_cast<double>(y_.size()) * std::log(sse);
    if (spec_.box_cox) value -= 2.0 * (par.lambda - 1.0) * log_sum_;
    if (!std::isfinite(value)) return kInadmissibleValue;
    const bool verify = full_check ||
                        check_ == ForecastabilityCheck::kEveryEvaluation ||
                        value < best_;
    if (verify && !admissibility::system_admissible(system, dense_)) {
      return kInadmissibleValue;
    }
    if (value < best_) best_ = value;
    return value;
  }

 private:
  static double log_sum(const std::vector<double> &y) {
    double out = 0;
    for (const double v : y) out += std::log(v);
    return out;
  }

  ModelSpec spec_;
  ForecastabilityCheck check_;
  double best_;  // best admissible value seen, for kOnImprovement
  std::vector<double> y_, seed_, scales_;
  double log_sum_;
  System system_;
  Parameters par_;
  std::vector<double> unscaled_;
  admissibility::RootScratch roots_;
  admissibility::SystemScratch dense_;
  std::vector<double> y_model_, x0_, state_, next_;
};

// ---- optimiser dispatch -----------------------------------------------------

enum class Optimizer { kNelderMead, kBFGS, kLBFGSB, kPatternSearch };

inline Optimizer optimizer_from_name(const std::string &name) {
  if (name == "nelder_mead") return Optimizer::kNelderMead;
  if (name == "bfgs") return Optimizer::kBFGS;
  if (name == "lbfgsb") return Optimizer::kLBFGSB;
  if (name == "pattern_search") return Optimizer::kPatternSearch;
  throw std::invalid_argument("tbats: unknown optimizer '" + name + "'");
}

// Initial simplex step of R's Nelder-Mead: a tenth of the largest scaled
// coordinate, or 0.1 when all are zero.
constexpr double kSimplexFraction = 0.1;

struct OptimizerSettings {
  Optimizer method = Optimizer::kNelderMead;
  // Nelder-Mead budget: iterations per squared parameter count, the
  // forecast package's 100 dim^2 evaluations.
  size_t iterations_per_dim_squared = 100;
  // relative tolerance on the objective, R's optim reltol
  double tolerance = 1e-8;
  // iteration budget of the quasi-Newton methods
  size_t max_iterations = 200;
  // Nelder-Mead restarts from the point found, with a fresh simplex. The
  // simplex collapses in narrow valleys before reaching the floor; on the
  // frozen fits two restarts recover the reference optimum within 0.01
  // percent where a single run was up to 1 percent short, at twice the
  // evaluations.
  size_t restarts = 2;
  // initial simplex step as a fraction of the largest scaled coordinate
  double simplex_fraction = kSimplexFraction;
  // budget for the evaluation-driven methods
  size_t max_evaluations = 20000;
};

struct OptimizationResult {
  double value;
  size_t evaluations;
  bool converged;
};

template <typename Objective>
OptimizationResult minimise(const OptimizerSettings &settings, Objective &f,
                            std::vector<double> &x) {
  const size_t dim = x.size();
  switch (settings.method) {
    case Optimizer::kNelderMead: {
      double largest = 0;
      for (const double v : x) largest = std::max(largest, std::abs(v));
      const double step = largest > 0 ? settings.simplex_fraction * largest
                                      : settings.simplex_fraction;
      const size_t max_iter = settings.iterations_per_dim_squared * dim * dim;
      nlsolver::NelderMead<Objective, double> solver(
          f, step, 1.0, 2.0, 0.5, 0.5, settings.tolerance, max_iter,
          max_iter, settings.restarts);
      const auto status = solver.minimize(x);
      const auto [calls, iters, value, grads, hess] = status.get_summary();
      return {value, calls, iters < max_iter};
    }
    case Optimizer::kBFGS: {
      nlsolver::BFGS<Objective, double> solver(f, {}, settings.max_iterations);
      const auto status = solver.minimize(x);
      const auto [calls, iters, value, grads, hess] = status.get_summary();
      return {value, calls, status.success()};
    }
    case Optimizer::kLBFGSB: {
      nlsolver::LBFGSB<Objective, double> solver(f, {}, settings.max_iterations);
      const auto status = solver.minimize(x);
      const auto [calls, iters, value, grads, hess] = status.get_summary();
      return {value, calls, status.success()};
    }
    case Optimizer::kPatternSearch: {
      nlsolver::PatternSearch<Objective, double> solver(
          f, 1.0, settings.max_evaluations);
      const auto status = solver.minimize(x);
      const auto [calls, iters, value, grads, hess] = status.get_summary();
      return {value, calls, status.success()};
    }
  }
  throw std::logic_error("tbats: unhandled optimizer");
}

// ===========================================================================
// 7. ARMA by exact likelihood, and order selection
// ===========================================================================
//
// Self-contained: an ARMA(p, q) with a mean, fitted by exact Gaussian maximum
// likelihood through the Kalman filter, exactly as R's arima() does for
// method = "ML" without differencing or regressors. Used to choose the ARMA
// orders of the error process from the residuals of a fit without ARMA
// terms; the coefficients are re-estimated inside the TBATS likelihood.
//
// [[CITATION]] Gardner, Harvey and Phillips (1980), "An algorithm for exact
// maximum likelihood estimation of autoregressive-moving average models by
// means of Kalman filtering", Applied Statistics 29(3), 311-322 (the state
// space form and the initial state covariance, algorithm AS 154).
// [[CITATION]] Jones (1980), "Maximum likelihood fitting of ARMA models to
// time series with missing observations", Technometrics 22(3), 389-395 (the
// partial autocorrelation parametrisation of a stationary AR part).

namespace arma {

// AR coefficients from unconstrained values: tanh maps each to a partial
// autocorrelation in (-1, 1), and the Durbin-Levinson recursion turns the
// partial autocorrelations into the coefficients of a stationary AR(p).
inline void ar_from_unconstrained(const double *raw, const size_t p,
                                  double *ar, std::vector<double> &work) {
  work.assign(p, 0.0);
  for (size_t j = 0; j < p; ++j) ar[j] = work[j] = std::tanh(raw[j]);
  for (size_t j = 0; j < p; ++j) {
    const double a = ar[j];
    for (size_t k = 0; k < j; ++k) work[k] -= a * ar[j - k - 1];
    for (size_t k = 0; k < j; ++k) ar[k] = work[k];
  }
}

// Kalman filter for a zero-mean ARMA(p, q) in the Gardner, Harvey and
// Phillips state space form of dimension r = max(p, q + 1): the first state
// is the series, the transition carries phi in its first column and a shift
// above the diagonal, and the innovation enters through (1, theta_1, ...,
// theta_{r-1}). Everything is in units of the innovation variance, which is
// concentrated out of the likelihood. All buffers are sized once.
class KalmanArma {
 public:
  KalmanArma(const size_t p, const size_t q)
      : p_(p),
        q_(q),
        r_(std::max(p, q + 1)),
        np_(r_ * (r_ + 1) / 2),
        nrbar_(np_ * (np_ - 1) / 2),
        phi_(p),
        theta_(q),
        a_(r_),
        P_(r_ * r_),
        Pn_(r_ * r_),
        anew_(r_),
        M_(r_),
        xnext_(np_),
        xrow_(np_),
        rbar_(nrbar_),
        thetab_(np_),
        V_(np_),
        Q_(r_ * r_) {}

  size_t p() const { return p_; }
  size_t q() const { return q_; }

  // Sets the coefficients and puts the filter at its stationary start: zero
  // state and the stationary covariance of the state as initial uncertainty.
  void reset(const double *ar, const double *ma) {
    std::copy(ar, ar + p_, phi_.begin());
    std::copy(ma, ma + q_, theta_.begin());
    std::fill(a_.begin(), a_.end(), 0.0);
    stationary_covariance(Pn_);
    for (size_t i = 0; i < r_; ++i) M_[i] = Pn_[i];
  }

  struct Result {
    double sigma2;     // concentrated innovation variance
    double objective;  // 0.5 (log sigma2 + mean log gain): -log L / n + const
    size_t used;       // observations entering the likelihood
  };

  // Runs the filter over the centred series. The one-step prediction error
  // e_t and its variance gain_t give -2 log L = sum log gain_t + n log(sum
  // e_t^2 / gain_t / n) up to constants.
  //
  // The covariance recursion converges to a fixed point for a stationary
  // model, after which gain and Kalman vector are constant; once two
  // successive prediction covariances agree to kSteadyStateTol they are
  // frozen and each step costs O(r) instead of O(r^2). The tolerance is far
  // below the 1e-8 at which the likelihood is compared with R's.
  Result likelihood(const double *y, const size_t n) {
    double ssq = 0, sumlog = 0;
    size_t used = n;
    bool steady = false;
    double gain = Pn_[0];
    for (size_t l = 0; l < n; ++l) {
      // state prediction: anew = T a
      for (size_t i = 0; i < r_; ++i) {
        double tmp = (i + 1 < r_) ? a_[i + 1] : 0.0;
        if (i < p_) tmp += phi_[i] * a_[0];
        anew_[i] = tmp;
      }
      // covariance prediction: Pn = T P T' + R R', skipped at l = 0 where Pn
      // holds the stationary covariance
      if (l > 0 && !steady) {
        double change = 0, scale = 0;
        for (size_t i = 0; i < r_; ++i) {
          const double vi = innovation_weight(i);
          for (size_t j = 0; j < r_; ++j) {
            double tmp = vi * innovation_weight(j);
            if (i < p_ && j < p_) tmp += phi_[i] * phi_[j] * P_[0];
            if (i + 1 < r_ && j + 1 < r_) tmp += P_[i + 1 + r_ * (j + 1)];
            if (i < p_ && j + 1 < r_) tmp += phi_[i] * P_[j + 1];
            if (j < p_ && i + 1 < r_) tmp += phi_[j] * P_[i + 1];
            change = std::max(change, std::abs(tmp - Pn_[i + r_ * j]));
            scale = std::max(scale, std::abs(tmp));
            Pn_[i + r_ * j] = tmp;
          }
        }
        steady = change <= kSteadyStateTol * (1.0 + scale);
        for (size_t i = 0; i < r_; ++i) M_[i] = Pn_[i];
        gain = M_[0];
      }
      // measurement update with the first state observed
      const double resid = y[l] - anew_[0];
      if (gain < kGainCutoff) {
        ssq += resid * resid / gain;
        sumlog += std::log(gain);
      } else {
        --used;
      }
      for (size_t i = 0; i < r_; ++i) a_[i] = anew_[i] + M_[i] * resid / gain;
      if (!steady) {
        for (size_t i = 0; i < r_; ++i) {
          for (size_t j = 0; j < r_; ++j) {
            P_[i + j * r_] = Pn_[i + j * r_] - M_[i] * M_[j] / gain;
          }
        }
      }
    }
    const double nu = static_cast<double>(used);
    const double sigma2 = ssq / nu;
    return {sigma2, 0.5 * (std::log(sigma2) + sumlog / nu), used};
  }

 private:
  // A prediction variance this large marks a diffuse initial state; R
  // excludes such observations from the likelihood, and so does this filter.
  static constexpr double kGainCutoff = 1e4;
  // relative change of the prediction covariance below which it is frozen
  static constexpr double kSteadyStateTol = 1e-13;

  double innovation_weight(const size_t i) const {
    if (i == 0) return 1.0;
    return (i - 1 < q_) ? theta_[i - 1] : 0.0;
  }

  // Givens-style inclusion of one row into the triangular system of
  // algorithm AS 154 (R's inclu2).
  void include_row(double ynext) {
    for (size_t i = 0; i < np_; ++i) xrow_[i] = xnext_[i];
    size_t ithisr = 0;
    for (size_t i = 0; i < np_; ++i) {
      if (xrow_[i] != 0.0) {
        const double xi = xrow_[i];
        const double di = Q_[i];
        const double dpi = di + xi * xi;
        Q_[i] = dpi;
        const double cbar = di / dpi;
        const double sbar = xi / dpi;
        for (size_t k = i + 1; k < np_; ++k) {
          const double xk = xrow_[k];
          const double rbthis = rbar_[ithisr];
          xrow_[k] = xk - xi * rbthis;
          rbar_[ithisr++] = cbar * rbthis + sbar * xk;
        }
        const double xk = ynext;
        ynext = xk - xi * thetab_[i];
        thetab_[i] = cbar * thetab_[i] + sbar * xk;
        if (di == 0.0) return;
      } else {
        ithisr += np_ - i - 1;
      }
    }
  }

  // Stationary covariance Q0 of the state (R's getQ0, algorithm AS 154): the
  // solution of Q0 = T Q0 T' + R R', formed as a linear system in the packed
  // upper triangle of Q0 and solved by orthogonal triangularisation. For a
  // pure moving average the solution is a back substitution.
  void stationary_covariance(std::vector<double> &out) {
    std::vector<double> &P = Q_;  // packed, then unpacked in place
    size_t ind = 0;
    for (size_t j = 0; j < r_; ++j) {
      const double vj = innovation_weight(j);
      for (size_t i = j; i < r_; ++i) V_[ind++] = innovation_weight(i) * vj;
    }
    if (r_ == 1) {
      out[0] = (p_ == 0) ? 1.0 : 1.0 / (1.0 - phi_[0] * phi_[0]);
      return;
    }
    if (p_ > 0) {
      std::fill(rbar_.begin(), rbar_.end(), 0.0);
      std::fill(P.begin(), P.end(), 0.0);
      std::fill(thetab_.begin(), thetab_.end(), 0.0);
      std::fill(xnext_.begin(), xnext_.end(), 0.0);
      const size_t npr = np_ - r_, npr1 = npr + 1;
      ind = 0;
      size_t ind1_next = 0;  // next free slot of the -1 entries
      size_t indj = npr;
      size_t ind2 = npr - 1;
      for (size_t j = 0; j < r_; ++j) {
        const double phij = (j < p_) ? phi_[j] : 0.0;
        xnext_[indj++] = 0.0;
        size_t indi = npr1 + j;
        for (size_t i = j; i < r_; ++i) {
          const double ynext = V_[ind++];
          const double phii = (i < p_) ? phi_[i] : 0.0;
          if (j + 1 != r_) {
            xnext_[indj] = -phii;
            if (i + 1 != r_) {
              xnext_[indi] -= phij;
              xnext_[ind1_next++] = -1.0;
            }
          }
          xnext_[npr] = -phii * phij;
          if (++ind2 >= np_) ind2 = 0;
          xnext_[ind2] += 1.0;
          include_row(ynext);
          xnext_[ind2] = 0.0;
          if (i + 1 != r_) {
            xnext_[indi++] = 0.0;
            xnext_[ind1_next - 1] = 0.0;
          }
        }
      }
      // back substitution through the triangular factor
      size_t ithisr = nrbar_ - 1;
      size_t im = np_ - 1;
      for (size_t i = 0; i < np_; ++i) {
        double bi = thetab_[im];
        size_t jm = np_ - 1;
        for (size_t j = 0; j < i; ++j) bi -= rbar_[ithisr--] * P[jm--];
        P[im--] = bi;
      }
      // move the r diagonal-block entries to the front
      ind = npr;
      for (size_t i = 0; i < r_; ++i) xnext_[i] = P[ind++];
      ind = np_ - 1;
      size_t ind1 = npr - 1;
      for (size_t i = 0; i < npr; ++i) P[ind--] = P[ind1--];
      for (size_t i = 0; i < r_; ++i) P[i] = xnext_[i];
    } else {
      size_t indn = np_;
      ind = np_;
      for (size_t i = 0; i < r_; ++i) {
        for (size_t j = 0; j <= i; ++j) {
          --ind;
          P[ind] = V_[ind];
          if (j != 0) P[ind] += P[--indn];
        }
      }
    }
    // unpack the triangle into the full r by r matrix
    ind = np_;
    for (size_t i = r_ - 1; i > 0; --i) {
      for (size_t j = r_ - 1; j >= i; --j) P[r_ * i + j] = P[--ind];
    }
    for (size_t i = 0; i + 1 < r_; ++i) {
      for (size_t j = i + 1; j < r_; ++j) P[i + r_ * j] = P[j + r_ * i];
    }
    std::copy(P.begin(), P.begin() + r_ * r_, out.begin());
  }

  size_t p_, q_, r_, np_, nrbar_;
  std::vector<double> phi_, theta_, a_, P_, Pn_, anew_, M_;
  std::vector<double> xnext_, xrow_, rbar_, thetab_, V_, Q_;
};

// Objective of the exact likelihood over (unconstrained AR values, MA
// coefficients, mean), as R's arima() optimises it.
class ArmaObjective {
 public:
  ArmaObjective(const std::vector<double> &y, const size_t p, const size_t q)
      : y_(y), filter_(p, q), centred_(y.size()), ar_(p), ma_(q), work_(p) {}

  size_t count() const { return filter_.p() + filter_.q() + 1; }

  double operator()(std::vector<double> &x) {
    const size_t p = filter_.p(), q = filter_.q();
    ar_from_unconstrained(x.data(), p, ar_.data(), work_);
    std::copy(x.begin() + p, x.begin() + p + q, ma_.begin());
    const double mean = x[p + q];
    for (size_t t = 0; t < y_.size(); ++t) centred_[t] = y_[t] - mean;
    filter_.reset(ar_.data(), ma_.data());
    const KalmanArma::Result result = filter_.likelihood(centred_.data(), y_.size());
    last_ = result;
    return std::isfinite(result.objective) ? result.objective
                                           : kInadmissibleValue;
  }
  const std::vector<double> &ar() const { return ar_; }
  const std::vector<double> &ma() const { return ma_; }
  const KalmanArma::Result &last() const { return last_; }

 private:
  std::vector<double> y_;
  KalmanArma filter_;
  std::vector<double> centred_, ar_, ma_, work_;
  KalmanArma::Result last_{0, 0, 0};
};

struct ArmaFit {
  size_t p, q;
  std::vector<double> ar, ma;
  double mean, sigma2, loglik, aic;
  size_t evaluations;
  bool converged;
};

// log(2 pi), for the Gaussian log-likelihood constant
constexpr double kLogTwoPi = 1.8378770664093453;

inline Optimizer other_quasi_newton(const Optimizer method) {
  return method == Optimizer::kBFGS ? Optimizer::kLBFGSB : Optimizer::kBFGS;
}

// Exact maximum likelihood fit of ARMA(p, q) with a mean, starting from zero
// coefficients and the sample mean. The AIC counts p + q + 1 coefficients
// plus the innovation variance, as R's arima() does. With `second_opinion`
// the other quasi-Newton method also runs from the start and the better
// fit is kept: each method stops at a local optimum short of the other's
// on some orders, and a polish from the stalled point does not escape it.
inline ArmaFit fit_arma(const std::vector<double> &y, const size_t p,
                        const size_t q, const OptimizerSettings &settings,
                        const bool second_opinion = true) {
  const double n = static_cast<double>(y.size());
  ArmaObjective objective(y, p, q);
  std::vector<double> x(p + q + 1, 0.0);
  double mean = 0;
  for (const double v : y) mean += v;
  x.back() = mean / n;
  const std::vector<double> start = x;
  OptimizationResult opt = minimise(settings, objective, x);
  if (second_opinion) {
    OptimizerSettings other = settings;
    other.method = other_quasi_newton(settings.method);
    std::vector<double> x_other = start;
    const OptimizationResult second = minimise(other, objective, x_other);
    opt.evaluations += second.evaluations;
    if (second.value < opt.value) {
      x = x_other;
      opt.value = second.value;
      opt.converged = second.converged;
    }
  }
  const double value = objective(x);  // leaves the coefficients in place
  const double loglik = -(n * value + 0.5 * n * (1.0 + kLogTwoPi));
  const double aic = -2.0 * loglik + 2.0 * static_cast<double>(p + q + 1) + 2.0;
  return {p, q, objective.ar(), objective.ma(), x.back(),
          objective.last().sigma2, loglik, aic, opt.evaluations, opt.converged};
}

struct ArmaOrder {
  size_t p, q;
  double aic;
};

// Optimiser for the ARMA fits. The likelihood is smooth in the transformed
// parameters, so a quasi-Newton method converges in a few hundred
// evaluations. The iteration budget is R's arima() default.
inline OptimizerSettings arma_optimizer_defaults() {
  OptimizerSettings settings;
  settings.method = Optimizer::kLBFGSB;
  settings.max_iterations = 100;
  return settings;
}

// AIC distance from the best order within which an order is refitted with
// the other quasi-Newton method. Each method stopped at a local optimum
// short of R's on one order of the frozen residual series, by 0.3 and 0.7
// AIC units, and R's own optimiser by 1.2 on another; a stall farther from
// the top than this window cannot change which orders are proposed.
constexpr double kSecondOpinionWindow = 3.0;

// Largest orders the selection considers, as in the forecast package.
constexpr size_t kMaxArmaOrder = 5;

// AIC of every order with p <= max_p, q <= max_q, in p-major order; an
// order whose fit fails is left out. The grid is fitted with one method;
// the orders within kSecondOpinionWindow of the best are then refitted with
// the other method and keep the better value.
inline std::vector<ArmaOrder> arma_order_grid(const std::vector<double> &y,
                                              const size_t max_p,
                                              const size_t max_q,
                                              const OptimizerSettings &settings) {
  std::vector<ArmaOrder> grid;
  double best = std::numeric_limits<double>::infinity();
  for (size_t p = 0; p <= max_p; ++p) {
    for (size_t q = 0; q <= max_q; ++q) {
      double aic;
      try {
        aic = fit_arma(y, p, q, settings, false).aic;
      } catch (const std::exception &) {
        continue;
      }
      if (!std::isfinite(aic)) continue;
      grid.push_back({p, q, aic});
      best = std::min(best, aic);
    }
  }
  OptimizerSettings other = settings;
  other.method = other_quasi_newton(settings.method);
  for (ArmaOrder &order : grid) {
    if (order.aic > best + kSecondOpinionWindow) continue;
    try {
      order.aic = std::min(order.aic, fit_arma(y, order.p, order.q, other, false).aic);
    } catch (const std::exception &) {
    }
  }
  return grid;
}

// The order with the smallest AIC, ties resolved by the first in p-major
// order.
inline ArmaOrder best_order(const std::vector<ArmaOrder> &grid) {
  ArmaOrder best{0, 0, std::numeric_limits<double>::infinity()};
  for (const auto &o : grid) {
    if (o.aic < best.aic) best = o;
  }
  return best;
}

// Stepwise order search: four starting orders, then the six neighbours
// (p +- 1, q +- 1 singly, and both together) of the best fitted so far,
// repeated until no neighbour improves, every order fitted once. Typically
// ten to fifteen fits instead of the full grid's thirty-six; the orders it
// visits are the ones returned, so the second opinion and the parsimonious
// candidate draw on them.
// [[CITATION]] Hyndman and Khandakar (2008), "Automatic time series
// forecasting: the forecast package for R", Journal of Statistical Software
// 27(3), section 3.2.
inline std::vector<ArmaOrder> arma_order_stepwise(const std::vector<double> &y,
                                                  const size_t max_p,
                                                  const size_t max_q,
                                                  const OptimizerSettings &settings) {
  std::vector<ArmaOrder> visited;
  auto aic_of = [&](const size_t p, const size_t q) -> double {
    for (const auto &o : visited) {
      if (o.p == p && o.q == q) return o.aic;
    }
    double aic = std::numeric_limits<double>::infinity();
    try {
      aic = fit_arma(y, p, q, settings, false).aic;
    } catch (const std::exception &) {
    }
    if (!std::isfinite(aic)) aic = std::numeric_limits<double>::infinity();
    visited.push_back({p, q, aic});
    return aic;
  };
  ArmaOrder best{0, 0, std::numeric_limits<double>::infinity()};
  auto consider = [&](const size_t p, const size_t q) {
    if (p > max_p || q > max_q) return;
    const double aic = aic_of(p, q);
    if (aic < best.aic) best = {p, q, aic};
  };
  consider(2, 2);
  consider(0, 0);
  consider(1, 0);
  consider(0, 1);
  while (true) {
    const ArmaOrder current = best;
    const size_t p = current.p, q = current.q;
    if (p > 0) consider(p - 1, q);
    consider(p + 1, q);
    if (q > 0) consider(p, q - 1);
    consider(p, q + 1);
    if (p > 0 && q > 0) consider(p - 1, q - 1);
    consider(p + 1, q + 1);
    if (best.p == current.p && best.q == current.q) break;
  }
  OptimizerSettings other = settings;
  other.method = other_quasi_newton(settings.method);
  for (ArmaOrder &order : visited) {
    if (!std::isfinite(order.aic) || order.aic > best.aic + kSecondOpinionWindow) {
      continue;
    }
    try {
      order.aic = std::min(order.aic, fit_arma(y, order.p, order.q, other, false).aic);
    } catch (const std::exception &) {
    }
  }
  std::vector<ArmaOrder> out;
  for (const auto &o : visited) {
    if (std::isfinite(o.aic)) out.push_back(o);
  }
  return out;
}

enum class ArmaSearch { kStepwise, kGrid };

inline std::vector<ArmaOrder> arma_orders(const std::vector<double> &y,
                                          const size_t max_p, const size_t max_q,
                                          const OptimizerSettings &settings,
                                          const ArmaSearch search) {
  return search == ArmaSearch::kGrid
             ? arma_order_grid(y, max_p, max_q, settings)
             : arma_order_stepwise(y, max_p, max_q, settings);
}

inline ArmaOrder select_arma_order(const std::vector<double> &y,
                                   const size_t max_p, const size_t max_q,
                                   const OptimizerSettings &settings,
                                   const ArmaSearch search = ArmaSearch::kGrid) {
  return best_order(arma_orders(y, max_p, max_q, settings, search));
}

// AIC difference below which two orders are not told apart. The residual
// fit only proposes orders; the coefficients are re-estimated inside the
// model's likelihood, where each costs two AIC units, so among near-ties the
// most parsimonious order is a candidate in its own right.
constexpr double kParsimonyWindow = 2.0;

// Among the orders within kParsimonyWindow of the best, the one with the
// fewest coefficients (lowest AIC among equals).
inline ArmaOrder most_parsimonious_near(const std::vector<ArmaOrder> &grid,
                                        const ArmaOrder &best) {
  ArmaOrder out = best;
  for (const auto &o : grid) {
    if (o.aic > best.aic + kParsimonyWindow) continue;
    const size_t size = o.p + o.q, out_size = out.p + out.q;
    if (size < out_size || (size == out_size && o.aic < out.aic)) out = o;
  }
  return out;
}

}  // namespace arma

// ---- fitting one specification --------------------------------------------

// Starting values of the forecast package. Dummy seasonal models with long
// periods start the level and trend smoothing near zero because their many
// states make large smoothing steps unstable.
constexpr double kInitialAlpha = 0.09;
constexpr double kInitialBeta = 0.05;
constexpr double kInitialAlphaLongPeriod = 1e-6;
constexpr double kInitialBetaLongPeriod = 5e-7;
constexpr double kLongPeriodTotal = 16;
constexpr double kInitialPhi = 0.999;
constexpr double kInitialGammaDummy = 0.001;
constexpr double kInitialGammaTrigonometric = 0.0;

inline Parameters initial_parameters(const ModelSpec &spec,
                                     const double init_lambda) {
  Parameters par;
  const size_t m = spec.seasonal.size();
  double total_period = 0;
  for (const auto &s : spec.seasonal) total_period += s.period;
  const bool dummy = spec.seasonal_type == SeasonalType::kDummy;
  const bool long_period = dummy && total_period > kLongPeriodTotal;
  par.alpha = long_period ? kInitialAlphaLongPeriod : kInitialAlpha;
  if (spec.trend) par.beta = long_period ? kInitialBetaLongPeriod : kInitialBeta;
  par.phi = spec.damping ? kInitialPhi : 1.0;
  par.gamma_one.assign(m, dummy ? kInitialGammaDummy : kInitialGammaTrigonometric);
  if (!dummy) par.gamma_two.assign(m, kInitialGammaTrigonometric);
  par.ar.assign(spec.p, 0.0);
  par.ma.assign(spec.q, 0.0);
  if (spec.box_cox) par.lambda = init_lambda;
  return par;
}

struct FitSettings {
  OptimizerSettings optimizer;
  bool bias_adjust = false;  // bias-adjusted fitted values after Box-Cox
  ForecastabilityCheck forecastability = ForecastabilityCheck::kAuto;
};

// One specification fitted to one series.
struct FittedSpec {
  FittedSpec(const ModelSpec &spec_, const Parameters &par,
             const std::vector<double> &seed, FilterResult &&result,
             std::vector<double> &&fitted_, const double variance_,
             const double neg2loglik_, const OptimizationResult &opt)
      : spec(spec_),
        parameters(par),
        seed_states(seed),
        filter(std::move(result)),
        fitted(std::move(fitted_)),
        variance(variance_),
        neg2loglik(neg2loglik_),
        aic(neg2loglik_ +
            2.0 * static_cast<double>(parameter_count(spec_) + seed.size())),
        evaluations(opt.evaluations),
        converged(opt.converged) {}
  ModelSpec spec;
  Parameters parameters;
  std::vector<double> seed_states;  // model scale
  FilterResult filter;              // errors, predictions, states; model scale
  std::vector<double> fitted;       // original scale
  double variance;                  // sum e_t^2 / n on the model scale
  double neg2loglik;
  double aic;  // neg2loglik + 2 (parameters + seed states), the paper's count
  size_t evaluations;
  bool converged;
};

// Seed states for `par`, estimated on the model scale of par.lambda and
// returned on the original scale when the specification has a Box-Cox
// transformation, so that they can be re-transformed with other lambdas.
inline std::vector<double> seed_states_for(const ModelSpec &spec,
                                           const Parameters &par,
                                           const std::vector<double> &y) {
  const System system(spec, par);
  std::vector<double> x0 = seed::seed_states(
      system, spec.box_cox ? box_cox::transform(y, par.lambda) : y);
  if (spec.box_cox) {
    for (double &v : x0) v = box_cox::inverse(v, par.lambda);
  }
  return x0;
}

inline FittedSpec fit_specific(const std::vector<double> &y,
                               const ModelSpec &spec, const double init_lambda,
                               const FitSettings &settings) {
  spec.validate();
  const Parameters start = initial_parameters(spec, init_lambda);
  const std::vector<double> seed = seed_states_for(spec, start, y);
  Likelihood objective(spec, y, seed, settings.forecastability);
  std::vector<double> x = objective.scale(start);
  const OptimizationResult opt = minimise(settings.optimizer, objective, x);
  const Parameters par = objective.unscale(x);
  const double neg2loglik = objective.evaluate(par);  // full check
  const std::vector<double> x0 = objective.seed_for(par);
  const System system(spec, par);
  FilterResult result =
      filter(system, spec.box_cox ? box_cox::transform(y, par.lambda) : y, x0);
  double sse = 0;
  for (const double e : result.errors) sse += e * e;
  const double variance = sse / static_cast<double>(y.size());
  std::vector<double> fitted = result.predictions;
  if (spec.box_cox) {
    for (double &v : fitted) {
      v = settings.bias_adjust
              ? box_cox::inverse_bias_adjusted(v, par.lambda, variance)
              : box_cox::inverse(v, par.lambda);
    }
  }
  return FittedSpec(spec, par, x0, std::move(result), std::move(fitted),
                    variance, neg2loglik, opt);
}

// ===========================================================================
// 8. Model search
// ===========================================================================

// What the search may vary. An axis with no value tries both settings; a
// pinned axis tries one. Damping without a trend is never tried.
struct SearchOptions {
  std::optional<bool> box_cox;  // disabled regardless when the series is not positive
  std::optional<bool> trend;
  std::optional<bool> damping;
  bool arma_errors = true;
  double box_cox_lower = 0.0;
  double box_cox_upper = 1.0;
  bool bias_adjust = false;
  OptimizerSettings optimizer;
  OptimizerSettings arma_optimizer = arma::arma_optimizer_defaults();
  size_t max_arma_order = arma::kMaxArmaOrder;
  arma::ArmaSearch arma_search = arma::ArmaSearch::kStepwise;
  ForecastabilityCheck forecastability = ForecastabilityCheck::kAuto;
};

namespace search {

inline std::vector<bool> axis(const std::optional<bool> &pinned) {
  if (pinned) return {*pinned};
  return {false, true};
}

inline bool any_true(const std::vector<bool> &values) {
  for (const bool v : values) {
    if (v) return true;
  }
  return false;
}

// Relative tolerance under which a series counts as constant (R's all.equal).
constexpr double kConstantTolerance = 1.5e-8;

inline bool is_constant(const std::vector<double> &y) {
  double scale = 0;
  for (const double v : y) scale += std::abs(v);
  scale = scale / static_cast<double>(y.size());
  const double tolerance = kConstantTolerance * std::max(scale, 1.0);
  for (const double v : y) {
    if (std::abs(v - y[0]) > tolerance) return false;
  }
  return true;
}

// The model of a constant series: a level equal to the series, no
// uncertainty. A member of the family with one state and nothing to
// estimate, so every consumer handles it like any other fit.
inline FittedSpec constant_series_fit(const std::vector<double> &y) {
  ModelSpec spec;
  spec.seasonal_type = SeasonalType::kDummy;
  Parameters par;
  par.alpha = 0.9999;  // the forecast package's value; the level never moves
  const std::vector<double> seed(1, y[0]);
  const System system(spec, par);
  FilterResult result = filter(system, y, seed);
  std::vector<double> fitted = result.predictions;
  const OptimizationResult opt{-std::numeric_limits<double>::infinity(), 0, true};
  return FittedSpec(spec, par, seed, std::move(result), std::move(fitted), 0.0,
                    -std::numeric_limits<double>::infinity(), opt);
}

inline double aic_or_infinity(const std::optional<FittedSpec> &fit) {
  return fit ? fit->aic : std::numeric_limits<double>::infinity();
}

// A fit that fails (singular seed regression, no admissible parameters) is
// a candidate with infinite AIC, as in the forecast package.
inline std::optional<FittedSpec> try_fit(const std::vector<double> &y,
                                         const ModelSpec &spec,
                                         const double init_lambda,
                                         const SearchOptions &options) {
  FitSettings settings;
  settings.optimizer = options.optimizer;
  settings.bias_adjust = options.bias_adjust;
  settings.forecastability = options.forecastability;
  try {
    FittedSpec fit = fit_specific(y, spec, init_lambda, settings);
    if (!std::isfinite(fit.aic)) return std::nullopt;
    return fit;
  } catch (const std::exception &) {
    return std::nullopt;
  }
}

inline void keep_better(std::optional<FittedSpec> &best,
                        std::optional<FittedSpec> &&candidate) {
  if (candidate && candidate->aic < aic_or_infinity(best)) {
    best = std::move(candidate);
  }
}

// ARMA errors: rank the orders on the errors of `fit`, refit with the best
// order and with the most parsimonious order near it, and keep whichever
// of the three fits has the lowest AIC. On short series the top orders are
// near-ties by the residual fit and the parsimonious one is often the
// better model once its coefficients are re-estimated.
inline std::optional<FittedSpec> with_arma_errors(
    const std::vector<double> &y, std::optional<FittedSpec> fit,
    const double init_lambda, const SearchOptions &options) {
  if (!fit || !options.arma_errors) return fit;
  std::vector<arma::ArmaOrder> grid;
  try {
    grid = arma::arma_orders(fit->filter.errors, options.max_arma_order,
                             options.max_arma_order, options.arma_optimizer,
                             options.arma_search);
  } catch (const std::exception &) {
    return fit;
  }
  if (grid.empty()) return fit;
  const arma::ArmaOrder best = arma::best_order(grid);
  const arma::ArmaOrder parsimonious = arma::most_parsimonious_near(grid, best);
  std::vector<arma::ArmaOrder> candidates{best};
  if (parsimonious.p != best.p || parsimonious.q != best.q) {
    candidates.push_back(parsimonious);
  }
  std::optional<FittedSpec> out = std::move(fit);
  for (const arma::ArmaOrder &order : candidates) {
    if (order.p == 0 && order.q == 0) continue;
    ModelSpec with_arma = out->spec;
    with_arma.p = order.p;
    with_arma.q = order.q;
    keep_better(out, try_fit(y, with_arma, init_lambda, options));
  }
  return out;
}

inline ModelSpec cell_spec(const SeasonalType type,
                           const std::vector<SeasonalPeriod> &seasonal,
                           const bool box_cox, const bool trend,
                           const bool damping, const SearchOptions &options) {
  ModelSpec spec;
  spec.box_cox = box_cox;
  spec.box_cox_lower = options.box_cox_lower;
  spec.box_cox_upper = options.box_cox_upper;
  spec.trend = trend;
  spec.damping = damping;
  spec.seasonal_type = type;
  spec.seasonal = seasonal;
  return spec;
}

// The Box-Cox parameter the search starts from: Guerrero's choice with
// blocks of the largest seasonal period (2 when there is none).
inline double starting_lambda(const std::vector<double> &y,
                              const std::vector<SeasonalPeriod> &seasonal,
                              const SearchOptions &options) {
  double largest = 0;
  for (const auto &s : seasonal) largest = std::max(largest, s.period);
  const size_t period = seasonal.empty() ? 2 : static_cast<size_t>(largest);
  return box_cox::guerrero_lambda(y, period, options.box_cox_lower,
                                  options.box_cox_upper);
}

inline bool positive(const std::vector<double> &y) {
  for (const double v : y) {
    if (v <= 0) return false;
  }
  return true;
}

}  // namespace search

// BATS: dummy seasonality with the given integer periods (none for a
// non-seasonal model). Every cell of the Box-Cox, trend and damping grid is
// fitted; with seasonal periods the cell's non-seasonal counterpart competes
// with it; the ARMA refinement is applied to the winner; the lowest AIC over
// the grid is the model.
inline FittedSpec search_bats(const std::vector<double> &y,
                              const std::vector<double> &periods,
                              SearchOptions options) {
  if (y.empty()) throw std::invalid_argument("tbats: empty series");
  std::vector<SeasonalPeriod> seasonal;
  for (const double m : periods) {
    if (m > 1.0) seasonal.push_back({m, 0});
  }
  if (search::is_constant(y)) return search::constant_series_fit(y);
  if (!search::positive(y)) options.box_cox = false;
  if (options.trend && !*options.trend) options.damping = false;
  const std::vector<bool> box_cox_axis = search::axis(options.box_cox);
  const std::vector<bool> trend_axis = search::axis(options.trend);
  const std::vector<bool> damping_axis = search::axis(options.damping);
  const double init_lambda =
      search::any_true(box_cox_axis)
          ? search::starting_lambda(y, seasonal, options)
          : 1.0;
  std::optional<FittedSpec> best;
  for (const bool box_cox : box_cox_axis) {
    for (const bool trend : trend_axis) {
      for (const bool damping : damping_axis) {
        if (damping && !trend) continue;
        std::optional<FittedSpec> cell = search::try_fit(
            y,
            search::cell_spec(SeasonalType::kDummy, seasonal, box_cox, trend,
                              damping, options),
            init_lambda, options);
        if (!seasonal.empty()) {
          std::optional<FittedSpec> plain = search::try_fit(
              y,
              search::cell_spec(SeasonalType::kDummy, {}, box_cox, trend,
                                damping, options),
              init_lambda, options);
          if (plain && plain->aic < search::aic_or_infinity(cell)) {
            cell = std::move(plain);
          }
        }
        search::keep_better(
            best, search::with_arma_errors(y, std::move(cell), init_lambda,
                                           options));
      }
    }
  }
  if (!best) throw std::runtime_error("tbats: unable to fit a BATS model");
  return *best;
}

namespace search {

// Harmonics per period chosen by AIC, in the most general cell of the grid.
// The maximum is floor((m - 1) / 2), lowered where a lower harmonic order
// would repeat the frequencies of an earlier period. Up to six harmonics the
// search starts at the maximum and steps down while the AIC improves; above
// six it compares five, six and seven and walks in the improving direction.
// Returns the fit with the chosen harmonics, or nothing if no fit succeeded.
inline std::optional<FittedSpec> choose_harmonics(
    const std::vector<double> &y, std::vector<SeasonalPeriod> &seasonal,
    const bool box_cox, const bool trend, const bool damping,
    const double init_lambda, const SearchOptions &options) {
  auto fit_with = [&](const std::vector<SeasonalPeriod> &s) {
    return try_fit(y,
                   cell_spec(SeasonalType::kTrigonometric, s, box_cox, trend,
                             damping, options),
                   init_lambda, options);
  };
  constexpr size_t kStepSearchLimit = 6;
  std::optional<FittedSpec> best = fit_with(seasonal);
  for (size_t i = 0; i < seasonal.size(); ++i) {
    const double m = seasonal[i].period;
    if (m == 2.0) continue;
    size_t max_k = static_cast<size_t>(std::floor((m - 1.0) / 2.0));
    if (i != 0) {
      for (size_t k = 2; k <= max_k; ++k) {
        if (std::fmod(m, static_cast<double>(k)) != 0.0) continue;
        const double latter = m / static_cast<double>(k);
        bool repeats = false;
        for (size_t j = 0; j < i; ++j) {
          if (std::fmod(seasonal[j].period, latter) == 0.0) repeats = true;
        }
        if (repeats) {
          max_k = k - 1;
          break;
        }
      }
    }
    if (max_k == 1) continue;
    size_t &k = seasonal[i].harmonics;
    if (max_k <= kStepSearchLimit) {
      k = max_k;
      double best_aic = std::numeric_limits<double>::infinity();
      while (true) {
        std::optional<FittedSpec> candidate = fit_with(seasonal);
        const double aic = aic_or_infinity(candidate);
        if (aic > best_aic) {
          ++k;
          break;
        }
        best = std::move(candidate);
        best_aic = aic;
        if (k == 1) break;
        --k;
      }
      continue;
    }
    std::vector<SeasonalPeriod> up = seasonal, down = seasonal;
    up[i].harmonics = kStepSearchLimit + 1;
    down[i].harmonics = kStepSearchLimit - 1;
    k = kStepSearchLimit;
    std::optional<FittedSpec> up_fit = fit_with(up);
    std::optional<FittedSpec> level_fit = fit_with(seasonal);
    std::optional<FittedSpec> down_fit = fit_with(down);
    const double up_aic = aic_or_infinity(up_fit);
    const double level_aic = aic_or_infinity(level_fit);
    const double down_aic = aic_or_infinity(down_fit);
    const double lowest = std::min({up_aic, level_aic, down_aic});
    if (lowest == down_aic) {
      best = std::move(down_fit);
      k = kStepSearchLimit - 1;
      while (true) {
        --k;
        std::optional<FittedSpec> candidate = fit_with(seasonal);
        if (aic_or_infinity(candidate) > aic_or_infinity(best)) {
          ++k;
          break;
        }
        best = std::move(candidate);
        if (k == 1) break;
      }
    } else if (lowest == level_aic) {
      best = std::move(level_fit);
    } else {
      best = std::move(up_fit);
      k = kStepSearchLimit + 1;
      while (true) {
        ++k;
        std::optional<FittedSpec> candidate = fit_with(seasonal);
        if (aic_or_infinity(candidate) > aic_or_infinity(best)) {
          --k;
          break;
        }
        best = std::move(candidate);
        if (k == max_k) break;
      }
    }
  }
  return best;
}

}  // namespace search

// TBATS: trigonometric seasonality with the given periods, which may be
// non-integer. A non-seasonal BATS search competes throughout. The
// harmonics are chosen once in the most general cell of the grid, then the
// grid is searched with the ARMA refinement per cell, and the lowest AIC
// including the non-seasonal candidate is the model.
inline FittedSpec search_tbats(const std::vector<double> &y,
                               const std::vector<double> &periods,
                               SearchOptions options) {
  if (y.empty()) throw std::invalid_argument("tbats: empty series");
  std::vector<SeasonalPeriod> seasonal;
  for (const double m : periods) {
    if (m > 1.0) seasonal.push_back({m, 1});
  }
  if (search::is_constant(y)) return search::constant_series_fit(y);
  if (!search::positive(y)) options.box_cox = false;
  FittedSpec non_seasonal = search_bats(y, {}, options);
  if (seasonal.empty()) return non_seasonal;
  if (options.trend && !*options.trend) options.damping = false;
  const std::vector<bool> box_cox_axis = search::axis(options.box_cox);
  const std::vector<bool> trend_axis = search::axis(options.trend);
  const std::vector<bool> damping_axis = search::axis(options.damping);
  const double init_lambda =
      search::any_true(box_cox_axis)
          ? search::starting_lambda(y, seasonal, options)
          : 1.0;
  const bool general_box_cox = search::any_true(box_cox_axis);
  const bool general_trend = search::any_true(trend_axis);
  const bool general_damping = search::any_true(damping_axis);
  std::optional<FittedSpec> general = search::choose_harmonics(
      y, seasonal, general_box_cox, general_trend, general_damping,
      init_lambda, options);
  std::optional<FittedSpec> best;
  if (non_seasonal.aic < search::aic_or_infinity(general)) {
    best = non_seasonal;
  } else {
    best = general;
  }
  for (const bool box_cox : box_cox_axis) {
    for (const bool trend : trend_axis) {
      for (const bool damping : damping_axis) {
        if (damping && !trend) continue;
        const bool is_general = box_cox == general_box_cox &&
                                trend == general_trend &&
                                damping == general_damping;
        std::optional<FittedSpec> cell =
            is_general ? general
                       : search::try_fit(
                             y,
                             search::cell_spec(SeasonalType::kTrigonometric,
                                               seasonal, box_cox, trend,
                                               damping, options),
                             init_lambda, options);
        search::keep_better(
            best, search::with_arma_errors(y, std::move(cell), init_lambda,
                                           options));
      }
    }
  }
  if (!best) throw std::runtime_error("tbats: unable to fit a TBATS model");
  return *best;
}

// ===========================================================================
// 9. The fitted model on the original scale
// ===========================================================================

// Quantile of the standard normal distribution, algorithm AS 241 (PPND16),
// accurate to about 1e-16 over the open unit interval.
// [[CITATION]] Wichura (1988), "Algorithm AS 241: The percentage points of
// the normal distribution", Applied Statistics 37(3), 477-484.
inline double normal_quantile(const double prob) {
  if (!(prob > 0.0 && prob < 1.0)) {
    throw std::invalid_argument("tbats: quantile probability must lie in (0, 1)");
  }
  const double q = prob - 0.5;
  if (std::abs(q) <= 0.425) {
    const double r = 0.180625 - q * q;
    return q *
           (((((((2509.0809287301226727 * r + 33430.575583588128105) * r +
                 67265.770927008700853) * r + 45921.953931549871457) * r +
               13731.693765509461125) * r + 1971.5909503065514427) * r +
             133.14166789178437745) * r + 3.387132872796366608) /
           (((((((5226.495278852545925 * r + 28729.085735721942674) * r +
                 39307.89580009271061) * r + 21213.794301586595867) * r +
               5394.1960214247511077) * r + 687.1870074920579083) * r +
             42.313330701600911252) * r + 1.0);
  }
  double r = q < 0 ? prob : 1.0 - prob;
  r = std::sqrt(-std::log(r));
  double value;
  if (r <= 5.0) {
    r -= 1.6;
    value = (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r +
                  0.24178072517745061177) * r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) * r +
              4.6303378461565452959) * r + 1.42343711074968357734) /
            (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r +
                  0.0151986665636164571966) * r + 0.14810397642748007459) * r +
                0.68976733498510000455) * r + 1.6763848301838038494) * r +
              2.05319162663775882187) * r + 1.0);
  } else {
    r -= 5.0;
    value = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +
                  0.0012426609473880784386) * r + 0.026532189526576123093) * r +
                0.29656057182850489123) * r + 1.7848265399172913358) * r +
              5.4637849111641143699) * r + 6.6579046435011037772) /
            (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r +
                  1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r +
                0.0148753612908506148525) * r + 0.13692988092273580531) * r +
              0.59983220655588793769) * r + 1.0);
  }
  return q < 0 ? -value : value;
}

struct Forecast {
  explicit Forecast(const size_t h)
      : mean(h), lower(h), upper(h), mean_model_scale(h), variance_model_scale(h) {}
  std::vector<double> mean, lower, upper;  // original scale
  std::vector<double> mean_model_scale, variance_model_scale;
};

// Forecasts h steps past the end of the fitted series with a central
// interval of `level` percent. On the model scale the h-step distribution
// is normal with the variance of section 4; mean and bounds are mapped back
// through the inverse Box-Cox transformation, the mean with the second
// order bias adjustment when asked for. A lower bound below zero is clamped
// to zero for lambda < 1, where the transformation has no inverse there.
inline Forecast forecast(const FittedSpec &fit, const size_t h,
                         const double level, const bool bias_adjust) {
  const System system(fit.spec, fit.parameters);
  const ForecastModelScale model =
      forecast_model_scale(system, fit.filter.final_state(), fit.variance, h);
  const double z = normal_quantile(0.5 + level / 200.0);
  Forecast out(h);
  out.mean_model_scale = model.mean;
  out.variance_model_scale = model.variance;
  for (size_t i = 0; i < h; ++i) {
    const double sd = std::sqrt(model.variance[i]);
    const double lower = model.mean[i] - z * sd, upper = model.mean[i] + z * sd;
    if (fit.spec.box_cox) {
      const double lambda = fit.parameters.lambda;
      out.mean[i] = bias_adjust ? box_cox::inverse_bias_adjusted(
                                      model.mean[i], lambda, model.variance[i])
                                : box_cox::inverse(model.mean[i], lambda);
      out.lower[i] = box_cox::inverse(lower, lambda);
      if (lambda < 1.0) out.lower[i] = std::max(out.lower[i], 0.0);
      out.upper[i] = box_cox::inverse(upper, lambda);
    } else {
      out.mean[i] = model.mean[i];
      out.lower[i] = lower;
      out.upper[i] = upper;
    }
  }
  return out;
}

// Sample paths on the original scale from the end of the fitted series, one
// column per path of the model-scale innovations (rows are steps), without
// bias adjustment.
inline Matrix simulate(const FittedSpec &fit, const Matrix &innovations) {
  const System system(fit.spec, fit.parameters);
  Matrix paths = simulate_model_scale(system, fit.filter.final_state(), innovations);
  if (fit.spec.box_cox) {
    for (double &v : paths.storage()) {
      v = box_cox::inverse(v, fit.parameters.lambda);
    }
  }
  return paths;
}

// Components at every observation, on the model scale: column 0 the level,
// then the slope when there is a trend, then one column per seasonal period
// holding the seasonal effect of that period (the sum of its harmonic states
// for a trigonometric block, the current dummy for a dummy block). Row t
// holds the components after observation t.
inline Matrix components(const FittedSpec &fit) {
  const StateLayout layout(fit.spec);
  const size_t n = fit.filter.states.cols();
  const size_t count = 1 + (layout.has_trend() ? 1 : 0) + layout.seasonal_count();
  Matrix out(n, count);
  for (size_t t = 0; t < n; ++t) {
    const double *x = fit.filter.states.column(t);
    size_t c = 0;
    out(t, c++) = x[StateLayout::kLevel];
    if (layout.has_trend()) out(t, c++) = x[StateLayout::kTrend];
    for (size_t i = 0; i < layout.seasonal_count(); ++i) {
      const double *block = x + layout.seasonal_offset(i);
      double effect = 0;
      if (fit.spec.seasonal_type == SeasonalType::kTrigonometric) {
        for (size_t j = 0; j < fit.spec.seasonal[i].harmonics; ++j) {
          effect += block[j];
        }
      } else {
        effect = block[0];
      }
      out(t, c++) = effect;
    }
  }
  return out;
}

// The same specification and parameters applied to another series: the
// seed states are re-estimated for it and nothing is optimised.
inline FittedSpec refit(const FittedSpec &fit, const std::vector<double> &y,
                        const bool bias_adjust) {
  if (y.empty()) throw std::invalid_argument("tbats: empty series");
  const ModelSpec &spec = fit.spec;
  const Parameters &par = fit.parameters;
  if (spec.box_cox) {
    for (const double v : y) {
      if (v <= 0) {
        throw std::invalid_argument(
            "tbats: a Box-Cox model cannot be refitted to a series that is "
            "not positive");
      }
    }
  }
  std::vector<double> seed = seed_states_for(spec, par, y);
  if (spec.box_cox) {
    for (double &v : seed) v = box_cox::transform(v, par.lambda);
  }
  const System system(spec, par);
  FilterResult result =
      filter(system, spec.box_cox ? box_cox::transform(y, par.lambda) : y, seed);
  double sse = 0;
  for (const double e : result.errors) sse += e * e;
  const double n = static_cast<double>(y.size());
  const double variance = sse / n;
  double neg2loglik = n * std::log(sse);
  if (spec.box_cox) {
    double log_sum = 0;
    for (const double v : y) log_sum += std::log(v);
    neg2loglik -= 2.0 * (par.lambda - 1.0) * log_sum;
  }
  std::vector<double> fitted = result.predictions;
  if (spec.box_cox) {
    for (double &v : fitted) {
      v = bias_adjust ? box_cox::inverse_bias_adjusted(v, par.lambda, variance)
                      : box_cox::inverse(v, par.lambda);
    }
  }
  const OptimizationResult opt{neg2loglik, 0, true};
  return FittedSpec(spec, par, seed, std::move(result), std::move(fitted),
                    variance, neg2loglik, opt);
}

}  // namespace tbats

#endif  // TBATS_TBATS_H_
