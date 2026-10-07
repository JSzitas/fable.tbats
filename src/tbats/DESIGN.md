# BATS and TBATS in dependency-free C++17: design

This document specifies the replacement of the fitting code in this package.
Today the package carries about 3100 lines of R vendored from the forecast
package plus six RcppArmadillo files, bound together by `.Call` and by
in-place mutation of R matrices. The replacement is one header,
`src/tbats/tbats.h`, with no dependency beyond the C++17 standard library and
two vendored header-only libraries (nlsolver for optimisation, tinyqr for
linear algebra). R sees it through one Rcpp glue file that converts arguments
and results and does nothing else. The header is written so that any C++
program can include it and fit, forecast and simulate these models without R.

## 1. What the models are and why they exist

BATS and TBATS are exponential smoothing state space models for a single
series with possibly several seasonal periods, with Box-Cox transformation,
damped trend and ARMA errors as optional parts. The acronym names the parts:
Box-Cox, ARMA errors, Trend, Seasonal; the leading T of TBATS stands for the
Trigonometric representation of seasonality.

They exist to forecast series whose seasonality is not a single integer
period: hourly data with daily and weekly cycles, daily data with weekly and
yearly cycles, or periods that are not integers at all (a 365.25 day year).
A seasonal dummy representation, which BATS uses, needs one state per
position in each period, so a weekly period on hourly data costs 168 states
and a yearly period on daily data costs 365. TBATS represents each seasonal
component by a few Fourier harmonics, so the same series costs a handful of
states per period, fits in a fraction of the time, and can hold non-integer
periods. The consequence for a user is that TBATS is the only model in the
fable ecosystem that fits multi-seasonal series automatically, including the
choice of how many harmonics each period needs.

[[CITATION]] De Livera, Hyndman and Snyder (2011), "Forecasting time series
with complex seasonal patterns using exponential smoothing", Journal of the
American Statistical Association 106(496), 1513-1527. All equations below are
from this paper unless stated otherwise.

## 2. Why rewrite rather than keep the vendored code

The vendored implementation works and is the reference for the port. It is
replaced because its structure, not its mathematics, is the liability:

- The state layout (which position in the state vector holds which component)
  is recomputed by arithmetic in more than a dozen places, in R and in C++.
  One such recomputation in the forecast code reads the trend coefficient
  from a field that does not exist, so the gain vector used for TBATS
  prediction intervals omits the trend and shifts the seasonal entries. A
  single layout registry makes that class of error impossible.
- Every matrix is dense, and during optimisation the matrices are updated
  in place entry by entry. That update skips the coupling of the second and
  later dummy seasonal blocks to the ARMA states, so a BATS model with two
  seasonal periods and ARMA errors is optimised against an objective that is
  not its likelihood; the frozen multi-seasonal BATS fit carries an AIC
  six units off its own errors for that reason. The transition matrix is in
  fact sparse with a rank-one coupling, so the port never forms it in the
  kernel and the per-step cost falls from the square of the state dimension
  to the dimension itself.
- The optimisation, the ARMA order selection and the linear algebra are R's,
  which ties the model to R. With the fitting in a self-contained header the
  same code can serve other hosts later.
- The R object built for fable is assembled by hand in four places, which is
  how refit ended up with empty residuals and the BATS refit ended up calling
  the TBATS fitter. One constructor on the R side removes both.

The forecast package's behaviour is not reproduced bit for bit: the optimiser
changes, so parameter estimates and model orders will differ within
tolerance. Section 8 defines what must match and how closely.

## 3. The model

Notation. The observed series is y_1, ..., y_n. The Box-Cox parameter is w;
the transformed series is y^(w) = (y^w - 1)/w for w != 0 and log y for w = 0.
The components at time t are the level l_t, the trend b_t with damping
parameter phi in [0.8, 1] (phi = 1 means no damping), one seasonal component
s_t^(i) per seasonal period m_i for i = 1..M, and a disturbance d_t that
follows an ARMA(p, q) process driven by innovations e_t ~ N(0, sigma^2) with
autoregressive coefficients ar_1..ar_p and moving average coefficients
ma_1..ma_q. The smoothing parameters are alpha (level), beta (trend) and,
per seasonal component, gamma_i (BATS) or the pair gamma1_i, gamma2_i
(TBATS).

Observation and component recursions (paper, section 2):

    y_t^(w) = l_{t-1} + phi b_{t-1} + sum_i s_{t-1}^(i) + d_t
    l_t     = l_{t-1} + phi b_{t-1} + alpha d_t
    b_t     = phi b_{t-1} + beta d_t
    d_t     = sum_{j=1}^p ar_j d_{t-j} + sum_{j=1}^q ma_j e_{t-j} + e_t

BATS seasonal component: s_t^(i) = s_{t-m_i}^(i) + gamma_i d_t, held as an
m_i-vector that shifts cyclically each step.

TBATS seasonal component with k_i harmonics, angles lambda_j = 2 pi j / m_i:

    s_t^(i)   = sum_{j=1}^{k_i} s_{j,t}
    s_{j,t}   =  s_{j,t-1} cos lambda_j + s*_{j,t-1} sin lambda_j + gamma1_i d_t
    s*_{j,t}  = -s_{j,t-1} sin lambda_j + s*_{j,t-1} cos lambda_j + gamma2_i d_t

Each harmonic is a 2 by 2 rotation; the component is 2 k_i states. For a
period of exactly 2 the single harmonic has angle pi, so the rotation is a
sign flip; the forecast package replaces its cosine by zero in that one case,
which lets the component decay to nothing within a step. The port keeps the
rotation. No frozen fixture has a period of 2, so this is the one place the
port knowingly departs from the reference kernel.

### 3.1 State space form and the structure the kernel exploits

Stack the states as x_t = (l_t, b_t, seasonal blocks, d_t..d_{t-p+1},
e_t..e_{t-q+1}). The dimension is 1 + [trend] + tau + p + q with tau the
total seasonal length (sum of 2 k_i for TBATS, sum of m_i for BATS). The
model is the innovations state space form

    y_t^(w) = w' x_{t-1} + e_t
    x_t     = F x_{t-1} + g e_t

with w the observation weights (1, phi, seasonal weights, ar, ma), g the
gains (alpha, beta, gamma entries, 1 at the head of the AR block, 1 at the
head of the MA block) and F the transition matrix of the paper.

F is a structural part plus a rank-one coupling. Define c as the vector
holding ar in the AR block and ma in the MA block, zeros elsewhere; u = w - c
as the observation weights without the ARMA part; g_d as g with the MA
block's 1 removed; and g_e as the unit vector at the head of the MA block.
Then

    d_t = c' x_{t-1} + e_t
    y_t^(w) = u' x_{t-1} + d_t
    x_t = F_struct x_{t-1} + g_d d_t + g_e e_t

where F_struct applies the level and trend update, one rotation (TBATS) or
cyclic shift (BATS) per seasonal block, and a one-position shift in the AR
and MA blocks with a zero head. The dense identities are F = F_struct + g_d c'
and g = g_d + g_e. The kernel never forms F: one step costs a handful of
multiply-adds per state. The matrix D = F - g w', which appears in the
admissibility check and in seed-state estimation, is formed densely only for
the eigenvalue test; the seed-state recursion uses v' D = v' F_struct +
(v' g_d) c' - (v' g) w', which is structured as well.

### 3.2 Seed states

The state x_0 before the first observation is estimated by least squares
(paper, section 4.2). Writing the recursion as x_t = D x_{t-1} + g y_t^(w),
the one-step error is affine in x_0: e_t = ytilde_t - wtilde_t' x_0, where
ytilde_t is the error sequence obtained with x_0 = 0 and wtilde_t' = w' D^{t-1}.
Regressing ytilde on wtilde without intercept gives x_0. The ARMA seed states
are fixed at zero and their columns dropped. For BATS the seasonal dummy
columns are collinear with the level within each block, and nested or
factor-sharing periods are collinear across blocks; one column is dropped per
block (more for nested periods, following the forecast package's mask) and
the dropped seasonal states are restored so that each block sums to zero.
TBATS harmonics need no such treatment.

### 3.3 Likelihood, information criterion and admissibility

With sigma^2 concentrated out, the profile log-likelihood is, up to constants,

    -2 log L = n log( sum_t e_t^2 ) - 2 (w - 1) sum_t log y_t

where the second term is the Jacobian of the Box-Cox transformation and is
absent without it (paper, equation 12). sigma^2 = sum e_t^2 / n. The AIC used
for every selection step is -2 log L + 2 (number of parameters + number of
seed states), seed states counted as parameters as the paper does.

A parameter vector is admissible when: w lies strictly inside the user's
bounds; phi lies in [0.8, 1]; the AR polynomial 1 - ar_1 z - ... - ar_p z^p
and the MA polynomial 1 + ma_1 z + ... + ma_q z^q have all roots of modulus
above 1.01 (stationarity and invertibility with a one percent margin); and
every eigenvalue of D has modulus below 1.01 (the forecastability condition
of [[CITATION]] Hyndman, Koehler, Ord and Snyder (2008), "Forecasting with
exponential smoothing: the state space approach", chapter 10). The objective
returns a large finite penalty for an inadmissible vector. It is finite
rather than infinite because the Nelder-Mead stopping rule takes a standard
deviation over the simplex values.

### 3.4 Optimisation

The parameters (w, alpha, phi, beta, gammas, ar, ma; each present only when
the specification has it) are packed into one vector and scaled per entry
before optimisation, with the scales of the forecast package kept as named
constants. The default optimiser is Nelder-Mead from nlsolver. The optimiser
is a runtime choice: an enumeration in the fit options is dispatched once,
at setup, to the nlsolver solver class; the likelihood functor is the same
for all of them (any solver that accepts a functor and a start vector can be
added by one case in that switch).

### 3.5 Model search

BATS(y, options): Box-Cox is disabled when any y <= 0. Each cell of the grid
Box-Cox in {off, on} x trend in {off, on} x damping in {off, on} (damping only
with trend; the user may pin any axis) is fitted without ARMA errors, the
ARMA orders are ranked on that fit's errors (section 4), the cell is
refitted with the best order and with the most parsimonious order within
two AIC units of it, and the lowest AIC of the three fits is kept. The
forecast package refits the best order only; on short series the top
orders are near-ties by the residual fit and the parsimonious one is often
the better model once its coefficients are re-estimated, which is what
decided the Lynx search in the port's favour. The best cell by AIC is the
model. A fit that fails, or whose parameters make the eigenvalue iteration
of the admissibility check diverge, is a candidate with infinite AIC.

TBATS(y, periods, options): a constant series yields the degenerate model
(level equal to the series, zero variance). Otherwise a non-seasonal BATS
candidate is fitted first. The harmonic counts k_i are then chosen by AIC
using the most general cell the options allow: the maximum is
floor((m_i - 1)/2), reduced where a lower harmonic order would duplicate the
frequencies of an earlier period; when that maximum is at most 6 the search
starts there and steps down while the AIC improves, otherwise it compares
k = 5, 6, 7 and walks in the improving direction. With k fixed, the grid of
section 3.5 runs as for BATS, with ARMA refinement per cell, and the best AIC
across the grid and the non-seasonal candidate is the model. The search is
sequential; fable parallelises across series.

### 3.6 Refit, forecast, simulation, components

Refit to a new series with the same parameters re-estimates the seed states
by the regression of section 3.2 and reruns the kernel; nothing is optimised.
(The forecast package reuses the old seed states instead, which is only
right when the new series starts at the same origin.) Refit with
re-estimation repeats the full search with the stored options.

Forecast means come from the kernel run forward from x_n with e_t = 0. The
h-step variance is sigma^2 (1 + sum_{j=1}^{h-1} c_j^2) with c_j = w' F^{j-1} g
(paper, section 5), computed by iterating the structured transition on g.
Means and interval bounds are mapped back through the inverse Box-Cox
transformation; the mean may be bias adjusted by the forecast package's
second-order expansion when the user asks for it.

Simulation runs the same kernel with supplied innovations, one column per
path, and maps the paths back without bias adjustment. Components are read
off the stored state history: level, trend, and each seasonal component as
its block's observation weights applied to the block.

## 4. ARMA order selection by exact likelihood

The ARMA orders are chosen on the errors of the fit without ARMA terms by
minimising AIC over p, q in 0..5, each candidate fitted by exact Gaussian
maximum likelihood with a mean term, as R's `arima` does. The coefficients
found here are only starting points: they are re-estimated inside the TBATS
likelihood. Exact likelihood is preferred over conditional sum of squares
because conditional sum of squares discards the first p observations and
treats the pre-sample moving-average terms as zero, which biases the order
choice on short series and on series with strong moving-average structure,
and because the same code is wanted for other models later.

The implementation is lifted from blaze's Kalman ARIMA and reduced to what is
needed: ARMA(p, q) with a mean, no differencing, no seasonal part, no
regressors. The state space form and the initial state covariance follow
[[CITATION]] Gardner, Harvey and Phillips (1980), "An algorithm for exact
maximum likelihood estimation of autoregressive-moving average models by
means of Kalman filtering", Applied Statistics 29(3), 311-322, as in R's
`arima`. Parameters are optimised in the partial-autocorrelation
transformed space so that every candidate is stationary and invertible
([[CITATION]] Jones (1980), "Maximum likelihood fitting of ARMA models to
time series with missing observations", Technometrics 22(3), 389-395).
Everything that tied the blaze code to its host is left behind: Eigen
vectors and matrices become column-major standard vectors, the cppoptlib
base class becomes an nlsolver functor, presence flags carried as template
booleans become plain data, the default scalar becomes double, the
two-phase construction (an empty object later filled by fit) becomes a
constructor, and console output is removed. The cost is one Kalman recursion
per objective evaluation instead of one conditional recursion; the forecast
package pays the same today.

## 5. Linear algebra: what tinyqr provides and what it must gain

tinyqr supplies the least-squares solve for seed states (`lm` on a
column-major design) and the QR decomposition behind it. Its eigenvalue
routine is unshifted QR iteration that reads the diagonal, which is correct
only for symmetric matrices: on a plane rotation it returns zeros, on a
matrix with eigenvalues 0.9 +- 0.4i it returns 0.9 twice, and on a real
upper-triangular matrix with eigenvalues 2 and 1 it returns 2.17 and 0.83.
The admissibility check needs the moduli of all eigenvalues of D, a real
non-symmetric matrix whose seasonal rotations produce complex pairs.

tinyqr therefore gains a general real eigenvalue routine before the port
depends on it: reduction to Hessenberg form by orthogonal similarity,
Francis double-shift QR iteration with deflation, and 2 by 2 blocks read
out as complex pairs ([[CITATION]] Golub and Van Loan, "Matrix
Computations", 4th edition, section 7.5). Its tests use spectra known in
closed form: rotations, complex pairs, triangular matrices, companion
matrices of polynomials with known roots, and random orthogonal similarity
transforms of diagonal matrices. The existing unshifted routine is either
corrected for non-symmetric input or documented as symmetric-only. The AR
and MA root checks of section 3.3 use the same routine on companion
matrices, so one eigenvalue primitive serves both tests. nlsolver bundles an
older tinyqr and uses two of its symbols (`lm` on vectors and `QRSolver`);
the updated tinyqr keeps both so one vendored copy serves both libraries.

## 6. Shape of the header

`tbats.h` is a single file in dependency order, each section usable without
the ones below it:

1. A column-major matrix type with its dimensions fixed at construction, and
   the few dense operations the admissibility check needs.
2. Box-Cox transform, inverse with optional bias adjustment, and Guerrero's
   choice of w by minimising the coefficient of variation of the transformed
   subseries over a bounded interval with Brent's method ([[CITATION]]
   Guerrero (1993), "Time-series analysis supported by power
   transformations", Journal of Forecasting 12(1), 37-48). The loss is
   blaze's; blaze drives it with a root finder, which is not the right tool
   for a minimum.
3. The specification (what is fixed before optimisation: Box-Cox and its
   bounds, trend, damping, seasonal periods with harmonics or dummies, ARMA
   orders), the state layout as the single registry of block offsets built
   from it, and the parameter set with pack, unpack and scaling. Absent
   parts are cheap fixed members, not optionals: no damping is phi = 1.
4. The system (u, c, g_d, g_e, the structural blocks) born from
   specification and parameters, and the one kernel step used by every
   consumer. D is a method that materialises on demand.
5. Seed-state estimation.
6. The likelihood functor with admissibility and scaling, and the
   optimiser dispatch.
7. ARMA exact likelihood and order selection, in its own namespace so it can
   be lifted out unchanged.
8. The fit of one specification, and the BATS and TBATS searches.
9. The fitted model: parameters, seed and final states, state history,
   errors, fitted values, variance, log-likelihood, AIC; with forecast,
   simulate, components and refit.

Errors are thrown at the point of failure. The header is templated on
nothing: the scalar is double throughout.

The header includes `nlsolver.h` by bare name, and nlsolver includes its
`tinyqr.h` from its own directory, so a consumer copies the three headers
side by side or adds the directory holding the two vendored ones to the
include path. In this package they live under `src/third_party` and the
build adds that directory.

## 7. The R side

The fable wrappers keep their specials and methods. The glue file converts
the series and options to C++, calls the fitter, and returns a list. Trimming
to the longest run of non-missing values stays in R. Seasonality detection
(`find_seasonalities`) stays in R. One R constructor builds the fitted object
for training and both refit paths, so residuals, fitted values and the
summary string come from one place; refit dispatches on the model class.
`generate()` calls the C++ simulation with the innovations fabletools or the
user supplies. The vendored R file, the Armadillo sources and the
RcppArmadillo dependency are deleted when the glue is switched over; the
build moves to C++17.

## 8. Verification

Exact agreement with the forecast package is not the goal: the optimiser
differs, so the estimates differ. The tests separate what is deterministic
given parameters from what depends on optimisation, and the deterministic
parts are compared against fixtures frozen from the current implementation
before any of it is deleted. Every comparison carries a tolerance.

Fixtures: for a handful of series (yearly, monthly, a multi-seasonal
synthetic series, a series with a zero to exercise the Box-Cox disabling,
and a constant series), and for each specification visited, the fixture
holds the parameter vector, the seed states, the dense w, F and g, the error
and fitted sequences, forecast means and standard deviations, the
log-likelihood and the AIC. For ARMA order selection the fixtures hold the
errors of the first fit and R's exact log-likelihood for a grid of orders.

What is tested, and what each test can see:

- Kernel against fixtures at fixed parameters: the structured step must
  reproduce the dense w, F, g (assembled from the structured parts), the
  errors and the fitted values. This catches any layout, sign or coupling
  error, including the trend-coefficient defect, because it is independent
  of optimisation. It cannot see an error in the likelihood formula.
- Seed states at fixed parameters: x_0 against the fixture. Different QR
  implementations give slightly different least-squares solutions on
  ill-conditioned designs, so the tolerance is looser than the kernel's.
- Likelihood at fixed parameters: -2 log L and AIC against the fixture. This
  is the test that sees a wrong likelihood or a wrong parameter count.
- Fit: for each fixture specification the port's optimised -2 log L must be
  no worse than the fixture's beyond a tolerance. A better value is
  recorded, not failed. This sees optimiser regressions and search-logic
  errors; it cannot see whether the search finds the global optimum, which
  the reference does not guarantee either.
- ARMA exact likelihood against R's `arima` log-likelihood on the fixture
  errors across the order grid. This sees initial-covariance and
  parameter-transform errors directly. The selected order is then checked
  to be within a small AIC margin of the port's own best, since near-ties
  may legitimately flip.
- Eigenvalues: tinyqr's own tests, section 5.
- Behaviour: the existing tests (simulation with zero innovations equals
  the forecast mean; replaying the fitted errors reproduces the data; path
  counts and keys; refit; glance) remain and are extended to a seasonal
  series. These see plumbing errors between R and C++.

Not covered by design: statistical quality against a ground truth, which
is a property of the model, not of the port.

## 9. Order of work

Each step lands with its tests passing before the next starts.

0. Freeze the fixtures from the current implementation.
1. tinyqr: eigenvalue routine, tests, vendor the updated copy here together
   with nlsolver.
2. Header sections 1 to 4; kernel and matrix tests against fixtures.
3. Sections 5, 6 and the single-specification fit; seed-state, likelihood
   and fit tests.
4. Section 7; ARMA tests.
5. Section 8; fit tests on the full searches.
6. Section 9, the glue, the R constructor and method changes; behaviour
   tests; delete the vendored R and Armadillo code; C++17 build; drop
   RcppArmadillo.
7. README and NEWS.

The draft headers `BATS.h` and `matrices.h` in this directory are superseded
by this design and are removed when `tbats.h` lands.
