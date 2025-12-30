# Statistical Analysis and Random Number Generation

Comprehensive toolkit for descriptive statistics, probability distributions, and high-quality random number generation.

## Overview

Statistical analysis is fundamental across all quantitative disciplines:
- **Science**: Experimental data analysis, hypothesis testing, error propagation
- **Engineering**: Quality control, reliability analysis, signal processing
- **Finance**: Risk assessment, portfolio optimization, Monte Carlo valuation
- **Machine Learning**: Data preprocessing, distribution fitting, sampling methods
- **Simulation**: Stochastic modeling, Monte Carlo integration, agent-based systems

This library provides **two complementary toolsets**:
1. **Descriptive Statistics**: Mean, variance, moments, skewness, kurtosis
2. **Probability & Random Generation**: Distributions (pdf/cdf/quantile) + high-quality random deviates

## Quick Reference

### Descriptive Statistics Functions

| Function | Purpose | Output | Use Case |
|----------|---------|--------|----------|
| **Mean** | Arithmetic mean | μ = Σxᵢ/n | Central tendency |
| **Variance** | Sample variance | σ² (unbiased) | Spread measurement |
| **StdDev** | Standard deviation | σ | Dispersion analysis |
| **Covariance** | Linear relationship | Cov(X,Y) | Joint variability |
| **Correlation** | Standardized covariance | ρ ∈ [-1,1] | Association strength |

### Probability Distributions

| Distribution | Parameters | Support | Key Use |
|--------------|------------|---------|---------|
| **Normal** | μ (mean), σ (std dev) | ℝ | Central limit, z-tests |
| **T-Distribution** | ν (degrees of freedom) | ℝ | Small-sample inference |
| **Chi-Square** | k (degrees of freedom) | [0, ∞) | Variance tests, goodness-of-fit |
| **F-Distribution** | d₁, d₂ (df) | [0, ∞) | ANOVA, variance ratios |

### Hypothesis Tests

| Test | Purpose | Assumptions | Output |
|------|---------|-------------|--------|
| **OneSampleTTest** | μ = μ₀? | Normal or n≥30 | t-statistic, p-value, CI |
| **TwoSampleTTest** | μ₁ = μ₂? | Independent, normal | Welch's t-test result |
| **PairedTTest** | μ_diff = 0? | Paired, normal diffs | More powerful for pairs |
| **ChiSquareGoodnessOfFit** | Fits distribution? | Expected freq ≥ 5 | χ² statistic, p-value |
| **ChiSquareTestOfIndependence** | Variables independent? | Expected counts ≥ 5 | χ² test on contingency table |
| **OneWayANOVA** | μ₁=μ₂=...=μₖ? | Normal, equal variance | F-statistic, p-value |

### Confidence Intervals

| Function | Estimates | Method | Use Case |
|----------|-----------|--------|----------|
| **ConfidenceIntervalMean** | Population μ | t-distribution | Single sample mean |
| **ConfidenceIntervalMeanDifference** | μ₁ - μ₂ | Welch's method | Compare two means |
| **ConfidenceIntervalProportion** | Population p | Normal approx | Binomial proportion |
| **ConfidenceIntervalProportionDifference** | p₁ - p₂ | Normal approx | Compare proportions |
| **ConfidenceIntervalPairedDifference** | μ_diff | t-distribution | Paired observations |

### Random Number Generators

| Generator | Distribution | Algorithm | Speed | Quality |
|-----------|--------------|-----------|-------|---------|
| **ExponentialDeviate** | Exponential(λ) | Inverse transform | ⭐⭐⭐⭐⭐ | Exact |
| **LogisticDeviate** | Logistic(μ,σ) | Inverse CDF | ⭐⭐⭐⭐⭐ | Exact |
| **NormalDeviateBoxMuller** | Normal(μ,σ) | Box-Muller pairs | ⭐⭐⭐⭐ | Exact |
| **NormalDeviate** | Normal(μ,σ) | **Leva ratio-of-uniforms** | ⭐⭐⭐⭐⭐ | Exact |
| **CauchyDeviate** | Cauchy(μ,σ) | Ratio method | ⭐⭐⭐⭐ | Exact |
| **GammaDeviate** | Gamma(α,β) | Marsaglia-Tsang | ⭐⭐⭐⭐ | High |
| **PoissonDeviate** | Poisson(λ) | Adaptive (direct/PTRS) | ⭐⭐⭐⭐ | Exact |
| **BinomialDeviate** | Binomial(n,p) | Adaptive (3 methods) | ⭐⭐⭐⭐ | Exact |

**All generators use Mersenne Twister 64-bit** (std::mt19937_64) with period 2¹⁹⁹³⁷−1.

## Mathematical Background

### Descriptive Statistics

**Sample Mean** (first moment):
```
μ = (1/n) · Σᵢ xᵢ
```

**Sample Variance** (second central moment):
```
σ² = [Σᵢ(xᵢ - μ)²] / (n-1)
```
Uses **Bessel's correction** (n-1) for unbiased estimator.

**Standard Deviation**:
```
σ = √(σ²)
```

**Mean Absolute Deviation**:
```
MAD = (1/n) · Σᵢ |xᵢ - μ|
```

**Skewness** (third standardized moment):
```
γ₁ = [Σᵢ(xᵢ - μ)³] / (n·σ³)
```
- γ₁ > 0: Right-skewed (long right tail)
- γ₁ < 0: Left-skewed (long left tail)
- γ₁ ≈ 0: Symmetric

**Kurtosis** (fourth standardized moment, excess):
```
γ₂ = [Σᵢ(xᵢ - μ)⁴] / (n·σ⁴) - 3
```
- γ₂ > 0: Heavy-tailed (leptokurtic)
- γ₂ < 0: Light-tailed (platykurtic)
- γ₂ ≈ 0: Normal-like tails (mesokurtic)

### Probability Theory

**Probability Density Function (PDF)**: f(x)
- For continuous X: P(a ≤ X ≤ b) = ∫ₐᵇ f(x)dx
- ∫₋∞^∞ f(x)dx = 1

**Cumulative Distribution Function (CDF)**: F(x) = P(X ≤ x)
```
F(x) = ∫₋∞ˣ f(t)dt
```

**Quantile Function (Inverse CDF)**: F⁻¹(p)
```
F⁻¹(p) = x  such that  F(x) = p
```

**Random Variate Generation**:
- **Inverse transform**: X = F⁻¹(U) where U ~ Uniform(0,1)
- **Rejection sampling**: Generate candidates, accept with probability ∝ f(x)
- **Transformation methods**: Box-Muller, ratio-of-uniforms

---

## Descriptive Statistics Functions

### Avg - Arithmetic Mean

**Purpose**: Compute sample mean (first moment).

**Function**:
```cpp
static Real Avg(const Vector<Real>& data)
```

**Algorithm**:
```
μ = (1/n) · Σᵢ₌₁ⁿ xᵢ
```

**Complexity**: O(n)

**Use When**:
- Need central tendency measure
- Data approximately symmetric
- No extreme outliers (otherwise consider median)

**Properties**:
- **Minimizes**: Sum of squared deviations Σ(xᵢ - μ)²
- **Affected by outliers**: Single extreme value can shift mean significantly
- **Linear**: E[aX + b] = a·E[X] + b

---

### AvgVar - Mean and Variance

**Purpose**: Compute both mean and sample variance in single pass.

**Function**:
```cpp
static void AvgVar(const Vector<Real>& data, Real& outAvg, Real& outVar)
```

**Algorithm** (two-pass for numerical stability):
```
Pass 1: μ = Avg(data)
Pass 2: 
  S = Σᵢ(xᵢ - μ)²
  ep = Σᵢ(xᵢ - μ)        // Correction term
  σ² = (S - ep²/n) / (n-1)
```

**Why Correction Term?**
Compensates for floating-point rounding errors in mean computation. Ensures variance never negative due to numerical precision.

**Bessel's Correction**: Division by (n-1) instead of n provides **unbiased estimator** of population variance.

**Complexity**: O(n) (two passes)

**Use When**:
- Need both mean and variance
- Numerical stability critical
- Sample size small (n < 100), where bias matters

---

### AvgStdDev - Mean and Standard Deviation

**Purpose**: Compute mean and standard deviation (square root of variance).

**Function**:
```cpp
static void AvgStdDev(const Vector<Real>& data, Real& outAvg, Real& outStdDev)
```

**Algorithm**:
```cpp
AvgVar(data, outAvg, var);
outStdDev = sqrt(var);
```

**Standard Deviation Properties**:
- **Same units** as data (variance has squared units)
- **68-95-99.7 rule** for normal distributions:
  - 68% within μ ± σ
  - 95% within μ ± 2σ
  - 99.7% within μ ± 3σ

**Use When**:
- Need interpretable spread measure
- Comparing variability across datasets
- Constructing confidence intervals

---

### Moments - Complete Distribution Shape

**Purpose**: Compute all four moments (mean, MAD, variance, skewness, kurtosis).

**Function**:
```cpp
static void Moments(const Vector<Real>& data, Real& ave, Real& adev, Real& sdev, 
                   Real& var, Real& skew, Real& curt)
```

**Parameters**:
- **ave**: Mean μ
- **adev**: Mean absolute deviation
- **sdev**: Standard deviation σ
- **var**: Variance σ²
- **skew**: Skewness γ₁
- **curt**: Excess kurtosis γ₂

**Algorithm**:
```
1. Compute mean: μ = Σxᵢ/n

2. Single-pass moment accumulation:
   For each xᵢ:
     s = xᵢ - μ
     adev += |s|
     var  += s²
     skew += s³
     curt += s⁴
     ep   += s      // Correction

3. Normalize:
   adev = adev / n
   var  = (var - ep²/n) / (n-1)
   sdev = √var
   skew = (skew/n) / (var · sdev)
   curt = (curt/n) / var² - 3
```

**Interpretation**:

**Skewness**:
- **γ₁ > 0.5**: Moderate right skew
- **γ₁ > 1.0**: Strong right skew
- **γ₁ < -0.5**: Moderate left skew
- **|γ₁| < 0.5**: Approximately symmetric

**Kurtosis**:
- **γ₂ > 0**: Heavy tails (more outliers than normal)
- **γ₂ < 0**: Light tails (fewer outliers than normal)
- **γ₂ ≈ 0**: Normal-like tail behavior

**Use When**:
- Complete distribution characterization needed
- Checking normality assumptions
- Identifying outliers and tail behavior
- Comparing distribution shapes

**Limitations**:
- Requires σ² > 0 (throws exception if all values identical)
- Sensitive to outliers (especially kurtosis)
- Sample size n ≥ 4 recommended for reliable kurtosis

---

## Probability Distributions

### Cauchy Distribution

**Mathematical Definition**:
```
f(x) = 1 / (π·σ·[1 + ((x-μ)/σ)²])
```

**Parameters**:
- **μ**: Location parameter (peak location, **not mean!**)
- **σ**: Scale parameter (half-width at half-maximum)

**Key Properties**:
- **No defined mean or variance** (integrals diverge)
- **Heavy tails**: P(|X| > x) ~ 1/x (power law)
- **Stable distribution**: Sums of Cauchy variates are Cauchy
- **Ratio interpretation**: If X, Y ~ Normal(0,1), then X/Y ~ Cauchy(0,1)

**Physics Applications**:
- Resonance line shapes (Lorentzian profile)
- Spectral line broadening
- Scattering theory

**Class**:
```cpp
struct CauchyDistribution {
    Real mu;    // Location
    Real sigma; // Scale (must be > 0)
    
    CauchyDistribution(Real location = 0.0, Real scale = 1.0);
    
    Real pdf(Real x) const;
    Real cdf(Real x) const;
    Real inverseCdf(Real p) const;  // Quantile function
};
```

**PDF**:
```cpp
Real pdf(Real x) const {
    Real z = (x - mu) / sigma;
    return 1.0 / (Constants::PI * sigma * (1.0 + z*z));
}
```

**CDF**:
```cpp
Real cdf(Real x) const {
    return 0.5 + std::atan2(x - mu, sigma) / Constants::PI;
}
```

**Quantile** (inverse CDF):
```cpp
Real inverseCdf(Real p) const {
    return mu + sigma * std::tan(Constants::PI * (p - 0.5));
}
```

**Special Values**:
- Median: μ
- Mode: μ
- IQR (interquartile range): 2σ

---

### Exponential Distribution

**Mathematical Definition**:
```
f(x) = λ · exp(-λx)    for x ≥ 0
```

**Parameters**:
- **λ**: Rate parameter (λ = 1/mean, must be > 0)

**Key Properties**:
- **Memoryless**: P(X > s+t | X > s) = P(X > t)
- **Mean**: 1/λ
- **Variance**: 1/λ²
- **Connection to Poisson**: If events occur at rate λ, waiting time ~ Exponential(λ)

**Applications**:
- Radioactive decay
- Component failure times (reliability engineering)
- Queue service times
- Time between arrivals (Poisson process)

**Class**:
```cpp
struct ExponentialDistribution {
    Real lambda;  // Rate parameter
    
    ExponentialDistribution(Real rate);
    
    Real pdf(Real x) const;
    Real cdf(Real x) const;
    Real inverseCdf(Real p) const;
    Real mean() const;      // = 1/λ
    Real variance() const;  // = 1/λ²
};
```

**PDF**:
```cpp
Real pdf(Real x) const {
    if (x < 0.0) throw StatisticsError("x must be non-negative");
    return lambda * std::exp(-lambda * x);
}
```

**CDF**:
```cpp
Real cdf(Real x) const {
    if (x < 0.0) throw StatisticsError("x must be non-negative");
    return 1.0 - std::exp(-lambda * x);
}
```

**Quantile** (inverse CDF):
```cpp
Real inverseCdf(Real p) const {
    if (p < 0.0 || p >= 1.0) throw StatisticsError("p must be in [0,1)");
    return -std::log(1.0 - p) / lambda;
}
```

**Special Values**:
- Mode: 0
- Median: ln(2)/λ ≈ 0.693/λ

---

### Logistic Distribution

**Mathematical Definition**:
```
f(x) = exp(-z) / (σ·(1 + exp(-z))²)
where z = (x - μ)/σ
```

**Parameters**:
- **μ**: Location parameter (mean and median)
- **σ**: Scale parameter (related to variance by Var = (πσ)²/3)

**Key Properties**:
- **Mean**: μ
- **Variance**: (πσ)²/3
- **CDF is logistic function**: F(x) = 1/(1 + exp(-(x-μ)/σ))
- **Symmetric** around μ
- **Heavier tails** than normal distribution

**Applications**:
- **Logistic regression**: Sigmoid activation function
- **Neural networks**: Neuron activation
- **Growth models**: Population dynamics
- **Epidemiology**: Disease spread curves

**Class**:
```cpp
struct LogisticDistribution {
    Real mu;     // Location (mean/median)
    Real sigma;  // Scale
    
    LogisticDistribution(Real location = 0.0, Real scale = 1.0);
    
    Real pdf(Real x) const;
    Real cdf(Real x) const;
    Real inverseCdf(Real p) const;
    Real mean() const;      // = μ
    Real variance() const;  // = (πσ)²/3
};
```

**PDF**:
```cpp
Real pdf(Real x) const {
    Real z = (x - mu) / sigma;
    Real exp_neg_z = std::exp(-z);
    Real denom = 1.0 + exp_neg_z;
    return exp_neg_z / (sigma * denom * denom);
}
```

**CDF** (numerically stable):
```cpp
Real cdf(Real x) const {
    Real z = (x - mu) / sigma;
    Real exp_z = std::exp(-std::abs(z));
    
    if (z >= 0.0)
        return 1.0 / (1.0 + exp_z);
    else
        return exp_z / (1.0 + exp_z);
}
```

**Quantile** (inverse CDF):
```cpp
Real inverseCdf(Real p) const {
    if (p <= 0.0 || p >= 1.0) throw StatisticsError("p must be in (0,1)");
    return mu + sigma * std::log(p / (1.0 - p));
}
```

**Special Values**:
- Mode: μ
- Median: μ
- IQR: 2σ·ln(3) ≈ 2.197σ

---

## Random Number Generators

### ExponentialDeviate - Exponential Random Numbers

**Algorithm**: Inverse transform method
```
X = -ln(U)/λ    where U ~ Uniform(0,1)
```

**Class**:
```cpp
class ExponentialDeviate {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;
    Real _lambda;
    
public:
    ExponentialDeviate(Real rate, uint64_t seed = std::random_device{}());
    
    Real generate();
    Real operator()();  // Alias for generate()
};
```

**Usage**:
```cpp
ExponentialDeviate expGen(0.5);  // Mean = 1/0.5 = 2.0

Real sample = expGen.generate();
Real sample2 = expGen();  // Same as generate()
```

**Performance**: **O(1)** per sample, very fast (one log, one division)

**Quality**: **Exact** - mathematically perfect transformation

**Thread Safety**: Each instance has own RNG state - safe if not shared

---

### LogisticDeviate - Logistic Random Numbers

**Algorithm**: Inverse CDF method
```
X = μ + σ·ln(U/(1-U))    where U ~ Uniform(0,1)
```

**Class**:
```cpp
class LogisticDeviate {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;
    Real _mu, _sigma;
    
public:
    LogisticDeviate(Real location, Real scale, uint64_t seed = std::random_device{}());
    
    Real generate();
    Real operator()();
};
```

**Implementation** (avoids U=0 and U=1):
```cpp
Real generate() {
    Real u;
    do {
        u = _uniform(_gen);
    } while (u == 0.0 || u == 1.0);  // Avoid log(0) and division by 0
    return _mu + _sigma * std::log(u / (1.0 - u));
}
```

**Performance**: **O(1)**, fast

**Quality**: **Exact** transformation

---

### NormalDeviateBoxMuller - Normal via Box-Muller

**Algorithm**: Box-Muller transform generates **pairs** of independent normals
```
Given U₁, U₂ ~ Uniform on unit circle (V₁²+V₂² < 1):
  R² = V₁² + V₂²
  fac = √(-2·ln(R²)/R²)
  X = V₁ · fac,  Y = V₂ · fac
Then X, Y ~ Normal(0,1) independently
```

**Class**:
```cpp
class NormalDeviateBoxMuller {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;  // Uniform(-1, 1)
    Real _mu, _sigma;
    Real _storedValue;  // Cache second value from pair
    bool _hasStored;
    
public:
    NormalDeviateBoxMuller(Real mean, Real stddev, uint64_t seed = std::random_device{}());
    
    Real generate();
    Real operator()();
};
```

**Implementation** (Marsaglia polar method variant):
```cpp
Real generate() {
    if (_hasStored) {
        _hasStored = false;
        return _mu + _sigma * _storedValue;  // Return cached value
    }
    
    Real v1, v2, rsq, fac;
    do {
        v1 = _uniform(_gen);   // Uniform(-1, 1)
        v2 = _uniform(_gen);
        rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);  // Rejection: keep if in unit circle
    
    fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
    _storedValue = v1 * fac;  // Cache for next call
    _hasStored = true;
    return _mu + _sigma * v2 * fac;
}
```

**Performance**: 
- Average **2.3 uniform samples** per normal (1.15 per pair)
- **Caches one value** - amortized cost very low

**Quality**: **Exact** - mathematically rigorous transformation

**Use When**:
- Need pairs of normals
- Classical method preferred
- Don't need absolute fastest single-value generation

---

### NormalDeviate - Normal via Leva's Method

**Algorithm**: Leva's ratio-of-uniforms with quick rejection tests

**Class**:
```cpp
class NormalDeviate {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;  // Uniform(0, 1)
    Real _mu, _sigma;
    
public:
    NormalDeviate(Real mean, Real stddev, uint64_t seed = std::random_device{}());
    
    Real generate();
    Real operator()();
};
```

**Implementation** (optimized rejection with quick tests):
```cpp
Real generate() {
    Real u, v, x, y, q;
    do {
        u = _uniform(_gen);
        v = 1.7156 * (_uniform(_gen) - 0.5);
        x = u - 0.449871;
        y = std::abs(v) + 0.386595;
        q = x*x + y*(0.19600*y - 0.25472*x);
    } while (q > 0.27597 && 
            (q > 0.27846 || v*v > -4.0*std::log(u)*u*u));
    return _mu + _sigma * v/u;
}
```

**Performance**: 
- **Very fast** - optimized for single values
- **No caching** overhead
- Average **~1.37 iterations** to accept

**Quality**: **Exact** method

**Use When**:
- Need single normal values
- **Prefer this over Box-Muller** for most cases
- Maximum efficiency

---

### CauchyDeviate - Cauchy Random Numbers

**Algorithm**: Ratio of uniforms on unit circle
```
If V₁, V₂ uniform on unit circle, then X = V₁/V₂ ~ Cauchy(0,1)
```

**Class**:
```cpp
class CauchyDeviate {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;  // Uniform(-1, 1)
    Real _mu, _sigma;
    
public:
    CauchyDeviate(Real location, Real scale, uint64_t seed = std::random_device{}());
    
    Real generate();
    Real operator()();
};
```

**Implementation**:
```cpp
Real generate() {
    Real v1, v2;
    do {
        v1 = _uniform(_gen);
        v2 = _uniform(_gen);
    } while (v1*v1 + v2*v2 >= 1.0 || v2 == 0.0);
    return _mu + _sigma * v1/v2;
}
```

**Performance**: 
- Average **π/2 ≈ 1.57 iterations** to accept
- **Fast** overall

**Quality**: **Exact** method

---

### GammaDeviate - Gamma Random Numbers

**Algorithm**: Marsaglia and Tsang's method (2000) with transformation for α < 1

**Mathematical Background**:
Gamma distribution Γ(α, β):
```
f(x) = (1/(β^α·Γ(α))) · x^(α-1) · exp(-x/β)    for x > 0
```
- **α**: Shape parameter
- **β**: Scale parameter (β = 1/rate)
- **Mean**: α·β
- **Variance**: α·β²

**Class**:
```cpp
class GammaDeviate {
private:
    NormalDeviate _normalGen;
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;
    Real _alpha, _originalAlpha;
    Real _beta;
    Real _a1, _a2;  // Precomputed constants
    
public:
    GammaDeviate(Real shape, Real scale, uint64_t seed = std::random_device{}());
    
    Real generate();
    Real operator()();
};
```

**Algorithm Details**:

**For α ≥ 1** (Marsaglia-Tsang):
```
1. Set d = α - 1/3,  c = 1/√(9d)
2. Repeat:
   a) X ~ Normal(0,1)
   b) v = (1 + c·X)³
   c) If v > 0:
      U ~ Uniform(0,1)
      Accept if U < 1 - 0.331·X⁴
      Or if log(U) < 0.5·X² + d·(1 - v + log(v))
3. Return d·v/β
```

**For α < 1** (transformation):
```
1. Generate Y ~ Gamma(α+1, β) using method above
2. U ~ Uniform(0,1)
3. Return Y · U^(1/α)
```

**Performance**: 
- **Efficient** - typically 1-2 iterations
- Optimal for α > 1

**Quality**: **High** - exact for α ≥ 1, near-exact for α < 1

**Use When**:
- Need gamma/chi-squared/Erlang distributions
- Bayesian prior distributions
- Queuing theory simulations

---

### PoissonDeviate - Poisson Random Integers

**Mathematical Background**:
Poisson distribution Pois(λ):
```
P(X = k) = (λ^k · exp(-λ)) / k!    for k = 0, 1, 2, ...
```
- **λ**: Mean rate (λ > 0)
- **Mean**: λ
- **Variance**: λ

**Class**:
```cpp
class PoissonDeviate {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;
    Real _lambda;
    Real _sqrtLambda, _logLambda, _lambdaExp;
    Real _previousLambda;  // For caching
    std::vector<Real> _logFactorial;  // Cache log(k!)
    
public:
    PoissonDeviate(Real mean, uint64_t seed = std::random_device{}());
    
    int generate();
    int generate(Real mean);  // Override mean for this call
    int operator()();
};
```

**Adaptive Algorithm**:

**For λ < 5** (Direct method / Knuth):
```
1. Set L = exp(-λ)
2. k = 0, p = 1
3. Repeat:
   U ~ Uniform(0,1)
   p *= U
   If p > L: k++
   Else: return k
```

**For λ ≥ 5** (PTRS - Poisson-Transformed Rejection Sampling):
```
Complex ratio-of-uniforms method with quick acceptance tests
Optimized for large λ - see implementation for details
Average ~1.4 iterations to accept
```

**Performance**: 
- **λ < 5**: O(λ) expected iterations
- **λ ≥ 5**: **O(1)** expected time - very fast
- Automatic method selection

**Quality**: **Exact** discrete distribution

**Use When**:
- Counting rare events
- Queue arrivals modeling
- Radioactive decay counts
- Genomics (read depth)

---

### BinomialDeviate - Binomial Random Integers

**Mathematical Background**:
Binomial distribution B(n, p):
```
P(X = k) = C(n,k) · p^k · (1-p)^(n-k)    for k = 0, 1, ..., n
```
- **n**: Number of trials
- **p**: Success probability
- **Mean**: n·p
- **Variance**: n·p·(1-p)

**Class**:
```cpp
class BinomialDeviate {
private:
    std::mt19937_64 _gen;
    std::uniform_real_distribution<Real> _uniform;
    Real _pp, _p, _pb;
    Real _expNP, _np, _glnp, _plog, _pclog, _sq;
    int _n;
    int _method;  // 0, 1, or 2
    
    // Method-specific data structures
    uint64_t _uz, _uo, _unfin, _diff, _rltp;
    int _pbits[5];
    Real _cdf[64];
    Real _logFactorial[1024];
    
public:
    BinomialDeviate(int trials, Real probability, uint64_t seed = std::random_device{}());
    
    int generate();
    int operator()();
};
```

**Adaptive Three-Method Algorithm**:

**Method 0** (n ≤ 64): **Bit comparison**
```
Precompute binary representation of p
For each trial, compare random bits
Very fast for small n
```

**Method 1** (n > 64, n·p < 30): **Inverse CDF**
```
Precompute CDF array
Binary search to find k where CDF[k-1] < U ≤ CDF[k]
```

**Method 2** (n·p ≥ 30): **Ratio-of-uniforms (BTRS)**
```
Similar to Poisson PTRS method
Quick acceptance tests for efficiency
Handles large n·p efficiently
```

**Probability Flipping**: Works with min(p, 1-p) for efficiency, flips result if p > 0.5

**Performance**: 
- **Method 0**: O(n) but very fast for n ≤ 64
- **Method 1**: O(log n) lookup after O(n) precomputation
- **Method 2**: **O(1)** expected - fastest for large n·p

**Quality**: **Exact** discrete distribution

**Use When**:
- Clinical trials (successes in n patients)
- Quality control (defects in n items)
- Genetics (allele counts)
- A/B testing simulations

---

## Integration with Other Modules

### Vector Class
All descriptive statistics functions operate on `Vector<Real>`:
```cpp
Vector<Real> data = {1.2, 3.4, 2.1, 4.5, 3.8};
Real mean = Statistics::Avg(data);
```

### Function Analysis
Combine with `RealFunctionAnalyzer`:
```cpp
// Generate random samples, analyze empirically
NormalDeviate norm(0.0, 1.0);
Vector<Real> samples(10000);
for (int i = 0; i < 10000; i++)
    samples[i] = norm.generate();

Real avg, var;
Statistics::AvgVar(samples, avg, var);
// Should be close to (0, 1)
```

### Monte Carlo Integration
Use random deviates for MC integration:
```cpp
// Estimate ∫₀¹ f(x)dx using importance sampling
ExponentialDistribution proposal(1.0);
ExponentialDeviate sampler(1.0);

Real sum = 0.0;
for (int i = 0; i < N; i++) {
    Real x = sampler.generate();
    if (x <= 1.0)
        sum += f(x) / proposal.pdf(x);
}
Real estimate = sum / N;
```

### Advanced Statistical Distributions

#### Normal Distribution
Standard normal and arbitrary μ, σ:
```cpp
NormalDistribution stdNorm(0.0, 1.0);
Real z = 1.96;  // 97.5th percentile
Real p = stdNorm.cdf(z);  // 0.975
Real x = stdNorm.inverseCdf(0.95);  // 1.645
```

#### T-Distribution (Student's t)
For small-sample inference:
```cpp
TDistribution t(10);  // 10 degrees of freedom
Real t_critical = t.inverseCdf(0.975);  // Two-tailed 95% CI
Real p_value = 1.0 - t.cdf(2.5);  // One-tailed p-value
```

#### Chi-Square Distribution
For variance tests and goodness-of-fit:
```cpp
ChiSquareDistribution chi2(5);  // df = 5
Real chi2_critical = chi2.inverseCdf(0.95);  // 11.07
Real p = chi2.cdf(10.0);  // Cumulative probability
```

#### F-Distribution
For ANOVA and variance ratio tests:
```cpp
FDistribution f(3, 20);  // df1=3, df2=20
Real f_critical = f.inverseCdf(0.95);  // 3.10
Real p = 1.0 - f.cdf(4.5);  // P(F > 4.5)
```

## Hypothesis Testing

### T-Tests

#### One-Sample t-Test
Test if sample mean differs from hypothesized value:
```cpp
Vector<Real> data = {12.5, 13.1, 11.8, 12.9, 12.3};
TTestResult result = OneSampleTTest(data, 12.0);  // H₀: μ = 12.0

// Result contains:
// - testStatistic: t-value
// - pValue: two-tailed p-value
// - degreesOfFreedom: n-1
// - sampleMean, sampleStdDev
// - confidenceInterval95: (lower, upper)
```

#### Two-Sample t-Test
Compare means of two independent groups:
```cpp
Vector<Real> group1 = {23, 25, 27, 24, 26};
Vector<Real> group2 = {18, 20, 19, 21, 17};
TTestResult result = TwoSampleTTest(group1, group2);

// Tests H₀: μ₁ = μ₂
// Uses Welch's approximation (unequal variances)
```

#### Paired t-Test
Compare paired observations (before/after, matched pairs):
```cpp
Vector<Real> before = {120, 135, 128, 142, 138};
Vector<Real> after  = {115, 130, 125, 135, 132};
TTestResult result = PairedTTest(before, after);

// Tests H₀: μ_diff = 0
// More powerful than two-sample when observations paired
```

### Chi-Square Tests

#### Goodness-of-Fit Test
Test if observed frequencies match expected distribution:
```cpp
Vector<Real> observed = {25, 30, 20, 25};
Vector<Real> expected = {25, 25, 25, 25};  // Uniform expected
ChiSquareTestResult result = ChiSquareGoodnessOfFit(observed, expected);

// χ² = Σ[(O-E)²/E]
// Tests H₀: data follows expected distribution
```

#### Test of Independence
Test if two categorical variables are independent:
```cpp
Matrix<Real> contingencyTable(2, 2);
contingencyTable(0,0) = 30; contingencyTable(0,1) = 10;  // Success: Treatment A, B
contingencyTable(1,0) = 15; contingencyTable(1,1) = 25;  // Failure: Treatment A, B

ChiSquareTestResult result = ChiSquareTestOfIndependence(contingencyTable);

// Tests H₀: variables are independent
// df = (rows-1)(cols-1)
```

### ANOVA (Analysis of Variance)

#### One-Way ANOVA
Compare means across multiple groups:
```cpp
std::vector<Vector<Real>> groups = {
    {23, 25, 27, 24, 26},  // Group 1
    {18, 20, 19, 21, 17},  // Group 2
    {30, 32, 31, 33, 29}   // Group 3
};
ANOVAResult result = OneWayANOVA(groups);

// Result contains:
// - fStatistic: F-value
// - pValue: probability under H₀
// - dfBetween, dfWithin: degrees of freedom
// - meanSquareBetween, meanSquareWithin: variance components
// - grandMean: overall mean across all groups

// Tests H₀: μ₁ = μ₂ = μ₃ = ...
```

**ANOVA Assumptions**:
- Independence of observations
- Normal distribution within each group
- Homogeneity of variance (Levene's test recommended)

## Confidence Intervals

### CI for Population Mean
Estimate population mean from sample:
```cpp
Vector<Real> sample = {10, 12, 14, 16, 18};
ConfidenceInterval ci = ConfidenceIntervalMean(sample, 0.95);

// Returns:
// - estimate: sample mean
// - lowerBound, upperBound: 95% CI
// - marginOfError: half-width of interval
// - confidenceLevel: 0.95
// - parameter: "Mean"

// Interpretation: 95% confident true μ is in [lower, upper]
```

### CI for Mean Difference (Two Samples)
Compare two population means:
```cpp
Vector<Real> group1 = {23, 25, 27, 24, 26};
Vector<Real> group2 = {18, 20, 19, 21, 17};
ConfidenceInterval ci = ConfidenceIntervalMeanDifference(group1, group2, 0.95);

// estimate: x̄₁ - x̄₂
// If CI doesn't contain 0, means are significantly different
```

### CI for Proportion
Estimate population proportion (binomial):
```cpp
int successes = 65;
int trials = 100;
ConfidenceInterval ci = ConfidenceIntervalProportion(successes, trials, 0.95);

// estimate: p̂ = 65/100 = 0.65
// Uses normal approximation (valid when np̂ ≥ 5 and n(1-p̂) ≥ 5)
```

### CI for Proportion Difference
Compare two population proportions:
```cpp
int s1 = 45, n1 = 100;  // Treatment A: 45% success
int s2 = 30, n2 = 100;  // Treatment B: 30% success
ConfidenceInterval ci = ConfidenceIntervalProportionDifference(s1, n1, s2, n2, 0.95);

// estimate: p̂₁ - p̂₂ = 0.15
// If CI doesn't contain 0, proportions differ significantly
```

### CI for Paired Difference
Confidence interval for paired observations:
```cpp
Vector<Real> before = {120, 135, 128, 142, 138};
Vector<Real> after  = {115, 130, 125, 135, 132};
ConfidenceInterval ci = ConfidenceIntervalPairedDifference(before, after, 0.95);

// estimate: mean(before - after)
// More precise than independent samples when data paired
```

**CI Interpretation**:
- 95% CI: If we repeated sampling infinitely, 95% of intervals would contain true parameter
- Wider CI = more uncertainty (smaller n, larger σ)
- CI that excludes null value (e.g., 0 for differences) indicates statistical significance

---

## Examples

### Example 1: Basic Descriptive Statistics

Analyze experimental measurements:

```cpp
#include "algorithms/Statistics.h"

void Example1() {
    // Measurement data (lengths in mm)
    Vector<Real> measurements = {
        10.2, 10.1, 10.3, 10.0, 10.2, 10.4, 10.1, 10.2, 10.3, 10.1
    };
    
    Real avg, var;
    Statistics::AvgVar(measurements, avg, var);
    
    std::cout << "Sample size: " << measurements.size() << "\n";
    std::cout << "Mean: " << avg << " mm\n";
    std::cout << "Variance: " << var << " mm²\n";
    std::cout << "Std dev: " << std::sqrt(var) << " mm\n";
    
    // Output:
    // Mean: 10.19 mm
    // Variance: 0.01433 mm²
    // Std dev: 0.1197 mm
}
```

### Example 2: Complete Moment Analysis

Characterize distribution shape:

```cpp
void Example2() {
    // Asymmetric data (incomes in $1000s)
    Vector<Real> incomes = {
        45, 52, 48, 55, 62, 58, 51, 49, 53, 47,
        150, 180, 95  // High earners create right skew
    };
    
    Real ave, adev, sdev, var, skew, curt;
    Statistics::Moments(incomes, ave, adev, sdev, var, skew, curt);
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Mean: " << ave << "\n";
    std::cout << "Std dev: " << sdev << "\n";
    std::cout << "Skewness: " << skew << "\n";
    std::cout << "Kurtosis: " << curt << "\n";
    
    // Interpretation
    if (skew > 0.5)
        std::cout << "Distribution is right-skewed (long right tail)\n";
    if (curt > 0)
        std::cout << "Distribution has heavy tails (outlier-prone)\n";
    
    // Output:
    // Mean: 68.692
    // Std dev: 40.158
    // Skewness: 1.834  (strong right skew!)
    // Kurtosis: 2.456  (heavy tails!)
}
```

### Example 3: Distribution Properties - Exponential

Verify exponential distribution properties:

```cpp
void Example3() {
    Real lambda = 0.5;  // Rate = 0.5, Mean = 2.0
    Statistics::ExponentialDistribution exp(lambda);
    
    std::cout << "Exponential(λ=" << lambda << ")\n";
    std::cout << "Mean: " << exp.mean() << "\n";      // 2.0
    std::cout << "Variance: " << exp.variance() << "\n";  // 4.0
    std::cout << "Median: " << exp.inverseCdf(0.5) << "\n";  // ln(2)/λ
    
    // PDF at various points
    std::cout << "\nPDF values:\n";
    for (Real x = 0.0; x <= 6.0; x += 1.0) {
        std::cout << "f(" << x << ") = " << exp.pdf(x) << "\n";
    }
    
    // CDF - probability mass concentration
    std::cout << "\nP(X ≤ 2.0) = " << exp.cdf(2.0) << "\n";  // 0.632
    std::cout << "P(X ≤ 4.0) = " << exp.cdf(4.0) << "\n";  // 0.865
    
    // Quantiles
    std::cout << "\n90th percentile: " << exp.inverseCdf(0.9) << "\n";
}
```

### Example 4: Logistic vs Normal Comparison

Compare logistic and normal distributions:

```cpp
void Example4() {
    Statistics::LogisticDistribution logistic(0.0, 1.0);
    
    // Equivalent normal: σ_normal = σ_logistic · π/√3
    Real sigma_normal = 1.0 * Constants::PI / std::sqrt(3.0);
    
    std::cout << "Logistic(0, 1) vs Normal(0, " << sigma_normal << ")\n\n";
    
    // Compare PDFs
    std::cout << "x\tLogistic PDF\tNormal PDF (approx)\n";
    for (Real x = -3.0; x <= 3.0; x += 1.0) {
        Real pdf_logistic = logistic.pdf(x);
        // Approximation using formula
        Real pdf_normal_approx = std::exp(-x*x/(2*sigma_normal*sigma_normal)) 
                                  / (sigma_normal * std::sqrt(2*Constants::PI));
        std::cout << x << "\t" << pdf_logistic << "\t" << pdf_normal_approx << "\n";
    }
    
    // Compare tail behavior
    std::cout << "\nTail comparison P(X > x):\n";
    std::cout << "x\tLogistic\tNormal\n";
    for (Real x = 1.0; x <= 5.0; x += 1.0) {
        Real tail_logistic = 1.0 - logistic.cdf(x);
        std::cout << x << "\t" << tail_logistic << "\t(Logistic heavier tails)\n";
    }
}
```

### Example 5: Random Sampling - Normal Distribution

Generate and verify normal samples:

```cpp
void Example5() {
    Real mu = 100.0, sigma = 15.0;
    Statistics::NormalDeviate normGen(mu, sigma, 42);  // Seed = 42
    
    // Generate 10,000 samples
    Vector<Real> samples(10000);
    for (int i = 0; i < 10000; i++)
        samples[i] = normGen.generate();
    
    // Analyze empirical distribution
    Real emp_mean, emp_var;
    Statistics::AvgVar(samples, emp_mean, emp_var);
    Real emp_sigma = std::sqrt(emp_var);
    
    std::cout << "Target:    μ = " << mu << ", σ = " << sigma << "\n";
    std::cout << "Empirical: μ = " << emp_mean << ", σ = " << emp_sigma << "\n";
    std::cout << "Error: μ = " << std::abs(emp_mean - mu) 
              << ", σ = " << std::abs(emp_sigma - sigma) << "\n";
    
    // Should be very close with n=10,000
    // Typical error: μ within 0.5, σ within 0.3
    
    // Check normality via moments
    Real ave, adev, sdev, var, skew, curt;
    Statistics::Moments(samples, ave, adev, sdev, var, skew, curt);
    
    std::cout << "\nShape parameters:\n";
    std::cout << "Skewness: " << skew << " (should be ~0)\n";
    std::cout << "Kurtosis: " << curt << " (should be ~0)\n";
}
```

### Example 6: Poisson Process Simulation

Simulate event counts in fixed intervals:

```cpp
void Example6() {
    Real lambda = 5.0;  // Average 5 events per interval
    Statistics::PoissonDeviate poisson(lambda, 123);
    
    // Simulate 100 intervals
    std::vector<int> counts(100);
    for (int i = 0; i < 100; i++)
        counts[i] = poisson.generate();
    
    // Analyze distribution of counts
    Vector<Real> countData(100);
    for (int i = 0; i < 100; i++)
        countData[i] = static_cast<Real>(counts[i]);
    
    Real avg, var;
    Statistics::AvgVar(countData, avg, var);
    
    std::cout << "Poisson(λ=" << lambda << ") simulation:\n";
    std::cout << "Empirical mean: " << avg << " (theoretical: " << lambda << ")\n";
    std::cout << "Empirical var:  " << var << " (theoretical: " << lambda << ")\n";
    
    // Count frequencies
    std::map<int, int> freq;
    for (int count : counts)
        freq[count]++;
    
    std::cout << "\nFrequency distribution:\n";
    std::cout << "k\tCount\tFreq\n";
    for (const auto& [k, n] : freq) {
        std::cout << k << "\t" << n << "\t" << (n/100.0) << "\n";
    }
}
```

### Example 7: Binomial Confidence Interval

Estimate proportion with confidence interval:

```cpp
void Example7() {
    int n = 100;      // Sample size
    Real p = 0.3;     // True proportion
    
    Statistics::BinomialDeviate binomial(n, p, 999);
    
    // Run 1000 experiments
    Vector<Real> proportions(1000);
    for (int i = 0; i < 1000; i++) {
        int successes = binomial.generate();
        proportions[i] = static_cast<Real>(successes) / n;
    }
    
    Real avg, var;
    Statistics::AvgVar(proportions, avg, var);
    Real se = std::sqrt(var);  // Standard error
    
    std::cout << "Binomial(" << n << ", " << p << ") proportion estimation:\n";
    std::cout << "Empirical mean: " << avg << " (true p: " << p << ")\n";
    std::cout << "Standard error: " << se << "\n";
    std::cout << "95% CI (approx): [" << (avg - 1.96*se) << ", " 
              << (avg + 1.96*se) << "]\n";
    
    // Theoretical SE = √(p(1-p)/n)
    Real theoretical_se = std::sqrt(p * (1-p) / n);
    std::cout << "Theoretical SE: " << theoretical_se << "\n";
}
```

### Example 8: Gamma Distribution - Waiting Times

Model waiting time for k events in Poisson process:

```cpp
void Example8() {
    // Gamma(k, θ) = waiting time for k events at rate λ = 1/θ
    Real k = 5.0;     // Shape (number of events)
    Real theta = 2.0; // Scale (1/rate)
    
    Statistics::GammaDeviate gamma(k, theta, 555);
    
    // Generate 5000 waiting times
    Vector<Real> waitingTimes(5000);
    for (int i = 0; i < 5000; i++)
        waitingTimes[i] = gamma.generate();
    
    Real avg, var;
    Statistics::AvgVar(waitingTimes, avg, var);
    
    std::cout << "Gamma(" << k << ", " << theta << ") simulation:\n";
    std::cout << "Theoretical mean: " << (k * theta) << "\n";
    std::cout << "Empirical mean:   " << avg << "\n";
    std::cout << "Theoretical var:  " << (k * theta * theta) << "\n";
    std::cout << "Empirical var:    " << var << "\n";
    
    // Shape analysis
    Real ave, adev, sdev, v, skew, curt;
    Statistics::Moments(waitingTimes, ave, adev, sdev, v, skew, curt);
    
    Real theoretical_skew = 2.0 / std::sqrt(k);  // = 2/√5 ≈ 0.894
    std::cout << "\nSkewness:\n";
    std::cout << "Theoretical: " << theoretical_skew << "\n";
    std::cout << "Empirical:   " << skew << "\n";
}
```

### Example 9: Cauchy Heavy Tails - Outlier Demo

Demonstrate Cauchy distribution's undefined mean:

```cpp
void Example9() {
    Statistics::CauchyDeviate cauchy(0.0, 1.0, 777);
    
    std::cout << "Cauchy(0, 1) - Sample mean convergence:\n";
    std::cout << "n\tSample Mean\t|Mean| > 10?\n";
    
    Vector<Real> samples(100000);
    for (int i = 0; i < 100000; i++)
        samples[i] = cauchy.generate();
    
    // Check convergence at different sample sizes
    std::vector<int> sizes = {100, 1000, 10000, 100000};
    for (int n : sizes) {
        Vector<Real> subset(n);
        for (int i = 0; i < n; i++)
            subset[i] = samples[i];
        
        Real avg = Statistics::Avg(subset);
        bool unstable = (std::abs(avg) > 10.0);
        
        std::cout << n << "\t" << avg << "\t" 
                  << (unstable ? "YES (unstable!)" : "no") << "\n";
    }
    
    std::cout << "\nCauchy mean does NOT converge with n!\n";
    std::cout << "Contrast with normal: mean error ~ 1/√n\n";
}
```

### Example 10: Multi-Distribution Comparison

Compare tail behavior across distributions:

```cpp
void Example10() {
    Statistics::CauchyDistribution cauchy(0.0, 1.0);
    Statistics::LogisticDistribution logistic(0.0, 1.0);
    
    // Compare P(X > x) for various thresholds
    std::cout << "Tail probabilities P(X > x):\n";
    std::cout << "x\tCauchy\t\tLogistic\tRatio\n";
    
    for (Real x = 1.0; x <= 10.0; x += 2.0) {
        Real tail_cauchy = 1.0 - cauchy.cdf(x);
        Real tail_logistic = 1.0 - logistic.cdf(x);
        Real ratio = tail_cauchy / tail_logistic;
        
        std::cout << std::fixed << std::setprecision(6);
        std::cout << x << "\t" << tail_cauchy << "\t" 
                  << tail_logistic << "\t" << ratio << "\n";
    }
    
    std::cout << "\nCauchy tails decay as 1/x (power law)\n";
    std::cout << "Logistic tails decay exponentially\n";
    std::cout << "Ratio grows unbounded as x → ∞\n";
}
```

---

## Best Practices

### Choosing Descriptive Statistics

**Decision Tree**:
```
Need central tendency only?
├─ Symmetric data → Avg()
└─ Skewed/outliers → Use median (external computation)

Need spread + center?
├─ Simple case → AvgStdDev()
└─ Need variance explicitly → AvgVar()

Need complete shape characterization?
└─ Moments() for skewness, kurtosis
```

### Distribution Selection

| Application | Recommended Distribution | Why |
|-------------|--------------------------|-----|
| **Time between events** | Exponential | Memoryless property |
| **Resonance/spectral** | Cauchy (Lorentzian) | Physical line shape |
| **Neural networks** | Logistic | Sigmoid activation |
| **Count data** | Poisson | Discrete, unbounded |
| **Success counts** | Binomial | Fixed trials |
| **General continuous** | Normal | Central limit theorem |
| **Waiting for k events** | Gamma | Sum of exponentials |

### Random Number Generator Selection

**Speed Priority**:
1. ExponentialDeviate (fastest - one log)
2. NormalDeviate (Leva - optimized)
3. LogisticDeviate (one log)
4. GammaDeviate, PoissonDeviate (adaptive)

**Theoretical Exactness**:
- All inverse transform methods: **Exact**
- Box-Muller: **Exact**
- Rejection samplers (Gamma, Poisson, Binomial): **Exact** (within FP precision)

**Thread Safety**:
- Each generator instance has **independent state**
- ✅ Safe if each thread has own instance
- ❌ NOT safe if shared across threads without locking
- **Recommendation**: Thread-local instances or lock-free queue of generators

### Seeding Strategy

**Reproducible Results** (testing, debugging):
```cpp
Statistics::NormalDeviate gen(0.0, 1.0, 42);  // Fixed seed
```

**Non-Reproducible** (production simulations):
```cpp
Statistics::NormalDeviate gen(0.0, 1.0);  // std::random_device{}
```

**Multiple Streams** (parallel simulations):
```cpp
std::vector<Statistics::NormalDeviate> generators;
for (int i = 0; i < numThreads; i++)
    generators.emplace_back(0.0, 1.0, baseSeed + i);
```

### Common Pitfalls

❌ **Pitfall 1**: Assuming Cauchy has finite mean
```cpp
// DON'T:
CauchyDeviate cauchy(0.0, 1.0);
Vector<Real> samples(10000);
for (int i = 0; i < 10000; i++)
    samples[i] = cauchy.generate();
Real mean = Statistics::Avg(samples);  // Meaningless! Won't converge

// DO: Use median or mode (both = μ)
```

❌ **Pitfall 2**: Using wrong variance formula
```cpp
// DON'T (biased):
Real variance = sumSquaredDiffs / n;  // Biased!

// DO (unbiased):
Real variance = sumSquaredDiffs / (n - 1);  // Bessel's correction
```

❌ **Pitfall 3**: Ignoring numerical stability in moments
```cpp
// DON'T (catastrophic cancellation):
Real var = (sumX2 / n) - (sumX / n)^2;  // Numerically unstable!

// DO (two-pass with correction):
Real mean = sumX / n;
Real var = (sum((x - mean)^2) - correction) / (n - 1);
```

❌ **Pitfall 4**: Reusing same seed in parallel
```cpp
// DON'T:
#pragma omp parallel for
for (int i = 0; i < N; i++) {
    NormalDeviate gen(0, 1, 42);  // All threads get SAME sequence!
    results[i] = gen.generate();
}

// DO:
std::vector<NormalDeviate> gens(omp_get_max_threads());
for (int i = 0; i < gens.size(); i++)
    gens[i] = NormalDeviate(0, 1, 42 + i);
    
#pragma omp parallel for
for (int i = 0; i < N; i++) {
    int tid = omp_get_thread_num();
    results[i] = gens[tid].generate();
}
```

✅ **Best Practice Pattern** (Statistical Analysis):
```cpp
// 1. Collect data
Vector<Real> data = CollectMeasurements();

// 2. Check for outliers/validity
if (data.size() < 3)
    throw std::runtime_error("Insufficient data");

// 3. Compute descriptive statistics
Real ave, adev, sdev, var, skew, curt;
Statistics::Moments(data, ave, adev, sdev, var, skew, curt);

// 4. Assess distribution shape
bool isSymmetric = (std::abs(skew) < 0.5);
bool hasOutliers = (std::abs(curt) > 1.0);

// 5. Choose appropriate analysis
if (isSymmetric && !hasOutliers) {
    // Normal approximation valid
    // Use mean ± 1.96·σ for 95% CI
} else {
    // Use robust methods (median, IQR)
    // Or fit appropriate distribution
}
```

---

## Performance Considerations

### Complexity Summary

| Operation | Time Complexity | Space Complexity |
|-----------|-----------------|------------------|
| **Avg** | O(n) | O(1) |
| **AvgVar** | O(n) (two-pass) | O(1) |
| **Moments** | O(n) (two-pass) | O(1) |
| **Generate exponential** | O(1) expected | O(1) |
| **Generate normal** | O(1) expected | O(1) |
| **Generate Poisson (λ<5)** | O(λ) expected | O(1) |
| **Generate Poisson (λ≥5)** | O(1) expected | O(λ) precompute |
| **Generate binomial** | O(1) to O(n) | O(min(n, 64)) |

### Optimization Tips

**Batch Processing**:
```cpp
// Generate many samples efficiently
NormalDeviate gen(0, 1);
Vector<Real> samples(1000000);
samples.Apply([&gen](Real x) { return gen.generate(); });
```

**Precomputation** (for repeated distribution queries):
```cpp
// If querying same distribution many times
ExponentialDistribution exp(0.5);
// Precompute common quantiles
std::vector<Real> percentiles(100);
for (int i = 1; i < 100; i++)
    percentiles[i] = exp.inverseCdf(i / 100.0);
```

**Cache Locality**:
```cpp
// Better: Process data in chunks
const int CHUNK_SIZE = 1024;
for (int i = 0; i < n; i += CHUNK_SIZE) {
    // Process chunk...
}
```

---

## Advanced Topics

### Copulas (Future Extension)

**Concept**: Model dependence structure separately from marginals
```cpp
// Pseudo-code for Gaussian copula
// Generate correlated uniforms → Transform to desired marginals
```

### Mixture Distributions

**Example**: Gaussian mixture model
```cpp
// Pseudo-implementation
class GaussianMixture {
    std::vector<NormalDeviate> components;
    std::vector<Real> weights;
    
    Real generate() {
        int k = SampleCategorical(weights);
        return components[k].generate();
    }
};
```

### Kernel Density Estimation

**Use statistics to build empirical distributions**:
```cpp
// From sample data, estimate PDF
class KernelDensity {
    Vector<Real> data;
    Real bandwidth;
    
    Real pdf(Real x) {
        Real sum = 0;
        for (Real xi : data)
            sum += GaussianKernel((x - xi) / bandwidth);
        return sum / (data.size() * bandwidth);
    }
};
```

---

## Summary

### Key Features

**Descriptive Statistics**:
- ✅ Mean, variance, standard deviation (with Bessel's correction)
- ✅ Covariance and correlation (Pearson's r)
- ✅ All four moments (mean, variance, skewness, kurtosis)
- ✅ Numerically stable two-pass algorithms

**Distributions** (7 types):
- ✅ Normal: Standard inference, z-tests
- ✅ T-Distribution: Small-sample inference
- ✅ Chi-Square: Variance tests, goodness-of-fit
- ✅ F-Distribution: ANOVA, variance ratios
- ✅ Cauchy: Heavy tails, no mean/variance
- ✅ Exponential: Memoryless, waiting times
- ✅ Logistic: Sigmoid CDF, ML applications
- ✅ All with pdf, cdf, inverseCdf (quantile)

**Hypothesis Tests** (6 types):
- ✅ One-sample t-test (μ = μ₀)
- ✅ Two-sample t-test (μ₁ = μ₂, Welch's method)
- ✅ Paired t-test (μ_diff = 0)
- ✅ Chi-square goodness-of-fit test
- ✅ Chi-square test of independence
- ✅ One-way ANOVA (multiple group comparison)

**Confidence Intervals** (5 types):
- ✅ Single population mean
- ✅ Difference of two means
- ✅ Single proportion (binomial)
- ✅ Difference of two proportions
- ✅ Paired difference mean

**Random Generators** (7 types):
- ✅ Mersenne Twister 64-bit engine (period 2¹⁹⁹³⁷−1)
- ✅ Exponential, Logistic, Cauchy (exact inverse transform)
- ✅ Normal: Box-Muller + Leva's ratio-of-uniforms
- ✅ Gamma: Marsaglia-Tsang (shape ≥ 1) + transformation (< 1)
- ✅ Poisson: Adaptive (Knuth direct + PTRS rejection)
- ✅ Binomial: Three-method adaptive (bits, CDF, BTRS)

### When to Use This Module

**Descriptive Statistics**:
- Experimental data analysis
- Quality control monitoring
- Distribution shape assessment
- Correlation analysis between variables

**Distributions**:
- Theoretical probability calculations
- Quantile/percentile lookups
- Critical value determination
- Risk analysis (tail probabilities)

**Hypothesis Tests**:
- A/B testing (t-tests, proportions)
- Quality control (goodness-of-fit)
- Categorical data analysis (independence)
- Multi-group comparison (ANOVA)

**Confidence Intervals**:
- Parameter estimation with uncertainty
- Effect size quantification
- Clinical trials and experiments
- Survey analysis (proportions)

**Random Generators**:
- Monte Carlo simulations
- Stochastic differential equations
- Agent-based modeling
- Bootstrap/resampling methods
- Synthetic data generation

### Statistical Workflow

1. **Explore**: Use descriptive statistics (Mean, StdDev, Correlation)
2. **Visualize**: Check distribution shape (skewness, kurtosis)
3. **Test**: Apply hypothesis tests (t-tests, ANOVA, chi-square)
4. **Estimate**: Compute confidence intervals for parameters
5. **Simulate**: Generate random data for validation/power analysis

### References

- **Probability**: Feller, "An Introduction to Probability Theory and Its Applications"
- **Hypothesis Testing**: Hogg & Tanis, "Probability and Statistical Inference"
- **ANOVA**: Montgomery, "Design and Analysis of Experiments"
- **Confidence Intervals**: Casella & Berger, "Statistical Inference"
- **Random Number Generation**: Press et al., "Numerical Recipes" Chapter 7
- **Statistical Moments**: Kendall & Stuart, "The Advanced Theory of Statistics"
- **Modern RNG Algorithms**: Marsaglia & Tsang (2000), Ahrens & Dieter (1982)
- **Mersenne Twister**: Matsumoto & Nishimura (1998)

---

**Part of MinimalMathLibrary** - `mml/algorithms/Statistics.h`
