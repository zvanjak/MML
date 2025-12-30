#ifndef MML_TEST_MATCHERS_H
#define MML_TEST_MATCHERS_H

#include "TestPrecision.h"
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace MML {
namespace Testing {
namespace Matchers {

// ============================================================================
// RealEquals - Exact equality matcher for Real values
// ============================================================================

class RealEqualsMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_expected;
    
public:
    explicit RealEqualsMatcher(Real expected) : m_expected(expected) {}
    
    bool match(Real const& actual) const override {
        // Handle NaN specially
        if (std::isnan(m_expected) && std::isnan(actual)) return true;
        if (std::isnan(m_expected) || std::isnan(actual)) return false;
        
        // Exact equality
        return actual == m_expected;
    }
    
    std::string describe() const override {
        std::ostringstream oss;
        oss << "equals " << std::setprecision(GetRealDigits() + 2) << m_expected
            << " (Real=" << GetRealPrecisionName() << ")";
        return oss.str();
    }
};

inline RealEqualsMatcher RealEquals(Real expected) {
    return RealEqualsMatcher(expected);
}

// ============================================================================
// RealApprox - Approximate equality with automatic tolerance scaling
// ============================================================================

class RealApproxMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_expected;
    Real m_epsilon;
    Real m_margin;
    
public:
    explicit RealApproxMatcher(Real expected)
        : m_expected(expected)
        , m_epsilon(Tolerance::Standard)
        , m_margin(Tolerance::Standard) {}
    
    RealApproxMatcher& epsilon(Real eps) {
        m_epsilon = eps;
        return *this;
    }
    
    RealApproxMatcher& margin(Real m) {
        m_margin = m;
        return *this;
    }
    
    bool match(Real const& actual) const override {
        return ApproxEqual(actual, m_expected, m_epsilon, m_margin);
    }
    
    std::string describe() const override {
        std::ostringstream oss;
        oss << "is approx " << std::setprecision(GetRealDigits() + 2) << m_expected;
        oss << " (epsilon=" << m_epsilon << ", margin=" << m_margin;
        oss << ", Real=" << GetRealPrecisionName() << ")";
        return oss.str();
    }
};

inline RealApproxMatcher RealApprox(Real expected) {
    return RealApproxMatcher(expected);
}

// ============================================================================
// RealWithinRel - Relative tolerance comparison with auto-scaling
// ============================================================================

class RealWithinRelMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_expected;
    Real m_tolerance;
    bool m_use_scaled_tolerance;
    
public:
    RealWithinRelMatcher(Real expected, Real tolerance, bool auto_scale = true)
        : m_expected(expected)
        , m_tolerance(auto_scale ? ScaleRelTolerance(tolerance) : tolerance)
        , m_use_scaled_tolerance(auto_scale) {}
    
    bool match(Real const& actual) const override {
        // Handle exact equality (including infinities)
        if (actual == m_expected) return true;
        
        // Handle NaN
        if (std::isnan(actual) || std::isnan(m_expected)) return false;
        
        // Relative error check
        Real diff = std::abs(actual - m_expected);
        Real max_val = std::max(std::abs(actual), std::abs(m_expected));
        
        // Avoid division by zero
        if (max_val < std::numeric_limits<Real>::min()) {
            return diff <= m_tolerance;
        }
        
        return (diff / max_val) <= m_tolerance;
    }
    
    std::string describe() const override {
        std::ostringstream oss;
        oss << "is within " << (m_tolerance * 100.0) << "% of "
            << std::setprecision(GetRealDigits() + 2) << m_expected;
        if (m_use_scaled_tolerance) {
            oss << " (tolerance auto-scaled for " << GetRealPrecisionName() << ")";
        }
        return oss.str();
    }
};

inline RealWithinRelMatcher RealWithinRel(Real expected, Real tolerance = Tolerance::Standard) {
    return RealWithinRelMatcher(expected, tolerance, true);
}

inline RealWithinRelMatcher RealWithinRelNoScale(Real expected, Real tolerance) {
    return RealWithinRelMatcher(expected, tolerance, false);
}

// ============================================================================
// RealWithinAbs - Absolute tolerance comparison with auto-scaling
// ============================================================================

class RealWithinAbsMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_expected;
    Real m_tolerance;
    bool m_use_scaled_tolerance;
    
public:
    RealWithinAbsMatcher(Real expected, Real tolerance, bool auto_scale = true)
        : m_expected(expected)
        , m_tolerance(auto_scale ? ScaleTolerance(tolerance) : tolerance)
        , m_use_scaled_tolerance(auto_scale) {}
    
    bool match(Real const& actual) const override {
        // Handle exact equality
        if (actual == m_expected) return true;
        
        // Handle NaN
        if (std::isnan(actual) || std::isnan(m_expected)) return false;
        
        // Absolute error check
        return std::abs(actual - m_expected) <= m_tolerance;
    }
    
    std::string describe() const override {
        std::ostringstream oss;
        oss << "is within " << m_tolerance << " of "
            << std::setprecision(GetRealDigits() + 2) << m_expected;
        if (m_use_scaled_tolerance) {
            oss << " (tolerance auto-scaled for " << GetRealPrecisionName() << ")";
        }
        return oss.str();
    }
};

inline RealWithinAbsMatcher RealWithinAbs(Real expected, Real tolerance = Tolerance::Standard) {
    return RealWithinAbsMatcher(expected, tolerance, true);
}

inline RealWithinAbsMatcher RealWithinAbsNoScale(Real expected, Real tolerance) {
    return RealWithinAbsMatcher(expected, tolerance, false);
}

// ============================================================================
// RealIsZero - Check if value is approximately zero
// ============================================================================

class RealIsZeroMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_tolerance;
    
public:
    explicit RealIsZeroMatcher(Real tolerance = Tolerance::Standard)
        : m_tolerance(ScaleTolerance(tolerance)) {}
    
    bool match(Real const& actual) const override {
        return ApproxZero(actual, m_tolerance);
    }
    
    std::string describe() const override {
        std::ostringstream oss;
        oss << "is approximately zero (tolerance=" << m_tolerance
            << ", Real=" << GetRealPrecisionName() << ")";
        return oss.str();
    }
};

inline RealIsZeroMatcher RealIsZero(Real tolerance = Tolerance::Standard) {
    return RealIsZeroMatcher(tolerance);
}

// ============================================================================
// RealIsPositive / RealIsNegative - Sign checks
// ============================================================================

class RealIsPositiveMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_tolerance;
    
public:
    explicit RealIsPositiveMatcher(Real tolerance = REAL(0.0))
        : m_tolerance(tolerance) {}
    
    bool match(Real const& actual) const override {
        return actual > m_tolerance;
    }
    
    std::string describe() const override {
        if (m_tolerance == REAL(0.0)) {
            return "is positive";
        } else {
            std::ostringstream oss;
            oss << "is positive (> " << m_tolerance << ")";
            return oss.str();
        }
    }
};

class RealIsNegativeMatcher : public Catch::Matchers::MatcherBase<Real> {
    Real m_tolerance;
    
public:
    explicit RealIsNegativeMatcher(Real tolerance = REAL(0.0))
        : m_tolerance(tolerance) {}
    
    bool match(Real const& actual) const override {
        return actual < -m_tolerance;
    }
    
    std::string describe() const override {
        if (m_tolerance == REAL(0.0)) {
            return "is negative";
        } else {
            std::ostringstream oss;
            oss << "is negative (< " << -m_tolerance << ")";
            return oss.str();
        }
    }
};

inline RealIsPositiveMatcher RealIsPositive(Real tolerance = REAL(0.0)) {
    return RealIsPositiveMatcher(tolerance);
}

inline RealIsNegativeMatcher RealIsNegative(Real tolerance = REAL(0.0)) {
    return RealIsNegativeMatcher(tolerance);
}

} // namespace Matchers
} // namespace Testing
} // namespace MML

// ============================================================================
// Convenience using declarations for test files
// ============================================================================

// Import matchers into global namespace for convenience in tests
using MML::Testing::Matchers::RealEquals;
using MML::Testing::Matchers::RealApprox;
using MML::Testing::Matchers::RealWithinRel;
using MML::Testing::Matchers::RealWithinRelNoScale;
using MML::Testing::Matchers::RealWithinAbs;
using MML::Testing::Matchers::RealWithinAbsNoScale;
using MML::Testing::Matchers::RealIsZero;
using MML::Testing::Matchers::RealIsPositive;
using MML::Testing::Matchers::RealIsNegative;

#endif // MML_TEST_MATCHERS_H
