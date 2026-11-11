///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ForwardAD.h                                                         ///
///  Description: Forward-mode Automatic Differentiation using Dual Numbers           ///
///               Computes exact derivatives via operator overloading                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_AD_FORWARD_AD_H
#define MML_AD_FORWARD_AD_H

#include <mml/MMLBase.h>
#include <mml/core/Derivation/ADArithmetic.h>
#include <algorithm>

#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>

namespace MML::AD
{
    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Dual number for forward-mode automatic differentiation
    /// 
    /// A dual number extends real numbers with an infinitesimal part ε where ε² = 0.
    /// A dual number is written as: a + b*ε
    /// 
    /// This property enables exact derivative computation:
    /// f(a + ε) = f(a) + f'(a)*ε
    /// 
    /// Usage:
    /// @code
    ///   Dual x(2.0, 1.0);  // x = 2, dx/dx = 1
    ///   Dual y = sin(x*x); // y.value = sin(4), y.deriv = cos(4)*2*2
    /// @endcode
    /// 
    /// Forward AD is efficient for functions f: R → R^m (one pass per input variable).
    /// For f: R^n → R with large n, consider reverse-mode AD instead.
    //////////////////////////////////////////////////////////////////////////////////////////
    template<typename T = Real>
    struct Dual
    {
        static_assert(Detail::HasFloatingLeaf<T>::value, "Dual requires a floating-point leaf type");

        T value;  ///< Function value f(x)
        T deriv;  ///< Derivative f'(x)

        /// Default constructor: zero dual number
        Dual() : value(T(0)), deriv(T(0)) {}

        /// Construct from value only (constant, derivative = 0)
        explicit Dual(T val) : value(val), deriv(T(0)) {}

        /// Construct with both value and derivative
        Dual(T val, T d) : value(val), deriv(d) {}

        /// Implicit conversion from scalar (treated as constant)
        template<typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
        Dual(U val) : value(static_cast<T>(val)), deriv(T(0)) {}

        //////////////////////////////////////////////////////////////////////////////////
        // Compound assignment operators
        //////////////////////////////////////////////////////////////////////////////////

        Dual& operator+=(const Dual& other) {
            value += other.value;
            deriv += other.deriv;
            return *this;
        }

        Dual& operator-=(const Dual& other) {
            value -= other.value;
            deriv -= other.deriv;
            return *this;
        }

        Dual& operator*=(const Dual& other) {
            // (a + a'ε)(b + b'ε) = ab + (a'b + ab')ε
            deriv = deriv * other.value + value * other.deriv;
            value *= other.value;
            return *this;
        }

        Dual& operator/=(const Dual& other) {
            const T quotient = value / other.value;
            const T tangent = Detail::QuotientTangent(value, deriv, other.value, other.deriv);
            value = quotient;
            deriv = tangent;
            return *this;
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Arithmetic operators
    //////////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    Dual<T> operator+(const Dual<T>& a, const Dual<T>& b) {
        return Dual<T>(a.value + b.value, a.deriv + b.deriv);
    }

    template<typename T>
    Dual<T> operator-(const Dual<T>& a, const Dual<T>& b) {
        return Dual<T>(a.value - b.value, a.deriv - b.deriv);
    }

    template<typename T>
    Dual<T> operator*(const Dual<T>& a, const Dual<T>& b) {
        // Product rule: (a,a') * (b,b') = (ab, a'b + ab')
        return Dual<T>(a.value * b.value, a.deriv * b.value + a.value * b.deriv);
    }

    template<typename T>
    Dual<T> operator/(const Dual<T>& a, const Dual<T>& b) {
        return Dual<T>(a.value / b.value, Detail::QuotientTangent(a.value, a.deriv, b.value, b.deriv));
    }

    template<typename T>
    Dual<T> operator-(const Dual<T>& a) {
        return Dual<T>(-a.value, -a.deriv);
    }

    template<typename T>
    Dual<T> operator+(const Dual<T>& a) {
        return a;
    }

    // Mixed operations with scalars
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator+(const Dual<T>& a, S b) { return a + Dual<T>(static_cast<T>(b)); }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator+(S a, const Dual<T>& b) { return Dual<T>(static_cast<T>(a)) + b; }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator-(const Dual<T>& a, S b) { return a - Dual<T>(static_cast<T>(b)); }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator-(S a, const Dual<T>& b) { return Dual<T>(static_cast<T>(a)) - b; }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator*(const Dual<T>& a, S b) { return a * Dual<T>(static_cast<T>(b)); }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator*(S a, const Dual<T>& b) { return Dual<T>(static_cast<T>(a)) * b; }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator/(const Dual<T>& a, S b) { return a / Dual<T>(static_cast<T>(b)); }
    
    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> operator/(S a, const Dual<T>& b) { return Dual<T>(static_cast<T>(a)) / b; }

    //////////////////////////////////////////////////////////////////////////////////////////
    // Comparison operators (compare values only)
    //////////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool operator==(const Dual<T>& a, const Dual<T>& b) { return a.value == b.value; }
    
    template<typename T>
    bool operator!=(const Dual<T>& a, const Dual<T>& b) { return a.value != b.value; }
    
    template<typename T>
    bool operator<(const Dual<T>& a, const Dual<T>& b) { return a.value < b.value; }
    
    template<typename T>
    bool operator<=(const Dual<T>& a, const Dual<T>& b) { return a.value <= b.value; }
    
    template<typename T>
    bool operator>(const Dual<T>& a, const Dual<T>& b) { return a.value > b.value; }
    
    template<typename T>
    bool operator>=(const Dual<T>& a, const Dual<T>& b) { return a.value >= b.value; }

    //////////////////////////////////////////////////////////////////////////////////////////
    // Mathematical functions
    //////////////////////////////////////////////////////////////////////////////////////////

		/// Square root: d/dx(sqrt(x)) = 1/(2*sqrt(x))
    template<typename T>
    Dual<T> sqrt(const Dual<T>& x) {
			using std::sqrt;
			T sq = sqrt(x.value);
			return Dual<T>(sq, x.deriv / (T(2) * sq));
    }

    /// Power function: d/dx(x^n) = n*x^(n-1)
    template<typename T>
    Dual<T> pow(const Dual<T>& base, const Dual<T>& exp) {
			Detail::RequirePowerDomain(Detail::Primal(base), Detail::Primal(exp), Detail::IsConstantJet(exp));
			return Detail::PowJetUnchecked(base, exp);
		}

    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> pow(const Dual<T>& base, S exp) {
			return pow(base, Dual<T>(T(exp)));
		}

    template<typename T, typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
    Dual<T> pow(S base, const Dual<T>& exp) {
				return pow(Dual<T>(T(base)), exp);
		}

    /// Exponential: d/dx(e^x) = e^x
    template<typename T>
    Dual<T> exp(const Dual<T>& x) {
			using std::exp;
			T ex = exp(x.value);
			return Dual<T>(ex, ex * x.deriv);
    }

    /// Natural logarithm: d/dx(ln(x)) = 1/x
    template<typename T>
    Dual<T> log(const Dual<T>& x) {
			using std::log;
			return Dual<T>(log(x.value), x.deriv / x.value);
		}

    /// Base-10 logarithm: d/dx(log10(x)) = 1/(x*ln(10))
    template<typename T>
    Dual<T> log10(const Dual<T>& x) {
			using std::log;
			using std::log10;
			static const T ln10 = log(T(10));
			return Dual<T>(log10(x.value), x.deriv / (x.value * ln10));
		}

    /// Sine: d/dx(sin(x)) = cos(x)
    template<typename T>
    Dual<T> sin(const Dual<T>& x) {
			using std::cos;
			using std::sin;
			return Dual<T>(sin(x.value), cos(x.value) * x.deriv);
		}

    /// Cosine: d/dx(cos(x)) = -sin(x)
    template<typename T>
    Dual<T> cos(const Dual<T>& x) {
			using std::cos;
			using std::sin;
			return Dual<T>(cos(x.value), -sin(x.value) * x.deriv);
		}

    /// Tangent: d/dx(tan(x)) = sec²(x) = 1/cos²(x)
    template<typename T>
    Dual<T> tan(const Dual<T>& x) {
			using std::cos;
			using std::tan;
			T c = cos(x.value);
			return Dual<T>(tan(x.value), x.deriv / (c * c));
		}

    /// Arc sine: d/dx(asin(x)) = 1/sqrt(1-x²)
    template<typename T>
    Dual<T> asin(const Dual<T>& x) {
			using std::asin;
			using std::sqrt;
			T denom = sqrt(T(1) - x.value * x.value);
			return Dual<T>(asin(x.value), x.deriv / denom);
		}

    /// Arc cosine: d/dx(acos(x)) = -1/sqrt(1-x²)
    template<typename T>
    Dual<T> acos(const Dual<T>& x) {
			using std::acos;
			using std::sqrt;
			T denom = sqrt(T(1) - x.value * x.value);
			return Dual<T>(acos(x.value), -x.deriv / denom);
		}

    /// Arc tangent: d/dx(atan(x)) = 1/(1+x²)
    template<typename T>
    Dual<T> atan(const Dual<T>& x) {
			using std::atan;
			return Dual<T>(atan(x.value), x.deriv / (T(1) + x.value * x.value));
		}

    /// Hyperbolic sine: d/dx(sinh(x)) = cosh(x)
    template<typename T>
    Dual<T> sinh(const Dual<T>& x) {
			using std::cosh;
			using std::sinh;
			return Dual<T>(sinh(x.value), cosh(x.value) * x.deriv);
		}

    /// Hyperbolic cosine: d/dx(cosh(x)) = sinh(x)
    template<typename T>
    Dual<T> cosh(const Dual<T>& x) {
			using std::cosh;
			using std::sinh;
			return Dual<T>(cosh(x.value), sinh(x.value) * x.deriv);
		}

    /// Hyperbolic tangent: d/dx(tanh(x)) = sech²(x) = 1 - tanh²(x)
    template<typename T>
    Dual<T> tanh(const Dual<T>& x) {
			using std::tanh;
			T th = tanh(x.value);
			return Dual<T>(th, (T(1) - th * th) * x.deriv);
    }

    /// Absolute value: d/dx(|x|) = sign(x) = x/|x|
    template<typename T>
    Dual<T> abs(const Dual<T>& x) {
        if (x.value >= T(0)) {
            return x;
        } else {
            return Dual<T>(-x.value, -x.deriv);
        }
    }

    /// Floor function (derivative is 0 almost everywhere)
    template<typename T>
    Dual<T> floor(const Dual<T>& x) {
			using std::floor;
			return Dual<T>(floor(x.value), T(0));
		}

    /// Ceiling function (derivative is 0 almost everywhere)
    template<typename T>
    Dual<T> ceil(const Dual<T>& x) {
			using std::ceil;
			return Dual<T>(ceil(x.value), T(0));
		}

    /// Two-argument arc tangent: atan2(y, x)
    template<typename T>
    Dual<T> atan2(const Dual<T>& y, const Dual<T>& x) {
        // d/dx atan2(y,x) = -y/(x²+y²)
        // d/dy atan2(y,x) = x/(x²+y²)
        T denom = x.value * x.value + y.value * y.value;
        T d = (x.value * y.deriv - y.value * x.deriv) / denom;
				using std::atan2;
				return Dual<T>(atan2(y.value, x.value), d);
		}

    /// Fused multiply-add: fma(a, b, c) = a*b + c
    template<typename T>
    Dual<T> fma(const Dual<T>& a, const Dual<T>& b, const Dual<T>& c) {
        return a * b + c;
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    // Stream output
    //////////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    std::ostream& operator<<(std::ostream& os, const Dual<T>& d) {
        os << "Dual(" << d.value << ", " << d.deriv << ")";
        return os;
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    // Utility functions for derivative computation
    //////////////////////////////////////////////////////////////////////////////////////////

    /// Compute derivative of f at x using forward-mode AD
    /// @tparam F Callable type (function, lambda, functor)
    /// @param f Function to differentiate
    /// @param x Point at which to compute derivative
    /// @return Pair of (f(x), f'(x))
    template<typename F, typename T = Real, typename = std::enable_if_t<Detail::HasFloatingLeaf<T>::value>>
    std::pair<T, T> CalcDerivativeForwardAD(F&& f, T x) {
        Dual<T> xd(x, T(1));  // dx/dx = 1
        Dual<T> result = f(xd);
        return {result.value, result.deriv};
    }

    /// Compute gradient of f: R^n -> R at point x using forward-mode AD
    /// Requires n evaluations of f
    template<typename F, typename T = Real, typename = std::enable_if_t<Detail::HasFloatingLeaf<T>::value>>
    std::vector<T> CalcGradientForwardAD(F&& f, const std::vector<T>& x) {
        std::vector<T> grad(x.size());
        std::vector<Dual<T>> xd(x.size());
        
        for (size_t i = 0; i < x.size(); ++i) {
            // Initialize with zero derivative
            for (size_t j = 0; j < x.size(); ++j) {
                xd[j] = Dual<T>(x[j], T(0));
            }
            // Set derivative seed for variable i
            xd[i].deriv = T(1);
            
            // Evaluate
            Dual<T> result = f(xd);
            grad[i] = result.deriv;
        }
        
        return grad;
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    // Type alias for convenience
    //////////////////////////////////////////////////////////////////////////////////////////

    using DualR = Dual<Real>;
    using DualD = Dual<double>;
    using DualF = Dual<float>;


    // ---- Deprecated 2.0 spellings (CONTRACT_POLICY section 11) ----
    template<typename F, typename T = Real, typename = std::enable_if_t<Detail::HasFloatingLeaf<T>::value>>
    [[deprecated("use CalcDerivativeForwardAD()")]]
    std::pair<T, T> derivative(F&& f, T x) { return CalcDerivativeForwardAD(std::forward<F>(f), x); }
    template<typename F, typename T = Real, typename = std::enable_if_t<Detail::HasFloatingLeaf<T>::value>>
    [[deprecated("use CalcGradientForwardAD()")]]
    std::vector<T> gradient(F&& f, const std::vector<T>& x) { return CalcGradientForwardAD(std::forward<F>(f), x); }

} // namespace MML::AD

#endif // MML_SYMBOLIC_FORWARD_AD_H
