///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        ADArithmetic.h                                                       ///
///  Description: Shared quotient-tangent arithmetic for automatic differentiation     ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                          ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_AD_ARITHMETIC_H
#define MML_AD_ARITHMETIC_H

#include <mml/MMLExceptions.h>

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace MML::AD
{
    template<typename T> struct Dual;
}

namespace MML::AD::Detail
{
    template<typename T> struct HasFloatingLeaf : std::is_floating_point<T> {};
    template<typename T> struct HasFloatingLeaf<Dual<T>> : HasFloatingLeaf<T> {};

    template<typename T>
    T Primal(const T& x)
    {
        return x;
    }

    template<typename T>
    auto Primal(const Dual<T>& x)
    {
        return Primal(x.value);
    }

    template<typename T>
    bool IsZero(const T& x)
    {
        return x == T(0);
    }

    template<typename T>
    bool IsZero(const Dual<T>& x)
    {
        return IsZero(x.value) && IsZero(x.deriv);
    }

    template<typename T>
    bool IsConstantJet(const T& x)
    {
        if constexpr (std::is_floating_point_v<T>)
            return true;
        else
            return IsZero(x.deriv) && IsConstantJet(x.value);
    }

    template<typename T>
    void RequirePowerDomain(T base, T exponent, bool fixedExponent)
    {
        if (!std::isfinite(base) || !std::isfinite(exponent))
            throw DomainError("AD pow: base and exponent must be finite");
        if (!fixedExponent && base <= T(0))
            throw DomainError("AD pow: variable exponent requires a positive base");
        if (base == T(0) && exponent < T(0))
            throw DomainError("AD pow: zero base requires a nonnegative exponent");
        if (base < T(0) && std::trunc(exponent) != exponent)
            throw DomainError("AD pow: negative base requires a fixed integer exponent");
    }

    // Validate the expression once; singular derivative intermediates retain IEEE results.
    template<typename T>
    T PowJetUnchecked(const T& base, const T& exponent)
    {
        if constexpr (std::is_floating_point_v<T>)
            return std::pow(base, exponent);
        else
        {
            using Value = std::decay_t<decltype(base.value)>;
            using std::log;
            if (IsZero(exponent.deriv))
            {
                if (IsZero(exponent.value))
                    return T(Value(1), Value(0));
                const Value p = PowJetUnchecked(base.value, exponent.value);
                if (IsZero(base.deriv))
                    return T(p, Value(0));
                return T(p, exponent.value * PowJetUnchecked(base.value, exponent.value - Value(1)) * base.deriv);
            }
            const Value p = PowJetUnchecked(base.value, exponent.value);
            if (IsZero(base.deriv))
                return T(p, p * log(base.value) * exponent.deriv);
            return T(p, p * (exponent.deriv * log(base.value) + exponent.value * base.deriv / base.value));
        }
    }

    template<typename T>
    bool IsFiniteJet(const T& x)
    {
        if constexpr (std::is_floating_point_v<T>)
            return std::isfinite(x);
        else
            return IsFiniteJet(x.value) && IsFiniteJet(x.deriv);
    }

    template<typename T>
    int JetExponent(const T& x)
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            int exponent;
            std::frexp(x, &exponent);
            return exponent;
        }
        else
        {
            if (IsZero(x.value)) return JetExponent(x.deriv);
            if (IsZero(x.deriv)) return JetExponent(x.value);
            return std::max(JetExponent(x.value), JetExponent(x.deriv));
        }
    }

    template<typename T>
    T ScaleJet(const T& x, int exponent)
    {
        if constexpr (std::is_floating_point_v<T>)
            return std::scalbn(x, exponent);
        else
            return T(ScaleJet(x.value, exponent), ScaleJet(x.deriv, exponent));
    }

    template<typename T>
    T FmaJet(const T& a, const T& b, const T& c)
    {
        if constexpr (std::is_floating_point_v<T>)
            return std::fma(a, b, c);
        else
            return T(FmaJet(a.value, b.value, c.value),
                FmaJet(a.deriv, b.value, FmaJet(a.value, b.deriv, c.deriv)));
    }

    template<typename T>
    T QuotientTangent(const T& a, const T& da, const T& b, const T& db)
    {
        if constexpr (HasFloatingLeaf<T>::value)
        {
            if (!IsFiniteJet(a) || !IsFiniteJet(da) || !IsFiniteJet(b)
                || !IsFiniteJet(db) || Primal(b) == 0)
                return (da * b - a * db) / (b * b);

            const bool firstZero = IsZero(da), secondZero = IsZero(a) || IsZero(db);
            if (firstZero && secondZero)
                return T(0);

            const int aExponent = JetExponent(a), daExponent = JetExponent(da);
            const int bExponent = JetExponent(b), dbExponent = JetExponent(db);
            const T aMantissa = ScaleJet(a, -aExponent);
            const T daMantissa = ScaleJet(da, -daExponent);
            const T bMantissa = ScaleJet(b, -bExponent);
            const T dbMantissa = ScaleJet(db, -dbExponent);
            const int firstExponent = daExponent + bExponent;
            const int secondExponent = aExponent + dbExponent;
            const int exponent = firstZero ? secondExponent
                : secondZero ? firstExponent : std::max(firstExponent, secondExponent);
            const T first = firstZero ? T(0) : ScaleJet(daMantissa, firstExponent - exponent);
            const T second = secondZero ? T(0) : ScaleJet(aMantissa, secondExponent - exponent);

            // Align products before cancellation; FMA retains the rounded product's residual.
            const T product = second * dbMantissa;
            const T residual = FmaJet(-second, dbMantissa, product);
            const T numerator = FmaJet(first, bMantissa, -product) + residual;
            if (IsZero(numerator))
                return T(0);
            int denominatorExponent;
            std::frexp(Primal(b), &denominatorExponent);
            const T denominator = ScaleJet(b, -denominatorExponent);
            return ScaleJet((numerator / denominator) / denominator, exponent - 2 * denominatorExponent);
        }
        else
        {
            return (da * b - a * db) / (b * b);
        }
    }
}

#endif // MML_AD_ARITHMETIC_H
