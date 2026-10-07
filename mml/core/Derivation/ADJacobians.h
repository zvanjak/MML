///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ADJacobians.h                                                       ///
///  Description: Jacobian and Hessian helpers using automatic differentiation        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_DERIVATION_AD_JACOBIANS_H
#define MML_DERIVATION_AD_JACOBIANS_H

#include <mml/core/Derivation/ForwardAD.h>
#include <mml/core/Derivation/ReverseAD.h>

#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>

#include <vector>

namespace MML::AD
{
    namespace Detail
    {
        template<typename T>
        std::vector<T> ToStdVector(const Vector<T>& point)
        {
            std::vector<T> values(static_cast<size_t>(point.size()));
            for (int i = 0; i < point.size(); ++i)
                values[static_cast<size_t>(i)] = point[i];
            return values;
        }
    }

    template<typename F, typename T = Real>
    Matrix<T> jacobianForwardAD(F&& func, const std::vector<T>& point)
    {
        const size_t n = point.size();

        std::vector<Dual<T>> testInput(n);
        for (size_t i = 0; i < n; ++i)
            testInput[i] = Dual<T>(point[i], T(0));

        auto testOutput = func(testInput);
        const size_t m = testOutput.size();
        Matrix<T> jac(static_cast<int>(m), static_cast<int>(n));

        for (size_t j = 0; j < n; ++j)
        {
            std::vector<Dual<T>> input(n);
            for (size_t k = 0; k < n; ++k)
                input[k] = Dual<T>(point[k], k == j ? T(1) : T(0));

            auto output = func(input);
            for (size_t i = 0; i < m; ++i)
                jac(static_cast<int>(i), static_cast<int>(j)) = output[i].deriv;
        }

        return jac;
    }

    template<typename F, typename T = Real>
    Matrix<T> jacobianForwardAD(F&& func, const Vector<T>& point)
    {
        return jacobianForwardAD(std::forward<F>(func), Detail::ToStdVector(point));
    }

    template<typename F>
    Matrix<Real> jacobianReverseAD(F&& func, const std::vector<Real>& point)
    {
        const size_t n = point.size();

        Tape testTape;
        std::vector<ADVar> testInput;
        testInput.reserve(n);
        for (Real xi : point)
            testInput.emplace_back(&testTape, xi);

        auto testOutput = func(testInput);
        const size_t m = testOutput.size();
        Matrix<Real> jac(static_cast<int>(m), static_cast<int>(n));

        for (size_t i = 0; i < m; ++i)
        {
            Tape tape;
            std::vector<ADVar> input;
            input.reserve(n);
            for (Real xi : point)
                input.emplace_back(&tape, xi);

            auto output = func(input);
            output[i].backward();

            for (size_t j = 0; j < n; ++j)
                jac(static_cast<int>(i), static_cast<int>(j)) = input[j].adjoint();
        }

        return jac;
    }

    template<typename F>
    Matrix<Real> jacobianReverseAD(F&& func, const Vector<Real>& point)
    {
        return jacobianReverseAD(std::forward<F>(func), Detail::ToStdVector(point));
    }

    template<typename F, typename T = Real>
    Matrix<T> hessianForwardAD(F&& func, const std::vector<T>& point)
    {
        using Dual2 = Dual<Dual<T>>;
        const size_t n = point.size();
        Matrix<T> hessian(static_cast<int>(n), static_cast<int>(n));

        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = i; j < n; ++j)
            {
                std::vector<Dual2> input(n);
                for (size_t k = 0; k < n; ++k)
                {
                    input[k] = Dual2(
                        Dual<T>(point[k], k == j ? T(1) : T(0)),
                        Dual<T>(k == i ? T(1) : T(0), T(0))
                    );
                }

                Dual2 output = func(input);
                const T value = output.deriv.deriv;
                hessian(static_cast<int>(i), static_cast<int>(j)) = value;
                if (i != j)
                    hessian(static_cast<int>(j), static_cast<int>(i)) = value;
            }
        }

        return hessian;
    }

    template<typename F, typename T = Real>
    Matrix<T> hessianForwardAD(F&& func, const Vector<T>& point)
    {
        return hessianForwardAD(std::forward<F>(func), Detail::ToStdVector(point));
    }
}

#endif // MML_DERIVATION_AD_JACOBIANS_H