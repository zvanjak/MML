///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ReverseAD.h                                                         ///
///  Description: Reverse-mode automatic differentiation using a computation tape      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_AD_REVERSE_AD_H
#define MML_AD_REVERSE_AD_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <cmath>
#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <vector>

namespace MML::AD
{
    enum class ADOp : uint8_t
    {
        CONSTANT,
        VAR,
        ADD,
        SUB,
        MUL,
        DIV,
        NEG,
        SIN,
        COS,
        TAN,
        ASIN,
        ACOS,
        ATAN,
        SINH,
        COSH,
        TANH,
        EXP,
        LOG,
        SQRT,
        POW,
        ABS,
    };

    struct TapeEntry
    {
        ADOp op;
        int arg1;
        int arg2;
        Real value;
        Real adjoint;
    };

    class Tape
    {
    private:
        std::vector<TapeEntry> entries_;
        bool backward_complete_ = false;

    public:
        void clear()
        {
            entries_.clear();
            backward_complete_ = false;
        }

        size_t size() const { return entries_.size(); }
        bool empty() const { return entries_.empty(); }

        TapeEntry& operator[](size_t i) { return entries_[i]; }
        const TapeEntry& operator[](size_t i) const { return entries_[i]; }

        int push(ADOp op, int arg1, int arg2, Real value)
        {
            const int idx = static_cast<int>(entries_.size());
            entries_.push_back({ op, arg1, arg2, value, 0.0 });
            backward_complete_ = false;
            return idx;
        }

        void backward(int outputIdx)
        {
            if (outputIdx < 0 || outputIdx >= static_cast<int>(entries_.size()))
                throw ArgumentError("Tape::backward: invalid output index");

            for (auto& entry : entries_)
                entry.adjoint = 0.0;

            entries_[outputIdx].adjoint = 1.0;

            for (int i = outputIdx; i >= 0; --i)
            {
                const auto& entry = entries_[i];
                const Real adj = entry.adjoint;
                if (adj == 0.0)
                    continue;

                switch (entry.op)
                {
                case ADOp::CONSTANT:
                case ADOp::VAR:
                    break;
                case ADOp::ADD:
                    entries_[entry.arg1].adjoint += adj;
                    entries_[entry.arg2].adjoint += adj;
                    break;
                case ADOp::SUB:
                    entries_[entry.arg1].adjoint += adj;
                    entries_[entry.arg2].adjoint -= adj;
                    break;
                case ADOp::MUL:
                    entries_[entry.arg1].adjoint += adj * entries_[entry.arg2].value;
                    entries_[entry.arg2].adjoint += adj * entries_[entry.arg1].value;
                    break;
                case ADOp::DIV:
                {
                    const Real denom = entries_[entry.arg2].value;
                    entries_[entry.arg1].adjoint += adj / denom;
                    entries_[entry.arg2].adjoint -= adj * entries_[entry.arg1].value / (denom * denom);
                    break;
                }
                case ADOp::NEG:
                    entries_[entry.arg1].adjoint -= adj;
                    break;
                case ADOp::SIN:
                    entries_[entry.arg1].adjoint += adj * std::cos(entries_[entry.arg1].value);
                    break;
                case ADOp::COS:
                    entries_[entry.arg1].adjoint -= adj * std::sin(entries_[entry.arg1].value);
                    break;
                case ADOp::TAN:
                {
                    const Real c = std::cos(entries_[entry.arg1].value);
                    entries_[entry.arg1].adjoint += adj / (c * c);
                    break;
                }
                case ADOp::ASIN:
                {
                    const Real a = entries_[entry.arg1].value;
                    entries_[entry.arg1].adjoint += adj / std::sqrt(1.0 - a * a);
                    break;
                }
                case ADOp::ACOS:
                {
                    const Real a = entries_[entry.arg1].value;
                    entries_[entry.arg1].adjoint -= adj / std::sqrt(1.0 - a * a);
                    break;
                }
                case ADOp::ATAN:
                {
                    const Real a = entries_[entry.arg1].value;
                    entries_[entry.arg1].adjoint += adj / (1.0 + a * a);
                    break;
                }
                case ADOp::SINH:
                    entries_[entry.arg1].adjoint += adj * std::cosh(entries_[entry.arg1].value);
                    break;
                case ADOp::COSH:
                    entries_[entry.arg1].adjoint += adj * std::sinh(entries_[entry.arg1].value);
                    break;
                case ADOp::TANH:
                {
                    const Real th = std::tanh(entries_[entry.arg1].value);
                    entries_[entry.arg1].adjoint += adj * (1.0 - th * th);
                    break;
                }
                case ADOp::EXP:
                    entries_[entry.arg1].adjoint += adj * entry.value;
                    break;
                case ADOp::LOG:
                    entries_[entry.arg1].adjoint += adj / entries_[entry.arg1].value;
                    break;
                case ADOp::SQRT:
                    entries_[entry.arg1].adjoint += adj / (2.0 * entry.value);
                    break;
                case ADOp::POW:
                {
                    const Real a = entries_[entry.arg1].value;
                    const Real b = entries_[entry.arg2].value;
                    entries_[entry.arg1].adjoint += adj * b * std::pow(a, b - 1.0);
                    if (a > 0.0)
                        entries_[entry.arg2].adjoint += adj * entry.value * std::log(a);
                    break;
                }
                case ADOp::ABS:
                {
                    const Real a = entries_[entry.arg1].value;
                    entries_[entry.arg1].adjoint += adj * (a >= 0.0 ? 1.0 : -1.0);
                    break;
                }
                }
            }

            backward_complete_ = true;
        }

        bool isBackwardComplete() const { return backward_complete_; }
    };

    class ADVar
    {
    private:
        Tape* tape_ = nullptr;
        int idx_ = -1;
        Real constant_value_ = 0.0;

        static Tape* commonTape(const ADVar& a, const ADVar& b)
        {
            if (a.tape_ && b.tape_ && a.tape_ != b.tape_)
                throw ArgumentError("ADVar: operands belong to different tapes");
            return a.tape_ ? a.tape_ : b.tape_;
        }

        int indexOn(Tape* tape) const
        {
            if (!tape)
                return -1;
            if (tape_ == tape)
                return idx_;
            if (tape_ != nullptr)
                throw ArgumentError("ADVar: operand belongs to a different tape");
            return tape->push(ADOp::CONSTANT, -1, -1, constant_value_);
        }

        struct TapeIndexTag {};

        ADVar(Tape* tape, int idx, TapeIndexTag) : tape_(tape), idx_(idx) {}

    public:
        static ADVar makeUnary(const ADVar& a, ADOp op, Real value)
        {
            if (!a.tape_)
                return ADVar(value);
            return ADVar(a.tape_, a.tape_->push(op, a.idx_, -1, value), TapeIndexTag{});
        }

        static ADVar makeBinary(const ADVar& a, const ADVar& b, ADOp op, Real value)
        {
            Tape* tape = commonTape(a, b);
            if (!tape)
                return ADVar(value);
            return ADVar(tape, tape->push(op, a.indexOn(tape), b.indexOn(tape), value), TapeIndexTag{});
        }

        explicit ADVar(Real value = 0.0) : constant_value_(value) {}

        ADVar(Tape* tape, Real value) : tape_(tape)
        {
            if (!tape_)
                throw ArgumentError("ADVar: input variable requires a tape");
            idx_ = tape_->push(ADOp::VAR, -1, -1, value);
        }

        static ADVar constant(Tape* tape, Real value)
        {
            if (!tape)
                return ADVar(value);
            return ADVar(tape, tape->push(ADOp::CONSTANT, -1, -1, value), TapeIndexTag{});
        }

        Real value() const
        {
            if (!tape_)
                return constant_value_;
            return (*tape_)[idx_].value;
        }

        Real adjoint() const
        {
            if (!tape_)
                return 0.0;
            return (*tape_)[idx_].adjoint;
        }

        int index() const { return idx_; }
        Tape* tape() const { return tape_; }

        void backward() const
        {
            if (!tape_)
                throw ArgumentError("ADVar::backward: constants are not recorded on a tape");
            tape_->backward(idx_);
        }

        friend ADVar operator+(const ADVar& a, const ADVar& b) { return makeBinary(a, b, ADOp::ADD, a.value() + b.value()); }
        friend ADVar operator-(const ADVar& a, const ADVar& b) { return makeBinary(a, b, ADOp::SUB, a.value() - b.value()); }
        friend ADVar operator*(const ADVar& a, const ADVar& b) { return makeBinary(a, b, ADOp::MUL, a.value() * b.value()); }
        friend ADVar operator/(const ADVar& a, const ADVar& b)
        {
            const Real denom = b.value();
            if (denom == 0.0)
                throw DivisionByZeroError("ADVar: division by zero");
            return makeBinary(a, b, ADOp::DIV, a.value() / denom);
        }
        friend ADVar operator-(const ADVar& a) { return makeUnary(a, ADOp::NEG, -a.value()); }
        friend ADVar operator+(const ADVar& a) { return a; }

        friend ADVar operator+(const ADVar& a, Real b) { return a + ADVar(b); }
        friend ADVar operator+(Real a, const ADVar& b) { return ADVar(a) + b; }
        friend ADVar operator-(const ADVar& a, Real b) { return a - ADVar(b); }
        friend ADVar operator-(Real a, const ADVar& b) { return ADVar(a) - b; }
        friend ADVar operator*(const ADVar& a, Real b) { return a * ADVar(b); }
        friend ADVar operator*(Real a, const ADVar& b) { return ADVar(a) * b; }
        friend ADVar operator/(const ADVar& a, Real b) { return a / ADVar(b); }
        friend ADVar operator/(Real a, const ADVar& b) { return ADVar(a) / b; }

        ADVar& operator+=(const ADVar& other) { *this = *this + other; return *this; }
        ADVar& operator-=(const ADVar& other) { *this = *this - other; return *this; }
        ADVar& operator*=(const ADVar& other) { *this = *this * other; return *this; }
        ADVar& operator/=(const ADVar& other) { *this = *this / other; return *this; }
        ADVar& operator+=(Real other) { *this = *this + other; return *this; }
        ADVar& operator-=(Real other) { *this = *this - other; return *this; }
        ADVar& operator*=(Real other) { *this = *this * other; return *this; }
        ADVar& operator/=(Real other) { *this = *this / other; return *this; }
    };

    inline ADVar sin(const ADVar& x) { return ADVar::makeUnary(x, ADOp::SIN, std::sin(x.value())); }
    inline ADVar cos(const ADVar& x) { return ADVar::makeUnary(x, ADOp::COS, std::cos(x.value())); }
    inline ADVar tan(const ADVar& x) { return ADVar::makeUnary(x, ADOp::TAN, std::tan(x.value())); }
    inline ADVar asin(const ADVar& x) { return ADVar::makeUnary(x, ADOp::ASIN, std::asin(x.value())); }
    inline ADVar acos(const ADVar& x) { return ADVar::makeUnary(x, ADOp::ACOS, std::acos(x.value())); }
    inline ADVar atan(const ADVar& x) { return ADVar::makeUnary(x, ADOp::ATAN, std::atan(x.value())); }
    inline ADVar sinh(const ADVar& x) { return ADVar::makeUnary(x, ADOp::SINH, std::sinh(x.value())); }
    inline ADVar cosh(const ADVar& x) { return ADVar::makeUnary(x, ADOp::COSH, std::cosh(x.value())); }
    inline ADVar tanh(const ADVar& x) { return ADVar::makeUnary(x, ADOp::TANH, std::tanh(x.value())); }
    inline ADVar exp(const ADVar& x) { return ADVar::makeUnary(x, ADOp::EXP, std::exp(x.value())); }
    inline ADVar log(const ADVar& x)
    {
        if (x.value() <= 0.0)
            throw DomainError("ADVar log: argument must be positive");
        return ADVar::makeUnary(x, ADOp::LOG, std::log(x.value()));
    }
    inline ADVar sqrt(const ADVar& x)
    {
        if (x.value() < 0.0)
            throw DomainError("ADVar sqrt: argument must be non-negative");
        return ADVar::makeUnary(x, ADOp::SQRT, std::sqrt(x.value()));
    }
    inline ADVar abs(const ADVar& x) { return ADVar::makeUnary(x, ADOp::ABS, std::abs(x.value())); }
    inline ADVar pow(const ADVar& base, const ADVar& exp) { return ADVar::makeBinary(base, exp, ADOp::POW, std::pow(base.value(), exp.value())); }
    inline ADVar pow(const ADVar& base, Real exp) { return pow(base, ADVar(exp)); }
    inline ADVar pow(Real base, const ADVar& exp) { return pow(ADVar(base), exp); }

    inline std::ostream& operator<<(std::ostream& os, const ADVar& v)
    {
        os << "ADVar(value=" << v.value() << ", adjoint=" << v.adjoint() << ")";
        return os;
    }

    template<typename F>
    std::vector<Real> gradientReverse(F&& f, const std::vector<Real>& x)
    {
        Tape tape;
        std::vector<ADVar> variables;
        variables.reserve(x.size());

        for (Real xi : x)
            variables.emplace_back(&tape, xi);

        ADVar result = f(variables);
        result.backward();

        std::vector<Real> grad(x.size());
        for (size_t i = 0; i < x.size(); ++i)
            grad[i] = variables[i].adjoint();

        return grad;
    }
}

#endif // MML_AD_REVERSE_AD_H
