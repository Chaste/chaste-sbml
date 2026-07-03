/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SBML_MATH_HPP_
#define SBML_MATH_HPP_

/**
 * SBML math functions.
 */

#include <cmath>

namespace sbmlmath
{
// Constants ==================================

// SBML Level 3 recommended avogadro value:
// https://sbml.org/documents/specifications/level-3/
inline constexpr double AVOGADRO = 6.02214179E23;

// Note: Avogadro value has been updated in the most recent SI Brochure.
// Bureau International des Poids et Mesures (2019):
// The International System of Units (SI), 9th edition

// Arithmetic =================================

// divide
inline double divide(double x, double y) { return x / y; }

// minus
inline double minus(double x, double y) { return x - y; }

// plus
template <typename... Args>
constexpr double plus(Args... args);

// times
template <typename... Args>
constexpr double times(Args... args);

// Logs and exponents =========================

// log
inline double log(double x) { return std::log(x); }
inline double log(double b, double x) { return std::log(x) / std::log(b); }

// root
inline double root(double n, double x) { return std::pow(x, 1.0 / n); }

// sqr
inline double sqr(double x) { return x * x; }

// Logical ====================================

// and_
template <typename... Args>
constexpr bool and_(Args... args);

// or_
template <typename... Args>
constexpr bool or_(Args... args);

// not_
inline bool not_(bool x) { return !x; }

// xor_
template <typename... Args>
constexpr bool xor_(Args... args);

// Relational =================================

// eq
template <typename... Args>
constexpr bool eq(double first, double second, Args... rest);

// geq
template <typename... Args>
constexpr bool geq(double first, double second, Args... rest);

// gt
template <typename... Args>
constexpr bool gt(double first, double second, Args... rest);

// leq
template <typename... Args>
constexpr bool leq(double first, double second, Args... rest);

// lt
template <typename... Args>
constexpr bool lt(double first, double second, Args... rest);

// neq
inline bool neq(double x, double y) { return x != y; }

// Trigonometry ===============================

// cot, coth, acot, acoth
inline double cot(double x) { return 1.0 / std::tan(x); }
inline double coth(double x) { return 1.0 / std::tanh(x); }
inline double acot(double x) { return std::atan(1.0 / x); }
inline double acoth(double x) { return std::atanh(1.0 / x); }

// csc, csch, acsc, acsch
inline double csc(double x) { return 1.0 / std::sin(x); }
inline double csch(double x) { return 1.0 / std::sinh(x); }
inline double acsc(double x) { return std::asin(1.0 / x); }
inline double acsch(double x) { return std::asinh(1.0 / x); }

// sec, sech, asec, asech
inline double sec(double x) { return 1.0 / std::cos(x); }
inline double sech(double x) { return 1.0 / std::cosh(x); }
inline double asec(double x) { return std::acos(1.0 / x); }
inline double asech(double x) { return std::acosh(1.0 / x); }

// Other functions ============================

// factorial
inline double factorial(double x) { return std::tgamma(x + 1.0); }

// max
template <typename... Args>
constexpr double max(double first, double second, Args... rest);

// min
template <typename... Args>
constexpr double min(double first, double second, Args... rest);

// piecewise
constexpr double piecewise(double otherwise);

constexpr double piecewise(double value, bool condition, double otherwise);

template <typename... Args>
constexpr double piecewise(double value, bool condition, Args... rest);

// quotient
inline double quotient(double numer, double denom) { return std::trunc(numer / denom); }

} // namespace sbmlmath

// Arithmetic (variadic templates) ============

// plus
template <typename... Args>
constexpr double sbmlmath::plus(Args... args)
{
    return (0.0 + ... + args);
}

// times
template <typename... Args>
constexpr double sbmlmath::times(Args... args)
{
    return (1.0 * ... * args);
}

// Logical (variadic templates) =================

// and_
template <typename... Args>
constexpr bool sbmlmath::and_(Args... args)
{
    return (true && ... && args);
}

// or_
template <typename... Args>
constexpr bool sbmlmath::or_(Args... args)
{
    return (false || ... || args);
}

// xor_
template <typename... Args>
constexpr bool sbmlmath::xor_(Args... args)
{
    return (false ^ ... ^ args);
}

// Relational (variadic templates) ==============

// eq
template <typename... Args>
constexpr bool sbmlmath::eq(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return (first == second) && sbmlmath::eq(second, rest...);
    }
    return first == second;
}

// geq
template <typename... Args>
constexpr bool sbmlmath::geq(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return (first >= second) && sbmlmath::geq(second, rest...);
    }
    return first >= second;
}

// gt
template <typename... Args>
constexpr bool sbmlmath::gt(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return (first > second) && sbmlmath::gt(second, rest...);
    }
    return first > second;
}

// leq
template <typename... Args>
constexpr bool sbmlmath::leq(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return (first <= second) && sbmlmath::leq(second, rest...);
    }
    return first <= second;
}

// lt
template <typename... Args>
constexpr bool sbmlmath::lt(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return (first < second) && sbmlmath::lt(second, rest...);
    }
    return first < second;
}

// Other (variadic templates) ===================

// max
template <typename... Args>
constexpr double sbmlmath::max(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return std::fmax(first, sbmlmath::max(second, rest...));
    }
    return std::fmax(first, second);
}

// min
template <typename... Args>
constexpr double sbmlmath::min(double first, double second, Args... rest)
{
    if constexpr (sizeof...(rest) > 0)
    {
        return std::fmin(first, sbmlmath::min(second, rest...));
    }
    return std::fmin(first, second);
}

// piecewise
constexpr double sbmlmath::piecewise(double otherwise)
{
    return otherwise;
}

constexpr double sbmlmath::piecewise(double value, bool condition, double otherwise)
{
    return condition ? value : otherwise;
}

template <typename... Args>
constexpr double sbmlmath::piecewise(double value, bool condition, Args... rest)
{
    return condition ? value : sbmlmath::piecewise(rest...);
}

#endif // SBML_MATH_HPP_
