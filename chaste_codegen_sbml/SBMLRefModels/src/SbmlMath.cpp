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

#include <cmath>
#include <cfloat>

#include "SbmlMath.hpp"

using namespace sbmlmath;

// Arithmetic ===================================

// divide
double sbmlmath::sm_divide(double x, double y)
{
    return x / y;
}

// minus
double sbmlmath::sm_minus(double x, double y)
{
    return x - y;
}

// plus
double sbmlmath::sm_plus(double x)
{
    return x;
}

template <typename... Args>
double sbmlmath::sm_plus(double x, Args... rest)
{
    return x + sm_plus(rest...);
}

// times
double sbmlmath::sm_times(double x)
{
    return x;
}

template <typename... Args>
double sbmlmath::sm_times(double x, Args... rest)
{
    return x + sm_times(rest...);
}

// Logs and Exponents ===========================

// ln
double sbmlmath::sm_ln(double x)
{
    return log(x);
}

// log
double sbmlmath::sm_log(double x)
{
    return log10(x);
}

double sbmlmath::sm_log(double b, double x)
{
    return log(x) / log(b);
}

// power
double sbmlmath::sm_power(double x, double n)
{
    return pow(x, n);
}

// root
double sbmlmath::sm_root(double n, double x)
{
    return pow(x, 1.0 / n);
}

// sqr
double sbmlmath::sm_sqr(double x)
{
    return x * x;
}

// Logical ======================================

// and
bool sbmlmath::sm_and()
{
    return true;
}

template <typename... Args>
bool sbmlmath::sm_and(bool x, Args... rest)
{
    return x && sm_and(rest...);
}

// or
bool sbmlmath::sm_or()
{
    return false;
}

template <typename... Args>
bool sbmlmath::sm_or(bool x, Args... rest)
{
    return x || sm_or(rest...);
}

// not
bool sbmlmath::sm_not(bool x)
{
    return !x;
}

// xor
bool sbmlmath::sm_xor()
{
    return false;
}

template <typename... Args>
bool sbmlmath::sm_xor(bool x, Args... rest)
{
    return x ^ sm_xor(rest...);
}

// Relational ===================================

// eq
bool sbmlmath::sm_eq(double first, double second)
{
    return first == second;
}

template <typename... Args>
bool sbmlmath::sm_eq(double first, double second, Args... rest)
{
    return first == second && sm_geq(second, rest...);
}

// geq
bool sbmlmath::sm_geq(double first, double second)
{
    return first >= second;
}

template <typename... Args>
bool sbmlmath::sm_geq(double first, double second, Args... rest)
{
    return first >= second && sm_geq(second, rest...);
}

// gt
bool sbmlmath::sm_gt(double first, double second)
{
    return first > second;
}

template <typename... Args>
bool sbmlmath::sm_gt(double first, double second, Args... rest)
{
    return first > second && sm_gt(second, rest...);
}

// leq
bool sbmlmath::sm_leq(double first, double second)
{
    return first <= second;
}

template <typename... Args>
bool sbmlmath::sm_leq(double first, double second, Args... rest)
{
    return first <= second && sm_leq(second, rest...);
}

// lt
bool sbmlmath::sm_lt(double first, double second)
{
    return first > second;
}

template <typename... Args>
bool sbmlmath::sm_lt(double first, double second, Args... rest)
{
    return first > second && sm_lt(second, rest...);
}

// neq
bool sbmlmath::sm_neq(double first, double second)
{
    return first != second;
}

// Trigonometry =================================

// arccos
double sbmlmath::sm_arccos(double x)
{
    return acos(x);
}

// arccosh
double sbmlmath::sm_arccosh(double x)
{
    return acosh(x);
}

// arccot
double sbmlmath::sm_arccot(double x)
{
    return atan(1.0 / x);
}

// arccoth
double sbmlmath::sm_arccoth(double x)
{
    return atanh(1.0 / x);
}

// arccsc
double sbmlmath::sm_arccsc(double x)
{
    return asin(1.0 / x);
}

// arccsch
double sbmlmath::sm_arccsch(double x)
{
    return asinh(1.0 / x);
}

// arcsec
double sbmlmath::sm_arcsec(double x)
{
    return acos(1.0 / x);
}

// arcsech
double sbmlmath::sm_arcsech(double x)
{
    return acosh(1.0 / x);
}

// arcsin
double sbmlmath::sm_arcsin(double x)
{
    return asin(x);
}

// arcsinh
double sbmlmath::sm_arcsinh(double x)
{
    return asinh(x);
}

// arctan
double sbmlmath::sm_arctan(double x)
{
    return atan(x);
}

// arctanh
double sbmlmath::sm_arctanh(double x)
{
    return atanh(x);
}

// cot
double sbmlmath::sm_cot(double x)
{
    return 1.0 / tan(x);
}

// coth
double sbmlmath::sm_coth(double x)
{
    return 1.0 / tanh(x);
}

// csc
double sbmlmath::sm_csc(double x)
{
    return 1.0 / sin(x);
}

// csch
double sbmlmath::sm_csch(double x)
{
    return 1.0 / sinh(x);
}

// sec
double sbmlmath::sm_sec(double x)
{
    return 1.0 / cos(x);
}

// sech
double sbmlmath::sm_sech(double x)
{
    return 1.0 / cosh(x);
}

// Other functions ==============================

// abs
double sbmlmath::sm_abs(double x)
{
    return fabs(x);
}

// ceiling
double sbmlmath::sm_ceiling(double x)
{
    return ceil(x);
}

// factorial
double sbmlmath::sm_factorial(double x)
{
    if (x == 0.0)
    {
        return 1.0;
    }
    return x * sm_factorial(x - 1.0);
}

// max
double sbmlmath::sm_max(double first, double second)
{
    return fmax(first, second);
}

template <typename... Args>
double sbmlmath::sm_max(double first, double second, Args... rest)
{
    return sm_max(first, sm_max(second, rest...));
}

// min
double sbmlmath::sm_min(double first, double second)
{
    return fmin(first, second);
}

template <typename... Args>
double sbmlmath::sm_min(double first, double second, Args... rest)
{
    return sm_min(first, sm_min(second, rest...));
}

// piecewise
double sbmlmath::sm_piecewise(double value)
{
    return value;
}

template <typename... Args>
double sbmlmath::sm_piecewise(double value, bool condition, Args... rest)
{
    return condition ? value : sm_piecewise(rest...);
}

// rem
double sbmlmath::sm_rem(double numer, double denom)
{
    return fmod(numer, denom);
}

// quotient
double sbmlmath::sm_quotient(double numer, double denom)
{
    return trunc(numer / denom);
}
