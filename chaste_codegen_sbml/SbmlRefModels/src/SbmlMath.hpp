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

#ifndef SBMLMATH_HPP_
#define SBMLMATH_HPP_

/**
 * SBML math functions.
 */

#include <cmath>

namespace sbmlmath
{
  // Constants ==================================

  // SBML Level 3 recommended avogadro value:
  // https://sbml.org/documents/specifications/level-3/
  constexpr double SM_AVOGADRO = 6.02214179E23;

  // Note: Avogadro value has been updated in the most recent SI Brochure.
  // Bureau International des Poids et Mesures (2019):
  // The International System of Units (SI), 9th edition

  // Arithmetic =================================

  // divide
  double sm_divide(double x, double y);

  // minus
  double sm_minus(double x, double y);

  // plus
  template <typename... Args>
  constexpr double sm_plus(Args... args)
  {
    return (0.0 + ... + args);
  }

  // times
  template <typename... Args>
  constexpr double sm_times(Args... args)
  {
    return (1.0 * ... * args);
  }

  // Logs and exponents =========================

  // log
  double sm_log(double x);
  double sm_log(double b, double x);

  // root
  double sm_root(double n, double x);

  // sqr
  double sm_sqr(double x);

  // Logical ====================================

  // and
  template <typename... Args>
  constexpr bool sm_and(Args... args)
  {
    return (true && ... && args);
  }

  // or
  template <typename... Args>
  constexpr bool sm_or(Args... args)
  {
    return (false || ... || args);
  }

  // not
  bool sm_not(bool x);

  // xor
  template <typename... Args>
  constexpr bool sm_xor(Args... args)
  {
    return (false ^ ... ^ args);
  }

  // Relational =================================

  // eq
  constexpr bool sm_eq(double first, double second)
  {
    return first == second;
  }

  template <typename... Args>
  constexpr bool sm_eq(double first, double second, Args... rest)
  {
    return sm_eq(first, second) && sm_eq(second, rest...);
  }

  // geq
  constexpr bool sm_geq(double first, double second)
  {
    return first >= second;
  }

  template <typename... Args>
  constexpr bool sm_geq(double first, double second, Args... rest)
  {
    return sm_geq(first, second) && sm_geq(second, rest...);
  }

  // gt
  constexpr bool sm_gt(double first, double second)
  {
    return first > second;
  }

  template <typename... Args>
  constexpr bool sm_gt(double first, double second, Args... rest)
  {
    return sm_gt(first, second) && sm_gt(second, rest...);
  }

  // leq
  constexpr bool sm_leq(double first, double second)
  {
    return first <= second;
  }

  template <typename... Args>
  constexpr bool sm_leq(double first, double second, Args... rest)
  {
    return sm_leq(first, second) && sm_leq(second, rest...);
  }

  // lt
  constexpr bool sm_lt(double first, double second)
  {
    return first < second;
  }

  template <typename... Args>
  constexpr bool sm_lt(double first, double second, Args... rest)
  {
    return sm_lt(first, second) && sm_lt(second, rest...);
  }

  // neq
  bool sm_neq(double x, double y);

  // Trigonometry ===============================

  // cot, coth, acot, acoth
  double sm_cot(double x);
  double sm_coth(double x);
  double sm_acot(double x);
  double sm_acoth(double x);

  // csc, csch, acsc, acsch
  double sm_csc(double x);
  double sm_csch(double x);
  double sm_acsc(double x);
  double sm_acsch(double x);

  // sec, sech, asec, asech
  double sm_sec(double x);
  double sm_sech(double x);
  double sm_asec(double x);
  double sm_asech(double x);

  // Other functions ============================

  // factorial
  double sm_factorial(double x);

  // max
  template <typename... Args>
  constexpr double sm_max(double first, double second)
  {
    return std::fmax(first, second);
  }

  template <typename... Args>
  constexpr double sm_max(double first, double second, Args... rest)
  {
    return std::fmax(first, sm_max(second, rest...));
  }

  // min
  template <typename... Args>
  constexpr double sm_min(double first, double second)
  {
    return std::fmin(first, second);
  }

  template <typename... Args>
  constexpr double sm_min(double first, double second, Args... rest)
  {
    return std::fmin(first, sm_min(second, rest...));
  }

  // piecewise
  constexpr double sm_piecewise(double value0, bool condition, double value1)
  {
    return condition ? value0 : value1;
  }

  template <typename... Args>
  constexpr double sm_piecewise(double value, bool condition, Args... rest)
  {
    return condition ? value : sm_piecewise(rest...);
  }

  // quotient
  double sm_quotient(double numer, double denom);

} // namespace sbmlmath

#endif // SBMLMATH_HPP_
