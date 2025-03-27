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

namespace sbmlmath
{

  // Arithmetic ===============================

  // divide
  double sm_divide(double x, double y);

  // minus
  double sm_minus(double x, double y);

  // plus
  double sm_plus(double x);

  template <typename... Args>
  double sm_plus(double x, Args... rest);

  // times
  double sm_times(double x);

  template <typename... Args>
  double sm_times(double x, Args... rest);

  // Logs and exponents =======================

  // ln
  double sm_ln(double x);

  // log
  double sm_log(double x);
  double sm_log(double b, double x);

  // power
  double sm_power(double x, double n);

  // root
  double sm_root(double n, double x);

  // sqr
  double sm_sqr(double x);

  // Logical ==================================

  // and
  bool sm_and();

  template <typename... Args>
  bool sm_and(bool x, Args... rest);

  // or
  bool sm_or();

  template <typename... Args>
  bool sm_or(bool x, Args... rest);

  // not
  bool sm_not(bool x);

  // xor
  bool sm_xor();

  template <typename... Args>
  bool sm_xor(bool x, Args... rest);

  // Relational ================================

  // eq
  bool sm_eq(double first, double second);

  template <typename... Args>
  bool sm_eq(double first, double second, Args... rest);

  // geq
  bool sm_geq(double first, double second);

  template <typename... Args>
  bool sm_geq(double first, double second, Args... rest);

  // gt
  bool sm_gt(double first, double second);

  template <typename... Args>
  bool sm_gt(double first, double second, Args... rest);

  // leq
  bool sm_leq(double first, double second);

  template <typename... Args>
  bool sm_leq(double first, double second, Args... rest);

  // lt
  bool sm_lt(double first, double second);

  template <typename... Args>
  bool sm_lt(double first, double second, Args... rest);

  // neq
  bool sm_neq(double first, double second);

  // Trigonometry =============================

  double sm_arccos(double x);
  double sm_arccosh(double x);

  double sm_arccot(double x);
  double sm_arccoth(double x);

  double sm_arccsc(double x);
  double sm_arccsch(double x);

  double sm_arcsec(double x);
  double sm_arcsech(double x);

  double sm_arcsin(double x);
  double sm_arcsinh(double x);

  double sm_arctan(double x);
  double sm_arctanh(double x);

  double sm_cot(double x);
  double sm_coth(double x);

  double sm_csc(double x);
  double sm_csch(double x);

  double sm_sec(double x);
  double sm_sech(double x);

  // Other functions ==========================

  // abs
  double sm_abs(double x);

  // ceiling
  double sm_ceiling(double x);

  // factorial
  double sm_factorial(double x);

  // max
  double sm_max(double first, double second);

  template <typename... Args>
  double sm_max(double first, double second, Args... rest);

  // min
  double sm_min(double first, double second);

  template <typename... Args>
  double sm_min(double first, double second, Args... rest);

  // piecewise
  double sm_piecewise(double value);

  template <typename... Args>
  double sm_piecewise(double value, bool condition, Args... rest);

  // rem
  double sm_rem(double numer, double denom);

  // quotient
  double sm_quotient(double x, double y);

  // TODO: Other MathML elements permitted in SBML Level 3
  // annotation, annotation-xml, apply, bvar, ci, cn, csymbol, degree,
  // exponentiale, false, implies, infinity, lambda, logbase, notanumber,
  // otherwise, piece, pi, quotient, rem, root, semantics, sep, true

}

#endif // SBMLMATH_HPP_
