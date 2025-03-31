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

#ifndef TESTSBMLMATH_HPP_
#define TESTSBMLMATH_HPP_

#include <cmath>

#include <cxxtest/TestSuite.h>

#include "SbmlMath.hpp"

using namespace sbmlmath;

class TestSbmlMath : public CxxTest::TestSuite
{
public:
  // Arithmetic =================================

  void TestDivide()
  {
    TS_ASSERT_EQUALS(sm_divide(1.0, 1.0), 1.0);
    TS_ASSERT_EQUALS(sm_divide(1.0, 2.0), 0.5);
    TS_ASSERT_EQUALS(sm_divide(2.0, 1.0), 2.0);
  }

  void TestMinus()
  {
    TS_ASSERT_EQUALS(sm_minus(1.0, 1.0), 0.0);
    TS_ASSERT_EQUALS(sm_minus(1.0, 2.0), -1.0);
    TS_ASSERT_EQUALS(sm_minus(2.0, 1.0), 1.0);
  }

  void TestPlus()
  {
    TS_ASSERT_EQUALS(sm_plus(), 0.0);
    TS_ASSERT_EQUALS(sm_plus(1.0), 1.0);
    TS_ASSERT_EQUALS(sm_plus(1.0, 2.0), 3.0);
    TS_ASSERT_EQUALS(sm_plus(1.0, 2.0, 3.0), 6.0);
  }

  void TestTimes()
  {
    TS_ASSERT_EQUALS(sm_times(), 1.0);
    TS_ASSERT_EQUALS(sm_times(1.0), 1.0);
    TS_ASSERT_EQUALS(sm_times(2.0), 2.0);
    TS_ASSERT_EQUALS(sm_times(1.0, 2.0), 2.0);
    TS_ASSERT_EQUALS(sm_times(1.0, 2.0, 3.0), 6.0);
  }

  // Logs and exponents =========================

  void TestLog()
  {
    TS_ASSERT_DELTA(sm_log(1.0), 0.0, 1e-6);
    TS_ASSERT_DELTA(sm_log(10.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_log(M_E), 0.434294, 1e-6);

    TS_ASSERT_DELTA(sm_log(2.0, 2.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_log(10.0, 2.0), std::log10(2.0), 1e-6);
    TS_ASSERT_DELTA(sm_log(M_E, 2.0), M_LN2, 1e-6);
  }

  void TestRoot()
  {
    TS_ASSERT_DELTA(sm_root(1.0, 1.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_root(2.0, 4.0), 2.0, 1e-6);
    TS_ASSERT_DELTA(sm_root(3.0, 8.0), 2.0, 1e-6);
    TS_ASSERT_DELTA(sm_root(2.0, 9.0), 3.0, 1e-6);
    TS_ASSERT_DELTA(sm_root(3.0, 27.0), 3.0, 1e-6);
  }

  void TestSqr()
  {
    TS_ASSERT_DELTA(sm_sqr(1.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_sqr(2.0), 4.0, 1e-6);
    TS_ASSERT_DELTA(sm_sqr(3.0), 9.0, 1e-6);
  }

  // Logical ===================================

  void TestAnd()
  {
    TS_ASSERT_EQUALS(sm_and(), true);
    TS_ASSERT_EQUALS(sm_and(true), true);
    TS_ASSERT_EQUALS(sm_and(false), false);
    TS_ASSERT_EQUALS(sm_and(true, true), true);
    TS_ASSERT_EQUALS(sm_and(true, false), false);
    TS_ASSERT_EQUALS(sm_and(false, true), false);
    TS_ASSERT_EQUALS(sm_and(false, false), false);
    TS_ASSERT_EQUALS(sm_and(true, true, true), true);
    TS_ASSERT_EQUALS(sm_and(false, true, true), false);
  };

  void TestOr()
  {
    TS_ASSERT_EQUALS(sm_or(), false);
    TS_ASSERT_EQUALS(sm_or(true), true);
    TS_ASSERT_EQUALS(sm_or(false), false);
    TS_ASSERT_EQUALS(sm_or(true, true), true);
    TS_ASSERT_EQUALS(sm_or(true, false), true);
    TS_ASSERT_EQUALS(sm_or(false, true), true);
    TS_ASSERT_EQUALS(sm_or(false, false), false);
    TS_ASSERT_EQUALS(sm_or(false, true, false), true);
    TS_ASSERT_EQUALS(sm_or(true, false, true), true);
  }

  void TestNot()
  {
    TS_ASSERT_EQUALS(sm_not(true), false);
    TS_ASSERT_EQUALS(sm_not(false), true);
  }

  void TestXor()
  {
    TS_ASSERT_EQUALS(sm_xor(), false);
    TS_ASSERT_EQUALS(sm_xor(true), true);
    TS_ASSERT_EQUALS(sm_xor(false), false);
    TS_ASSERT_EQUALS(sm_xor(true, true), false);
    TS_ASSERT_EQUALS(sm_xor(true, false), true);
    TS_ASSERT_EQUALS(sm_xor(false, true), true);
    TS_ASSERT_EQUALS(sm_xor(false, false), false);
    TS_ASSERT_EQUALS(sm_xor(false, true, false), true);
    TS_ASSERT_EQUALS(sm_xor(true, false, true), false);
  }

  // Relational =================================
  void TestEq()
  {
    TS_ASSERT_EQUALS(sm_eq(1.0, 1.0), true);
    TS_ASSERT_EQUALS(sm_eq(1.0, 2.0), false);
    TS_ASSERT_EQUALS(sm_eq(2.0, 1.0), false);
  }

  void TestGeq()
  {
    TS_ASSERT_EQUALS(sm_geq(1.0, 1.0), true);
    TS_ASSERT_EQUALS(sm_geq(1.0, 2.0), false);
    TS_ASSERT_EQUALS(sm_geq(2.0, 1.0), true);
  }

  void TestGt()
  {
    TS_ASSERT_EQUALS(sm_gt(1.0, 1.0), false);
    TS_ASSERT_EQUALS(sm_gt(1.0, 2.0), false);
    TS_ASSERT_EQUALS(sm_gt(2.0, 1.0), true);
  }

  void TestLeq()
  {
    TS_ASSERT_EQUALS(sm_leq(1.0, 1.0), true);
    TS_ASSERT_EQUALS(sm_leq(1.0, 2.0), true);
    TS_ASSERT_EQUALS(sm_leq(2.0, 1.0), false);
  }

  void TestLt()
  {
    TS_ASSERT_EQUALS(sm_lt(1.0, 1.0), false);
    TS_ASSERT_EQUALS(sm_lt(1.0, 2.0), true);
    TS_ASSERT_EQUALS(sm_lt(2.0, 1.0), false);
  }

  void TestNeq()
  {
    TS_ASSERT_EQUALS(sm_neq(1.0, 1.0), false);
    TS_ASSERT_EQUALS(sm_neq(1.0, 2.0), true);
    TS_ASSERT_EQUALS(sm_neq(2.0, 1.0), true);
  }

  // Trigonometry ===============================

  void TestArccot()
  {
    TS_ASSERT_DELTA(sm_arccot(0.0), M_PI / 2.0, 1e-6);
    TS_ASSERT_DELTA(sm_arccot(1.0), M_PI / 4.0, 1e-6);
  }

  void TestArccoth()
  {
    TS_ASSERT_DELTA(sm_arccoth(2.0), 0.549306, 1e-6);
    TS_ASSERT_DELTA(sm_arccoth(10.0), 0.100335, 1e-6);
  }

  void TestArccsc()
  {
    TS_ASSERT_DELTA(sm_arccsc(1.0), 1.570796, 1e-6);
    TS_ASSERT_DELTA(sm_arccsc(10.0), 0.100167, 1e-6);
  }

  void TestArccsch()
  {
    TS_ASSERT_DELTA(sm_arccsch(1.0), 0.881373, 1e-6);
    TS_ASSERT_DELTA(sm_arccsch(10.0), 0.099834, 1e-6);
  }

  void TestArcsec()
  {
    TS_ASSERT_DELTA(sm_arcsec(1.0), 0.0, 1e-6);
    TS_ASSERT_DELTA(sm_arcsec(10.0), 1.470629, 1e-6);
  }

  void TestArcsech()
  {
    TS_ASSERT_DELTA(sm_arcsech(0.1), 2.993223, 1e-6);
    TS_ASSERT_DELTA(sm_arcsech(1.0), 0.0, 1e-6);
  }

  void TestCot()
  {
    TS_ASSERT_DELTA(sm_cot(M_PI / 2.0), 0.0, 1e-6);
    TS_ASSERT_DELTA(sm_cot(M_PI / 4.0), 1.0, 1e-6);
  }

  void TestCoth()
  {
    TS_ASSERT_DELTA(sm_coth(1.0), 1.313035, 1e-6);
    TS_ASSERT_DELTA(sm_coth(10.0), 1.0, 1e-6);
  }

  void TestCsc()
  {
    TS_ASSERT_DELTA(sm_csc(M_PI / 2.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_csc(M_PI / 4.0), std::sqrt(2), 1e-6);
  }

  void TestCsch()
  {
    TS_ASSERT_DELTA(sm_csch(1.0), 0.850918, 1e-6);
    TS_ASSERT_DELTA(sm_csch(10.0), 0.000091, 1e-6);
  }

  void TestSec()
  {
    TS_ASSERT_DELTA(sm_sec(0.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_sec(M_PI), -1.0, 1e-6);
  }

  void TestSech()
  {
    TS_ASSERT_DELTA(sm_sech(0.0), 1.0, 1e-6);
    TS_ASSERT_DELTA(sm_sech(10.0), 0.000091, 1e-6);
  }

  // Other functions ============================
  
  void TestFactorial()
  {
    TS_ASSERT_EQUALS(sm_factorial(0.0), 1.0);
    TS_ASSERT_EQUALS(sm_factorial(1.0), 1.0);
    TS_ASSERT_EQUALS(sm_factorial(2.0), 2.0);
    TS_ASSERT_EQUALS(sm_factorial(3.0), 6.0);
  }

  void TestMax()
  {
    TS_ASSERT_EQUALS(sm_max(1.0, 2.0), 2.0);
    TS_ASSERT_EQUALS(sm_max(2.0, 1.0), 2.0);
  }

  void TestMin()
  {
    TS_ASSERT_EQUALS(sm_min(1.0, 2.0), 1.0);
    TS_ASSERT_EQUALS(sm_min(2.0, 1.0), 1.0);
  }

  void TestPiecewise()
  {
    TS_ASSERT_EQUALS(sm_piecewise(1.0, true, 2.0), 1.0);
    TS_ASSERT_EQUALS(sm_piecewise(1.0, false, 2.0), 2.0);
    TS_ASSERT_EQUALS(sm_piecewise(1.0, false, 2.0, true, 3.0), 2.0);
    TS_ASSERT_EQUALS(sm_piecewise(1.0, false, 2.0, false, 3.0), 3.0);
  }

  void TestQuotient()
  {
    TS_ASSERT_EQUALS(sm_quotient(1.0, 1.0), 1.0);
    TS_ASSERT_EQUALS(sm_quotient(1.0, 2.0), 0.0);
    TS_ASSERT_EQUALS(sm_quotient(2.0, 1.0), 2.0);
  }
};

#endif // TESTSBMLMATH_HPP_
