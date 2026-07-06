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

#ifndef TEST_SBML_MATH_HPP_
#define TEST_SBML_MATH_HPP_

#include <cmath>

#include <cxxtest/TestSuite.h>

#include "SbmlMath.hpp"

namespace sm = sbmlmath;

class TestSbmlMath : public CxxTest::TestSuite
{
public:
    // Arithmetic =================================

    void TestDivide()
    {
        TS_ASSERT_EQUALS(sm::divide(1.0, 1.0), 1.0);
        TS_ASSERT_EQUALS(sm::divide(1.0, 2.0), 0.5);
        TS_ASSERT_EQUALS(sm::divide(2.0, 1.0), 2.0);

        TS_ASSERT_EQUALS(sm::divide(1, 2), 0.5);
        TS_ASSERT_EQUALS(sm::divide(2, 1), 2);
    }

    void TestMinus()
    {
        TS_ASSERT_EQUALS(sm::minus(1.0, 1.0), 0.0);
        TS_ASSERT_EQUALS(sm::minus(1.0, 2.0), -1.0);
        TS_ASSERT_EQUALS(sm::minus(2.0, 1.0), 1.0);

        TS_ASSERT_EQUALS(sm::minus(1, 2), -1);
        TS_ASSERT_EQUALS(sm::minus(2, 1), 1);
    }

    void TestPlus()
    {
        TS_ASSERT_EQUALS(sm::plus(), 0.0);
        TS_ASSERT_EQUALS(sm::plus(1.0), 1.0);
        TS_ASSERT_EQUALS(sm::plus(1.0, 2.0), 3.0);
        TS_ASSERT_EQUALS(sm::plus(1.0, 2.0, 3.0), 6.0);

        TS_ASSERT_EQUALS(sm::plus(1, 2), 3);
        TS_ASSERT_EQUALS(sm::plus(2, 1), 3);
    }

    void TestTimes()
    {
        TS_ASSERT_EQUALS(sm::times(), 1.0);
        TS_ASSERT_EQUALS(sm::times(1.0), 1.0);
        TS_ASSERT_EQUALS(sm::times(2.0), 2.0);
        TS_ASSERT_EQUALS(sm::times(1.0, 2.0), 2.0);
        TS_ASSERT_EQUALS(sm::times(1.0, 2.0, 3.0), 6.0);

        TS_ASSERT_EQUALS(sm::times(1, 2), 2);
        TS_ASSERT_EQUALS(sm::times(2, 1), 2);
    }

    // Logs and exponents =========================

    void TestLog()
    {
        TS_ASSERT_DELTA(sm::log(1.0), 0.0, 1e-6);
        TS_ASSERT_DELTA(sm::log(10.0), M_LN10, 1e-6);
        TS_ASSERT_DELTA(sm::log(M_E), 1.0, 1e-6);

        TS_ASSERT_DELTA(sm::log(2.0, 2.0), 1.0, 1e-6);
        TS_ASSERT_DELTA(sm::log(10.0, 2.0), std::log10(2.0), 1e-6);
        TS_ASSERT_DELTA(sm::log(M_E, 2.0), M_LN2, 1e-6);
    }

    void TestRoot()
    {
        TS_ASSERT_DELTA(sm::root(1.0, 1.0), 1.0, 1e-6);
        TS_ASSERT_DELTA(sm::root(2.0, 4.0), 2.0, 1e-6);
        TS_ASSERT_DELTA(sm::root(3.0, 8.0), 2.0, 1e-6);
        TS_ASSERT_DELTA(sm::root(2.0, 9.0), 3.0, 1e-6);
        TS_ASSERT_DELTA(sm::root(3.0, 27.0), 3.0, 1e-6);
    }

    void TestSqr()
    {
        TS_ASSERT_DELTA(sm::sqr(1.0), 1.0, 1e-6);
        TS_ASSERT_DELTA(sm::sqr(2.0), 4.0, 1e-6);
        TS_ASSERT_DELTA(sm::sqr(3.0), 9.0, 1e-6);
    }

    // Logical ===================================

    void TestAnd()
    {
        TS_ASSERT_EQUALS(sm::and_(), true);
        TS_ASSERT_EQUALS(sm::and_(true), true);
        TS_ASSERT_EQUALS(sm::and_(false), false);
        TS_ASSERT_EQUALS(sm::and_(true, true), true);
        TS_ASSERT_EQUALS(sm::and_(true, false), false);
        TS_ASSERT_EQUALS(sm::and_(false, true), false);
        TS_ASSERT_EQUALS(sm::and_(false, false), false);
        TS_ASSERT_EQUALS(sm::and_(true, true, true), true);
        TS_ASSERT_EQUALS(sm::and_(false, true, true), false);
    };

    void TestOr()
    {
        TS_ASSERT_EQUALS(sm::or_(), false);
        TS_ASSERT_EQUALS(sm::or_(true), true);
        TS_ASSERT_EQUALS(sm::or_(false), false);
        TS_ASSERT_EQUALS(sm::or_(true, true), true);
        TS_ASSERT_EQUALS(sm::or_(true, false), true);
        TS_ASSERT_EQUALS(sm::or_(false, true), true);
        TS_ASSERT_EQUALS(sm::or_(false, false), false);
        TS_ASSERT_EQUALS(sm::or_(false, true, false), true);
        TS_ASSERT_EQUALS(sm::or_(true, false, true), true);
    }

    void TestNot()
    {
        TS_ASSERT_EQUALS(sm::not_(true), false);
        TS_ASSERT_EQUALS(sm::not_(false), true);
    }
    void TestImplies()
    {
        // implies(a, b) == (not a) or b
        TS_ASSERT_EQUALS(sm::implies(true, true), true);
        TS_ASSERT_EQUALS(sm::implies(true, false), false);
        TS_ASSERT_EQUALS(sm::implies(false, true), true);
        TS_ASSERT_EQUALS(sm::implies(false, false), true);
    }

    void TestXor()
    {
        TS_ASSERT_EQUALS(sm::xor_(), false);
        TS_ASSERT_EQUALS(sm::xor_(true), true);
        TS_ASSERT_EQUALS(sm::xor_(false), false);
        TS_ASSERT_EQUALS(sm::xor_(true, true), false);
        TS_ASSERT_EQUALS(sm::xor_(true, false), true);
        TS_ASSERT_EQUALS(sm::xor_(false, true), true);
        TS_ASSERT_EQUALS(sm::xor_(false, false), false);
        TS_ASSERT_EQUALS(sm::xor_(false, true, false), true);
        TS_ASSERT_EQUALS(sm::xor_(true, false, true), false);
    }

    // Relational =================================
    void TestEq()
    {
        TS_ASSERT_EQUALS(sm::eq(1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::eq(2.0, 2.0), true);
        TS_ASSERT_EQUALS(sm::eq(1.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::eq(2.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(2.0, 2.0, 2.0), true);
        TS_ASSERT_EQUALS(sm::geq(1.0, 2.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(1.0, 2.0, 1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 2.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 1.0, 2.0), false);

        TS_ASSERT_EQUALS(sm::eq(2, 2), true);
        TS_ASSERT_EQUALS(sm::eq(1, 2), false);
    }

    void TestGeq()
    {
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(1.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(3.0, 2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(2.0, 3.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::geq(3.0, 1.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::geq(1.0, 1.0, 1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(4.0, 3.0, 2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(3.0, 2.0, 2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::geq(3.0, 2.0, 1.0, 1.5), false);

        TS_ASSERT_EQUALS(sm::geq(2, 1), true);
        TS_ASSERT_EQUALS(sm::geq(1, 2), false);
    }

    void TestGt()
    {
        TS_ASSERT_EQUALS(sm::gt(1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::gt(2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::gt(1.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::gt(1.0, 1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::gt(3.0, 2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::gt(2.0, 3.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::gt(3.0, 1.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::gt(1.0, 1.0, 1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::gt(4.0, 3.0, 2.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::gt(3.0, 2.0, 2.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::gt(3.0, 2.0, 1.0, 1.5), false);

        TS_ASSERT_EQUALS(sm::gt(1, 1), false);
        TS_ASSERT_EQUALS(sm::gt(2, 1), true);
    }

    void TestLeq()
    {
        TS_ASSERT_EQUALS(sm::leq(1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::leq(1.0, 2.0), true);
        TS_ASSERT_EQUALS(sm::leq(2.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::leq(1.0, 1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::leq(1.0, 2.0, 3.0), true);
        TS_ASSERT_EQUALS(sm::leq(1.0, 3.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::leq(2.0, 1.0, 3.0), false);
        TS_ASSERT_EQUALS(sm::leq(1.0, 1.0, 1.0, 1.0), true);
        TS_ASSERT_EQUALS(sm::leq(1.0, 2.0, 3.0, 4.0), true);
        TS_ASSERT_EQUALS(sm::leq(1.0, 2.0, 2.0, 3.0), true);
        TS_ASSERT_EQUALS(sm::leq(1.0, 2.0, 3.0, 2.5), false);

        TS_ASSERT_EQUALS(sm::leq(1, 1), true);
        TS_ASSERT_EQUALS(sm::leq(2, 1), false);
    }

    void TestLt()
    {
        TS_ASSERT_EQUALS(sm::lt(1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::lt(1.0, 2.0), true);
        TS_ASSERT_EQUALS(sm::lt(2.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::lt(1.0, 1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::lt(1.0, 2.0, 3.0), true);
        TS_ASSERT_EQUALS(sm::lt(1.0, 3.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::lt(2.0, 1.0, 3.0), false);
        TS_ASSERT_EQUALS(sm::lt(1.0, 1.0, 1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::lt(1.0, 2.0, 3.0, 4.0), true);
        TS_ASSERT_EQUALS(sm::lt(1.0, 2.0, 2.0, 3.0), false);
        TS_ASSERT_EQUALS(sm::lt(1.0, 2.0, 3.0, 2.5), false);

        TS_ASSERT_EQUALS(sm::lt(1, 1), false);
        TS_ASSERT_EQUALS(sm::lt(1, 2), true);
    }

    void TestNeq()
    {
        TS_ASSERT_EQUALS(sm::neq(1.0, 1.0), false);
        TS_ASSERT_EQUALS(sm::neq(2.0, 2.0), false);
        TS_ASSERT_EQUALS(sm::neq(1.0, 2.0), true);
        TS_ASSERT_EQUALS(sm::neq(2.0, 1.0), true);

        TS_ASSERT_EQUALS(sm::neq(1, 1), false);
        TS_ASSERT_EQUALS(sm::neq(1, 2), true);
    }

    // Trigonometry ===============================

    void TestArccot()
    {
        TS_ASSERT_DELTA(sm::acot(0.0), M_PI / 2.0, 1e-6);
        TS_ASSERT_DELTA(sm::acot(1.0), M_PI / 4.0, 1e-6);
    }

    void TestArccoth()
    {
        TS_ASSERT_DELTA(sm::acoth(2.0), 0.549306, 1e-6);
        TS_ASSERT_DELTA(sm::acoth(10.0), 0.100335, 1e-6);
    }

    void TestArccsc()
    {
        TS_ASSERT_DELTA(sm::acsc(1.0), 1.570796, 1e-6);
        TS_ASSERT_DELTA(sm::acsc(10.0), 0.100167, 1e-6);
    }

    void TestArccsch()
    {
        TS_ASSERT_DELTA(sm::acsch(1.0), 0.881373, 1e-6);
        TS_ASSERT_DELTA(sm::acsch(10.0), 0.099834, 1e-6);
    }

    void TestArcsec()
    {
        TS_ASSERT_DELTA(sm::asec(1.0), 0.0, 1e-6);
        TS_ASSERT_DELTA(sm::asec(10.0), 1.470629, 1e-6);
    }

    void TestArcsech()
    {
        TS_ASSERT_DELTA(sm::asech(0.1), 2.993223, 1e-6);
        TS_ASSERT_DELTA(sm::asech(1.0), 0.0, 1e-6);
    }

    void TestCot()
    {
        TS_ASSERT_DELTA(sm::cot(M_PI / 2.0), 0.0, 1e-6);
        TS_ASSERT_DELTA(sm::cot(M_PI / 4.0), 1.0, 1e-6);
    }

    void TestCoth()
    {
        TS_ASSERT_DELTA(sm::coth(1.0), 1.313035, 1e-6);
        TS_ASSERT_DELTA(sm::coth(10.0), 1.0, 1e-6);
    }

    void TestCsc()
    {
        TS_ASSERT_DELTA(sm::csc(M_PI / 2.0), 1.0, 1e-6);
        TS_ASSERT_DELTA(sm::csc(M_PI / 4.0), std::sqrt(2), 1e-6);
    }

    void TestCsch()
    {
        TS_ASSERT_DELTA(sm::csch(1.0), 0.850918, 1e-6);
        TS_ASSERT_DELTA(sm::csch(10.0), 0.000091, 1e-6);
    }

    void TestSec()
    {
        TS_ASSERT_DELTA(sm::sec(0.0), 1.0, 1e-6);
        TS_ASSERT_DELTA(sm::sec(M_PI), -1.0, 1e-6);
    }

    void TestSech()
    {
        TS_ASSERT_DELTA(sm::sech(0.0), 1.0, 1e-6);
        TS_ASSERT_DELTA(sm::sech(10.0), 0.000091, 1e-6);
    }

    // Other functions ============================

    void TestFactorial()
    {
        TS_ASSERT_EQUALS(sm::factorial(0.0), 1.0);
        TS_ASSERT_EQUALS(sm::factorial(1.0), 1.0);
        TS_ASSERT_EQUALS(sm::factorial(2.0), 2.0);
        TS_ASSERT_EQUALS(sm::factorial(3.0), 6.0);
        TS_ASSERT_EQUALS(sm::factorial(10.0), 3628800.0);

        TS_ASSERT_EQUALS(sm::factorial(0), 1);
        TS_ASSERT_EQUALS(sm::factorial(1), 1);
    }

    void TestMax()
    {
        TS_ASSERT_EQUALS(sm::max(2.0), 2.0);  // single argument
        TS_ASSERT_EQUALS(sm::max(1.0, 2.0), 2.0);
        TS_ASSERT_EQUALS(sm::max(2.0, 1.0), 2.0);
        TS_ASSERT_EQUALS(sm::max(1.0, 2.0, 3.0), 3.0);
        TS_ASSERT_EQUALS(sm::max(1.0, 2.0, 3.0, 4.0), 4.0);
        TS_ASSERT_EQUALS(sm::max(1.0, 2.0, 3.0, 5.0, 4.0), 5.0);
        TS_ASSERT_EQUALS(sm::max(1.0, 2.0, 3.0, 6.0, 4.0, 5.0), 6.0);

        TS_ASSERT_EQUALS(sm::max(1, 2), 2);
        TS_ASSERT_EQUALS(sm::max(2, 1), 2);
    }

    void TestMin()
    {
        TS_ASSERT_EQUALS(sm::min(2.0), 2.0);  // single argument
        TS_ASSERT_EQUALS(sm::min(1.0, 2.0), 1.0);
        TS_ASSERT_EQUALS(sm::min(2.0, 1.0), 1.0);
        TS_ASSERT_EQUALS(sm::min(3.0, 2.0, 1.0), 1.0);
        TS_ASSERT_EQUALS(sm::min(4.0, 3.0, 2.0, 1.0), 1.0);
        TS_ASSERT_EQUALS(sm::min(4.0, 3.0, 2.0, 1.0, 5.0), 1.0);
        TS_ASSERT_EQUALS(sm::min(5.0, 4.0, 3.0, 1.0, 2.0, 6.0), 1.0);

        TS_ASSERT_EQUALS(sm::min(1, 2), 1);
        TS_ASSERT_EQUALS(sm::min(2, 1), 1);
    }

    void TestPiecewise()
    {
        TS_ASSERT_EQUALS(sm::piecewise(1.0, true, 2.0), 1.0);
        TS_ASSERT_EQUALS(sm::piecewise(1.0, false, 2.0), 2.0);
        TS_ASSERT_EQUALS(sm::piecewise(1.0, false, 2.0, true, 3.0), 2.0);
        TS_ASSERT_EQUALS(sm::piecewise(1.0, false, 2.0, false, 3.0), 3.0);

        TS_ASSERT_EQUALS(sm::piecewise(1, true, 2), 1);
        TS_ASSERT_EQUALS(sm::piecewise(1, false, 2), 2);
    }

    void TestQuotient()
    {
        TS_ASSERT_EQUALS(sm::quotient(1.0, 1.0), 1.0);
        TS_ASSERT_EQUALS(sm::quotient(1.0, 2.0), 0.0);
        TS_ASSERT_EQUALS(sm::quotient(2.0, 1.0), 2.0);

        TS_ASSERT_EQUALS(sm::quotient(1, 1), 1);
        TS_ASSERT_EQUALS(sm::quotient(1, 2), 0);
    }
};

#endif // TEST_SBML_MATH_HPP_
