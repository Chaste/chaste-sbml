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

#ifndef TESTVANLEEUWEN2007SBML_HPP_
#define TESTVANLEEUWEN2007SBML_HPP_

#include <iostream>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "Debug.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "Timer.hpp"

#include "VanLeeuwen2007OdeSystemAndSrnModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestVanLeeuwen2007SBML : public AbstractCellBasedTestSuite
{
public:
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "van_leeuwen_2007_ode.arch";

        {
            std::vector<double> state_variables;
            state_variables.push_back(3.0);
            state_variables.push_back(4.0);
            state_variables.push_back(5.0);
            state_variables.push_back(6.0);
            state_variables.push_back(7.0);
            state_variables.push_back(8.0);
            state_variables.push_back(9.0);
            state_variables.push_back(10.0);
            state_variables.push_back(11.0);
            state_variables.push_back(12.0);
            state_variables.push_back(13.0);
            state_variables.push_back(14.0);
            state_variables.push_back(15.0);
            state_variables.push_back(16.0);

            VanLeeuwen2007OdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 14u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.067, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[1], 0.67, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[2], 3.25, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[3], 0.45, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[4], 0.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[5], 10.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[6], 18.14, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[7], 25.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[8], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[9], 0.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[10], 0.48, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[11], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[12], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[13], 1.0, 1e-3);

            double var1 = ode_system.GetStateVariable(0);
            double var2 = ode_system.GetStateVariable(1);
            double var3 = ode_system.GetStateVariable(2);
            double var4 = ode_system.GetStateVariable(3);
            double var5 = ode_system.GetStateVariable(4);
            double var6 = ode_system.GetStateVariable(5);
            double var7 = ode_system.GetStateVariable(6);
            double var8 = ode_system.GetStateVariable(7);
            double var9 = ode_system.GetStateVariable(8);
            double var10 = ode_system.GetStateVariable(9);
            double var11 = ode_system.GetStateVariable(10);
            double var12 = ode_system.GetStateVariable(11);
            double var13 = ode_system.GetStateVariable(12);
            double var14 = ode_system.GetStateVariable(13);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);
            TS_ASSERT_DELTA(var4, 6.0, 1e-3);
            TS_ASSERT_DELTA(var5, 7.0, 1e-3);
            TS_ASSERT_DELTA(var6, 8.0, 1e-3);
            TS_ASSERT_DELTA(var7, 9.0, 1e-3);
            TS_ASSERT_DELTA(var8, 10.0, 1e-3);
            TS_ASSERT_DELTA(var9, 11.0, 1e-3);
            TS_ASSERT_DELTA(var10, 12.0, 1e-3);
            TS_ASSERT_DELTA(var11, 13.0, 1e-3);
            TS_ASSERT_DELTA(var12, 14.0, 1e-3);
            TS_ASSERT_DELTA(var13, 15.0, 1e-3);
            TS_ASSERT_DELTA(var14, 16.0, 1e-3);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem *const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem *p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 14u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.067, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[1], 0.67, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[2], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[3], 0.45, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[4], 0.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[5], 10.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[6], 18.14, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[7], 25.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[8], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[9], 0.0, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[10], 0.48, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[11], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[12], 2.54, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[13], 1.0, 1e-3);

            double var1 = p_ode_system->GetStateVariable(0);
            double var2 = p_ode_system->GetStateVariable(1);
            double var3 = p_ode_system->GetStateVariable(2);
            double var4 = p_ode_system->GetStateVariable(3);
            double var5 = p_ode_system->GetStateVariable(4);
            double var6 = p_ode_system->GetStateVariable(5);
            double var7 = p_ode_system->GetStateVariable(6);
            double var8 = p_ode_system->GetStateVariable(7);
            double var9 = p_ode_system->GetStateVariable(8);
            double var10 = p_ode_system->GetStateVariable(9);
            double var11 = p_ode_system->GetStateVariable(10);
            double var12 = p_ode_system->GetStateVariable(11);
            double var13 = p_ode_system->GetStateVariable(12);
            double var14 = p_ode_system->GetStateVariable(13);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);
            TS_ASSERT_DELTA(var4, 6.0, 1e-3);
            TS_ASSERT_DELTA(var5, 7.0, 1e-3);
            TS_ASSERT_DELTA(var6, 8.0, 1e-3);
            TS_ASSERT_DELTA(var7, 9.0, 1e-3);
            TS_ASSERT_DELTA(var8, 10.0, 1e-3);
            TS_ASSERT_DELTA(var9, 11.0, 1e-3);
            TS_ASSERT_DELTA(var10, 12.0, 1e-3);
            TS_ASSERT_DELTA(var11, 13.0, 1e-3);
            TS_ASSERT_DELTA(var12, 14.0, 1e-3);
            TS_ASSERT_DELTA(var13, 15.0, 1e-3);
            TS_ASSERT_DELTA(var14, 16.0, 1e-3);

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        VanLeeuwen2007OdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.067); // X
        initial_conditions.push_back(0.67);  // D
        initial_conditions.push_back(2.54);  // C_o
        initial_conditions.push_back(0.45);  // C_u
        initial_conditions.push_back(0.00);  // C_c
        initial_conditions.push_back(10.0);  // A
        initial_conditions.push_back(18.14); // C_A
        initial_conditions.push_back(25.0);  // T
        initial_conditions.push_back(2.54);  // C_oT
        initial_conditions.push_back(0.00);  // C_cT
        initial_conditions.push_back(0.48);  // Y
        initial_conditions.push_back(2.54);  // C_F
        initial_conditions.push_back(2.54);  // C_T
        initial_conditions.push_back(1.0);   // drag

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct (default wnt-level=0)
        TS_ASSERT_DELTA(derivs[0], -0.0500, 1e-5);
        TS_ASSERT_DELTA(derivs[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2], -1.1103, 1e-5);
        TS_ASSERT_DELTA(derivs[3], 0.0703, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5], -1.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 1.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[9], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[10], 0.00344, 1e-5);
        TS_ASSERT_DELTA(derivs[11], -1.1103, 1e-5);
        TS_ASSERT_DELTA(derivs[12], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[13], 0.0, 1e-5);

        // Change Wnt level and check still OK
        ode_system.SetParameter("wnt_level", 1);
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        TS_ASSERT_DELTA(derivs[0], -10.1, 1e-5);
        TS_ASSERT_DELTA(derivs[1], -6.7, 1e-5);
        TS_ASSERT_DELTA(derivs[2], -1.1103, 1e-5);
        TS_ASSERT_DELTA(derivs[3], 0.0703, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5], -1.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 1.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[9], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[10], 0.00344, 1e-5);
        TS_ASSERT_DELTA(derivs[11], -1.1103, 1e-5);
        TS_ASSERT_DELTA(derivs[12], 0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[13], 0.0, 1e-5);
    }

    void TestOdeWithChasteSolver()
    {
        try
        {
            VanLeeuwen2007OdeSystem ode_system;

            // Solve system using RK4 solver

            double dt = 0.0001;

            // RK4 solver solution worked out
            RungeKutta4IvpOdeSolver rk4_solver;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            Timer::Reset();
            OdeSolution solutions = rk4_solver.Solve(&ode_system, state_variables, 0.0, 10.0, dt, dt);
            Timer::Print("1. VanLeeuwen RK4");

            unsigned end = solutions.rGetSolutions().size() - 1;

            // The following code provides nice output for gnuplot
            // use the command
            // plot "VanLeeuwen.dat" u 1:2
            // or
            // plot "VanLeeuwen.dat" u 1:3 etc. for the various proteins...

            // OutputFileHandler handler("");
            // out_stream file=handler.OpenOutputFile("VanLeeuwen.dat");
            // for (unsigned i=0; i<=end; i++)
            // {
            //     (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\t" << solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] << "\t" << solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\t" << solutions.rGetSolutions()[i][6] << "\t" << solutions.rGetSolutions()[i][7] << "\t" << solutions.rGetSolutions()[i][8] << "\t" << solutions.rGetSolutions()[i][9] << "\t" << solutions.rGetSolutions()[i][10] << "\n" << std::flush;
            // }
            // file->close();

            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.067, 2e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.67, 1e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 2.54, 1e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 0.45, 2e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 0.0, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 10.0, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 18.14, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 25.0, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 2.54, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9], 0.0, 1e-5);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10], 0.48, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11], 2.54, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12], 2.54, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13], 1.0, 1e-3);
        }
        catch (Exception &e)
        {
            throw e;
        }
        catch (...)
        {
            exit(EXIT_FAILURE);
        }
    }

    void TestOdeWithCvodeSolver()
    {
        try
        {
            VanLeeuwen2007OdeSystem ode_system;

            double end_time = 1000;
            double h_value = 0.01;

            CvodeAdaptor solver;
            OdeSolution solutions;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            Timer::Reset();
            solutions = solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 0.1);
            Timer::Print("1. VanLeeuwen CVODE");

            unsigned end = solutions.rGetSolutions().size() - 1;

            // The following code provides nice output for gnuplot
            // use the command
            // plot "VanLeeuwen.dat" u 1:2
            // or
            // plot "VanLeeuwen.dat" u 1:3 etc. for the various proteins...

            // OutputFileHandler handler("");
            // out_stream file=handler.OpenOutputFile("VanLeeuwen.dat");
            // for (unsigned i=0; i<=end; i++)
            // {
            //     (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\t" << solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] << "\t" << solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\n" << std::flush;
            // }
            // file->close();

            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.067, 2e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.67, 1e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 2.54, 1e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 0.45, 2e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 0.0, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 10.0, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6], 18.14, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7], 25.0, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8], 2.54, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9], 0.0, 1e-5);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10], 0.48, 5e-2);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11], 2.54, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12], 2.54, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13], 1.0, 1e-3);
        }
        catch (Exception &e)
        {
            throw e;
        }
        catch (...)
        {
            exit(EXIT_FAILURE);
        }
    }
};

#endif // TESTVANLEEUWEN2007SBML_HPP_
