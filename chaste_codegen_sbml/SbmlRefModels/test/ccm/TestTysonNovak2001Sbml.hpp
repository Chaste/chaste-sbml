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

#ifndef TESTTYSONNOVAK2001SBML_HPP_
#define TESTTYSONNOVAK2001SBML_HPP_

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"
#include "Timer.hpp"

#include "TysonNovak2001SbmlOdeSystemAndCellCycleModel.hpp"
#include <boost/serialization/export.hpp>

// #include "PetscTools.hpp"
// #include "PetscSetupAndFinalize.hpp"
#include "FakePetscSetup.hpp"

class TestTysonNovak2001Sbml : public AbstractCellBasedTestSuite
{
public:
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak_2001_ode.arch";

        {
            std::vector<double> state_variables = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};

            TysonNovak2001SbmlOdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 11u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.001, 1e-6);  // CycBt
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6);    // CycB
            TS_ASSERT_DELTA(initial_conditions[2], 3.25, 1e-6);   // Cdc20a
            TS_ASSERT_DELTA(initial_conditions[3], 0.0, 1e-6);    // Trimer
            TS_ASSERT_DELTA(initial_conditions[4], 0.001, 1e-6);  // Cdh1
            TS_ASSERT_DELTA(initial_conditions[5], 0.5, 1e-6);    // m
            TS_ASSERT_DELTA(initial_conditions[6], 0.001, 1e-6);  // Cdc20t
            TS_ASSERT_DELTA(initial_conditions[7], 0.001, 1e-6);  // IEP
            TS_ASSERT_DELTA(initial_conditions[8], 0.0, 1e-6);    // Mad
            TS_ASSERT_DELTA(initial_conditions[9], 0.001, 1e-6);  // CKIt
            TS_ASSERT_DELTA(initial_conditions[10], 0.001, 1e-6); // SK

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem *const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem *p_ode_system = nullptr;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 11u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.001, 1e-6);  // CycBt
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6);    // CycB
            TS_ASSERT_DELTA(initial_conditions[2], 0.001, 1e-6);  // Cdc20a
            TS_ASSERT_DELTA(initial_conditions[3], 0.0, 1e-6);    // Trimer
            TS_ASSERT_DELTA(initial_conditions[4], 0.001, 1e-6);  // Cdh1
            TS_ASSERT_DELTA(initial_conditions[5], 0.5, 1e-6);    // m
            TS_ASSERT_DELTA(initial_conditions[6], 0.001, 1e-6);  // Cdc20t
            TS_ASSERT_DELTA(initial_conditions[7], 0.001, 1e-6);  // IEP
            TS_ASSERT_DELTA(initial_conditions[8], 0.0, 1e-6);    // Mad
            TS_ASSERT_DELTA(initial_conditions[9], 0.001, 1e-6);  // CKIt
            TS_ASSERT_DELTA(initial_conditions[10], 0.001, 1e-6); // SK

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

            TS_ASSERT_DELTA(var1, 1.0, 1e-6);
            TS_ASSERT_DELTA(var2, 2.0, 1e-6);
            TS_ASSERT_DELTA(var3, 3.0, 1e-6);
            TS_ASSERT_DELTA(var4, 4.0, 1e-6);
            TS_ASSERT_DELTA(var5, 5.0, 1e-6);
            TS_ASSERT_DELTA(var6, 6.0, 1e-6);
            TS_ASSERT_DELTA(var7, 7.0, 1e-6);
            TS_ASSERT_DELTA(var8, 8.0, 1e-6);
            TS_ASSERT_DELTA(var9, 9.0, 1e-6);
            TS_ASSERT_DELTA(var10, 10.0, 1e-6);
            TS_ASSERT_DELTA(var11, 11.0, 1e-6);

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        TysonNovak2001SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.001);
        initial_conditions.push_back(0.0);
        initial_conditions.push_back(0.001);
        initial_conditions.push_back(0.0);
        initial_conditions.push_back(0.001);
        initial_conditions.push_back(0.5);
        initial_conditions.push_back(0.001);
        initial_conditions.push_back(0.001);
        initial_conditions.push_back(0.0);
        initial_conditions.push_back(0.001);
        initial_conditions.push_back(0.001);

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct
        TS_ASSERT_DELTA(derivs[0], 0.0399, 1e-4);
        TS_ASSERT_DELTA(derivs[1], 0.0009, 1e-4);
        TS_ASSERT_DELTA(derivs[2], -0.0001, 1e-4);
        TS_ASSERT_DELTA(derivs[3], 0.0, 1e-4);
        TS_ASSERT_DELTA(derivs[4], 0.9710, 1e-4);
        TS_ASSERT_DELTA(derivs[5], 0.0023, 1e-4);
        TS_ASSERT_DELTA(derivs[6], 0.0049, 1e-4);
        TS_ASSERT_DELTA(derivs[7], -0.0, 1e-4);
        TS_ASSERT_DELTA(derivs[8], 1.0, 1e-4);
        TS_ASSERT_DELTA(derivs[9], 0.9997, 1e-4);
        TS_ASSERT_DELTA(derivs[10], 0.0278, 1e-4);
    }

    void TestOdeWithChasteSolver()
    {
        TysonNovak2001SbmlOdeSystem ode_system;

        // Solve system using backward Euler solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        const double timescale = 60.0;
        double dt = 0.1 / timescale;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(11);

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        Timer::Reset();
        OdeSolution solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 500.0, 0.01, 0.01);
        // OdeSolution solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 75.8350 / timescale, dt, dt);
        Timer::Print("1. Tyson Novak Backward Euler");

        // If you run it up to about 75min the ODE will stop, anything less and it will not and this test will fail
        TS_ASSERT_EQUALS(backward_euler_solver.StoppingEventOccurred(), true);

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste solvers
        // and CVODE.
        std::vector<double>& solution = solutions.rGetSolutions().back();
        TS_ASSERT_DELTA(solution[0], 0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(solution[1], 0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(solution[2], 1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(solution[3], 1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(solution[4], 0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(solution[5], 0.95328206604519, 2e-2);
        TS_ASSERT_DELTA(solution[6], 0.95328206604519, 2e-2);
        TS_ASSERT_DELTA(solution[7], 0.95328206604519, 2e-2);
        TS_ASSERT_DELTA(solution[8], 0.95328206604519, 2e-2);
        TS_ASSERT_DELTA(solution[9], 0.95328206604519, 2e-2);
        TS_ASSERT_DELTA(solution[10], 0.95328206604519, 2e-2);
    }

    void TestOdeWithCvodeSolver()
    {
        try
        {
            TysonNovak2001SbmlOdeSystem ode_system;

            double end_time = 500.0;
            double h_value = 0.01;

            CvodeAdaptor solver;
            OdeSolution solutions;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            Timer::Reset();
            solutions = solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 0.1);
            Timer::Print("1. Tyson Novak CVODE");

            // The following code provides nice output for gnuplot
            // use the command
            // plot "tysonnovak_2001.dat" u 1:2 etc. for the various species...
            // or
            // plot "tysonnovak_2001.dat" u 1:2, "" u 1:3, "" u 1:4 etc. for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("tysonnovak_2001.dat");
            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i++)
            {
                (*file) << solutions.rGetTimes()[i];
                for (unsigned j = 0; j < solutions.rGetSolutions()[i].size(); j++)
                {
                    (*file) << "\t" << solutions.rGetSolutions()[i][j];
                }
                (*file) << "\n"
                        << std::flush;
            }
            file->close();

            std::vector<double>& solution = solutions.rGetSolutions().back();
            TS_ASSERT_DELTA(solution[0], 0.10000000000000, 1e-2);
            TS_ASSERT_DELTA(solution[1], 0.98913684535843, 1e-2);
            TS_ASSERT_DELTA(solution[2], 1.54216806705641, 1e-1);
            TS_ASSERT_DELTA(solution[3], 1.40562614481544, 1e-1);
            TS_ASSERT_DELTA(solution[4], 0.67083371879876, 1e-2);
            TS_ASSERT_DELTA(solution[5], 0.95328206604519, 2e-2);
            TS_ASSERT_DELTA(solution[6], 0.95328206604519, 2e-2);
            TS_ASSERT_DELTA(solution[7], 0.95328206604519, 2e-2);
            TS_ASSERT_DELTA(solution[8], 0.95328206604519, 2e-2);
            TS_ASSERT_DELTA(solution[9], 0.95328206604519, 2e-2);
            TS_ASSERT_DELTA(solution[10], 0.95328206604519, 2e-2);
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

#endif // TESTTYSONNOVAK2001SBML_HPP_
