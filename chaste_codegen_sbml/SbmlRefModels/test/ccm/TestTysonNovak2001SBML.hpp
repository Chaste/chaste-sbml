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

#include <iostream>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "PetscTools.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "Timer.hpp"

#include "TysonNovak2001OdeSystemAndCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestTysonNovak2001SBML : public AbstractCellBasedTestSuite
{
public:
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tn_ode.arch";

        {
            TysonNovak2001OdeSystem ode_system;

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 11u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.001, 1e-6); // CycBt
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6); // CycB
            TS_ASSERT_DELTA(initial_conditions[2], 3.25, 1e-6); // Cdc20a
            TS_ASSERT_DELTA(initial_conditions[3], 0.0, 1e-6); // Trimer
            TS_ASSERT_DELTA(initial_conditions[4], 0.001, 1e-6); // Cdh1
            TS_ASSERT_DELTA(initial_conditions[5], 0.5, 1e-6); // m
            TS_ASSERT_DELTA(initial_conditions[6], 0.001, 1e-6); // Cdc20t
            TS_ASSERT_DELTA(initial_conditions[7], 0.001, 1e-6); // IEP
            TS_ASSERT_DELTA(initial_conditions[8], 0.0, 1e-6); // Mad
            TS_ASSERT_DELTA(initial_conditions[9], 0.001, 1e-6); // CKIt
            TS_ASSERT_DELTA(initial_conditions[10], 0.001, 1e-6); // SK

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
            TS_ASSERT_EQUALS(initial_conditions.size(), 11u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.001, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 0.001, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[3], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[4], 0.001, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[5], 0.5, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[6], 0.001, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[7], 0.001, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[8], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[9], 0.001, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[10], 0.001, 1e-6);

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        TysonNovak2001OdeSystem ode_system;

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
        TS_ASSERT_DELTA(derivs[0], 0.0399, 1e-6);
        TS_ASSERT_DELTA(derivs[1], 0.1198, 1e-6);
        TS_ASSERT_DELTA(derivs[2], -0.0001, 1e-6);
        TS_ASSERT_DELTA(derivs[3], 0.0799, 1e-6);
        TS_ASSERT_DELTA(derivs[4], 0.9710, 1e-6);
        TS_ASSERT_DELTA(derivs[5], 0.0023, 1e-6);
        TS_ASSERT_DELTA(derivs[6], 0.0049, 1e-6);
        TS_ASSERT_DELTA(derivs[7], -0.0000, 1e-6);
        TS_ASSERT_DELTA(derivs[8], 0.9997, 1e-6);
        TS_ASSERT_DELTA(derivs[9], 0.0278, 1e-6);
        TS_ASSERT_DELTA(derivs[10], 0.0, 1e-6);
    }

    void TestOdeWithChasteSolver()
    {
        TysonNovak2001OdeSystem ode_system;

        // Solve system using backward Euler solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        double dt = 0.1 / 60.0;

        // Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(6);

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        Timer::Reset();
        OdeSolution solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 75.8350 / 60.0, dt, dt);
        Timer::Print("1. Tyson Novak Backward Euler");

        // If you run it up to about 75min the ODE will stop, anything less and it will not and this test will fail
        TS_ASSERT_EQUALS(backward_euler_solver.StoppingEventOccurred(), true);

        unsigned end = solutions.rGetSolutions().size() - 1;

        // The following code provides nice output for gnuplot
        // use the command
        // plot "tyson_novak.dat" u 1:2
        // or
        // plot "tyson_novak.dat" u 1:3 etc. for the various proteins...

        //        OutputFileHandler handler("");
        //        out_stream file=handler.OpenOutputFile("tyson_novak.dat");
        //        for (unsigned i=0; i<=end; i++)
        //        {
        //            (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\t" << solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] << "\t" << solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\n" << std::flush;
        //        }
        //        file->close();

        ColumnDataWriter writer("TysonNovak", "TysonNovak");
        if (PetscTools::AmMaster()) // if master process
        {
            int step_per_row = 1;
            int time_var_id = writer.DefineUnlimitedDimension("Time", "s");

            std::vector<int> var_ids;
            for (unsigned i = 0; i < ode_system.rGetStateVariableNames().size(); i++)
            {
                var_ids.push_back(writer.DefineVariable(ode_system.rGetStateVariableNames()[i],
                                                        ode_system.rGetStateVariableUnits()[i]));
            }
            writer.EndDefineMode();

            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i += step_per_row)
            {
                writer.PutVariable(time_var_id, solutions.rGetTimes()[i]);
                for (unsigned j = 0; j < var_ids.size(); j++)
                {
                    writer.PutVariable(var_ids[j], solutions.rGetSolutions()[i][j]);
                }
                writer.AdvanceAlongUnlimitedDimension();
            }
            writer.Close();
        }
        PetscTools::Barrier();

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste solvers
        // and CVODE.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 0.95328206604519, 2e-2);
    }
};

#endif // TESTTYSONNOVAK2001SBML_HPP_
