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
        // Solve system using backward Euler solver
        TysonNovak2001SbmlOdeSystem ode_system;

        BackwardEulerIvpOdeSolver backward_euler_solver(11);

        const double time_step = 0.01;
        const double sampling_interval = 0.01;

        const double run_length = 200.0;
        double start_time = 0.0;
        double end_time = start_time + run_length;

        std::vector<double> initial_conditions;
        std::vector<OdeSolution> ode_solutions;
        OdeSolution ode_solution;

        // Repeatedly run ODE until it stops, then start again with updated initial conditions
        for (unsigned i = 0; i < 5; i++)
        {
            initial_conditions = ode_system.GetInitialConditions();

            Timer::Reset();
            ode_solution = backward_euler_solver.Solve(&ode_system, initial_conditions, start_time, end_time, time_step, sampling_interval);
            Timer::Print(std::to_string(i) + ". Tyson Novak Backward Euler");

            ode_solutions.push_back(ode_solution);

            // ODE should have stopped
            TS_ASSERT_EQUALS(backward_euler_solver.StoppingEventOccurred(), true);

            // Update for next run
            start_time = ode_solution.rGetTimes().back();
            end_time = start_time + run_length;
        }

        // The following code provides nice output for gnuplot
        // use the command
        // plot "tysonnovak_2001.dat" u 1:2 etc. for the various species...
        // or
        // plot "tysonnovak_2001.dat" u 1:2, "" u 1:3, "" u 1:4 etc. for all species

        // std::vector<std::vector<double>> solutions;
        // for (const auto &ode_solution : ode_solutions)
        // {
        //     solutions.insert(solutions.end(), ode_solution.rGetSolutions().begin(), ode_solution.rGetSolutions().end());
        // }

        // std::vector<double> times;
        // for (const auto &ode_solution : ode_solutions)
        // {
        //     times.insert(times.end(), ode_solution.rGetTimes().begin(), ode_solution.rGetTimes().end());
        // }

        // OutputFileHandler handler("");
        // out_stream file = handler.OpenOutputFile("tysonnovak_2001.dat");
        // for (unsigned i = 0; i < solutions.size(); i++)
        // {
        //     (*file) << times[i];
        //     for (unsigned j = 0; j < solutions[i].size(); j++)
        //     {
        //         (*file) << "\t" << solutions[i][j];
        //     }
        //     (*file) << "\n"
        //             << std::flush;
        // }
        // file->close();

        // Values calculated using roadrunner
        std::vector<double> &solution = ode_solutions.back().rGetSolutions().back();
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
            // Solve system using CVODE solver
            TysonNovak2001SbmlOdeSystem ode_system;

            CvodeAdaptor cvode_solver;

            const double max_step = 0.01;
            const double sampling_interval = 0.01;

            const double run_length = 200.0;
            double start_time = 0.0;
            double end_time = start_time + run_length;

            std::vector<double> initial_conditions;
            std::vector<OdeSolution> ode_solutions;
            OdeSolution ode_solution;

            // Repeatedly run ODE until it stops, then start again with updated initial conditions
            for (unsigned i = 0; i < 5; i++)
            {
                initial_conditions = ode_system.GetInitialConditions();

                Timer::Reset();
                ode_solution = cvode_solver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
                Timer::Print(std::to_string(i) + ". Tyson Novak CVODE");

                ode_solutions.push_back(ode_solution);

                // ODE should have stopped
                TS_ASSERT_EQUALS(cvode_solver.StoppingEventOccurred(), true);

                // Update for next run
                start_time = ode_solution.rGetTimes().back();
                end_time = start_time + run_length;
            }

            // The following code provides nice output for gnuplot
            // use the command
            // plot "tysonnovak_2001.dat" u 1:2 etc. for the various species...
            // or
            // plot "tysonnovak_2001.dat" u 1:2, "" u 1:3, "" u 1:4 etc. for all species

            std::vector<std::vector<double>> solutions;
            for (const auto &ode_solution : ode_solutions)
            {
                solutions.insert(solutions.end(), ode_solution.rGetSolutions().begin(), ode_solution.rGetSolutions().end());
            }

            std::vector<double> times;
            for (const auto &ode_solution : ode_solutions)
            {
                times.insert(times.end(), ode_solution.rGetTimes().begin(), ode_solution.rGetTimes().end());
            }

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("tysonnovak_2001.dat");
            for (unsigned i = 0; i < solutions.size(); i++)
            {
                (*file) << times[i];
                for (unsigned j = 0; j < solutions[i].size(); j++)
                {
                    (*file) << "\t" << solutions[i][j];
                }
                (*file) << "\n"
                        << std::flush;
            }
            file->close();

            // Values calculated using roadrunner
            std::vector<double> &solution = ode_solutions.back().rGetSolutions().back();
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
