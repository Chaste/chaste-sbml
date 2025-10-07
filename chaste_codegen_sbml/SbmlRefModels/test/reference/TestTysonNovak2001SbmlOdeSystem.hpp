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

#ifndef TEST_TYSON_NOVAK_2001_SBML_ODE_SYSTEM_HPP_
#define TEST_TYSON_NOVAK_2001_SBML_ODE_SYSTEM_HPP_

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"
#include "SbmlTestHelpers.hpp"
#include "Timer.hpp"

#include "TysonNovak2001SbmlOdeSystem.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace sth = sbmltesthelpers;

class TestTysonNovak2001SbmlOdeSystem : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 8u;
    const unsigned NUM_DERIVED_QUANTITIES = 3u;

    std::vector<double> default_initial_conditions = {
        0.001, // CycBt
        0.001, // Cdc20a
        0.001, // Cdh1
        0.5,   // m
        0.001, // Cdc20t
        0.001, // IEP
        0.001, // CKIt
        0.001  // SK
    };

    void RunOdeWithSolver(AbstractIvpOdeSolver& rSolver, const std::string solverName)
    {
        try
        {
            // Solve system using solver
            TysonNovak2001SbmlOdeSystem ode_system;

            const double max_step = 0.01;
            const double sampling_interval = 0.01;

            const double run_length = 200.0;
            double start_time = 0.0;
            double end_time = start_time + run_length;

            std::vector<double> initial_conditions;
            OdeSolution ode_solution;

            std::vector<double> times;
            std::vector<std::vector<double> > solutions;
            std::vector<std::vector<double> > derived_quantities(NUM_DERIVED_QUANTITIES);

            std::vector<double> expected_stop_times = {
                103.80, // Cell division
                250.32  // Cell division
            };

            initial_conditions = ode_system.GetInitialConditions();

            // Run ODE until it stops, then start again with updated initial conditions
            for (unsigned i = 0; i < 2; i++)
            {
                Timer::Reset();
                ode_solution = rSolver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
                Timer::Print(std::to_string(i) + ". Tyson Novak 2001 (" + solverName + ")");

                // ODE should have stopped
                TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), true);

                // Check stopping times
                TS_ASSERT_DELTA(ode_solution.rGetTimes().back(), expected_stop_times[i], 1e-2);

                // Collate solutions, times and derived quantities from all runs
                solutions.insert(solutions.end(), ode_solution.rGetSolutions().begin(), ode_solution.rGetSolutions().end());
                times.insert(times.end(), ode_solution.rGetTimes().begin(), ode_solution.rGetTimes().end());

                ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);
                std::vector<double> dq_cycb = ode_solution.GetAnyVariable("CycB");
                derived_quantities[0].insert(derived_quantities[0].end(), dq_cycb.begin(), dq_cycb.end());

                std::vector<double> dq_trimer = ode_solution.GetAnyVariable("Trimer");
                derived_quantities[1].insert(derived_quantities[1].end(), dq_trimer.begin(), dq_trimer.end());

                std::vector<double> dq_mad = ode_solution.GetAnyVariable("Mad");
                derived_quantities[2].insert(derived_quantities[2].end(), dq_mad.begin(), dq_mad.end());

                // Update initial conditions and time for next run
                start_time = ode_solution.rGetTimes().back();
                end_time = start_time + run_length;

                ode_system.SetStateVariables(ode_solution.rGetSolutions().back());
                ode_system.AdjustParameters(start_time);
                initial_conditions = ode_system.GetStateVariables();
            }

            // Compare solution stats with Tellurium values
            std::vector<std::vector<std::vector<double> > > expected_stats = {
                // {min, max, mean, stddev, q1, q2, q3}, {min_tol, max_tol, ...}
                { { 0.001000, 0.649287, 0.167091, 0.207068, 0.038556, 0.038778, 0.278560 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-2 } }, // CycBt
                { { 0.000001, 0.316459, 0.022664, 0.066668, 0.000011, 0.000388, 0.001993 }, { 1e-4, 1e-2, 1e-3, 1e-3, 1e-6, 1e-5, 1e-4 } }, // Cdc20a
                { { 0.001000, 0.999652, 0.655525, 0.452098, 0.014797, 0.991079, 0.997064 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-3, 1e-2, 1e-2 } }, // Cdh1
                { { 0.404845, 0.812511, 0.610497, 0.109946, 0.523156, 0.606383, 0.701869 }, { 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } }, // m
                { { 0.001000, 1.511564, 0.314996, 0.454442, 0.049908, 0.053844, 0.382421 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-2 } }, // Cdc20t
                { { 0.000635, 0.555370, 0.194751, 0.180816, 0.005462, 0.142140, 0.334785 }, { 1e-5, 1e-2, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2 } }, // IEP
                { { 0.001000, 0.971928, 0.305766, 0.288546, 0.046320, 0.212166, 0.513310 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-3, 1e-2, 1e-2 } }, // CKIt
                { { 0.001000, 0.446822, 0.071050, 0.094264, 0.017666, 0.034057, 0.074092 }, { 1e-4, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3 } }, // SK
                { { 0.000033, 0.629050, 0.127871, 0.212146, 0.000081, 0.000225, 0.222942 }, { 1e-6, 1e-2, 1e-2, 1e-2, 1e-6, 1e-5, 1e-2 } }, // CycB
                { { 0.000382, 0.139442, 0.039220, 0.013253, 0.037456, 0.038483, 0.038573 }, { 1e-5, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3 } }, // Trimer
                { { 1.000000, 1.000000, 1.000000, 0.000000, 1.000000, 1.000000, 1.000000 }, { 1e-2, 1e-2, 1e-2, 1e-6, 1e-2, 1e-2, 1e-2 } }, // Mad
            };

            for (unsigned i = 0; i < ODE_SIZE + NUM_DERIVED_QUANTITIES; i++)
            {
                std::string var_name;
                std::vector<double> values;

                if (i < ODE_SIZE)
                {
                    var_name = ode_system.rGetStateVariableNames()[i];
                    for (unsigned j = 0; j < times.size(); j++)
                    {
                        values.push_back(solutions[j][i]);
                    }
                }
                else
                {
                    var_name = ode_system.rGetDerivedQuantityNames()[i - ODE_SIZE];
                    values = derived_quantities[i - ODE_SIZE];
                }

                const char* var_cname = var_name.c_str();
                double exp_val;
                double tol;

                double min_val = sth::min(values);
                exp_val = expected_stats[i][0][0];
                tol = expected_stats[i][1][0];
                TSM_ASSERT_DELTA(var_cname, min_val, exp_val, tol);

                double max_val = sth::max(values);
                exp_val = expected_stats[i][0][1];
                tol = expected_stats[i][1][1];
                TSM_ASSERT_DELTA(var_cname, max_val, exp_val, tol);

                double mean_val = sth::mean(values);
                exp_val = expected_stats[i][0][2];
                tol = expected_stats[i][1][2];
                TSM_ASSERT_DELTA(var_cname, mean_val, exp_val, tol);

                double std_val = sth::stdev(values);
                exp_val = expected_stats[i][0][3];
                tol = expected_stats[i][1][3];
                TSM_ASSERT_DELTA(var_cname, std_val, exp_val, tol);

                double q1_val = sth::quantile(values, 0.25);
                exp_val = expected_stats[i][0][4];
                tol = expected_stats[i][1][4];
                TSM_ASSERT_DELTA(var_cname, q1_val, exp_val, tol);

                double q2_val = sth::quantile(values, 0.5);
                exp_val = expected_stats[i][0][5];
                tol = expected_stats[i][1][5];
                TSM_ASSERT_DELTA(var_cname, q2_val, exp_val, tol);

                double q3_val = sth::quantile(values, 0.75);
                exp_val = expected_stats[i][0][6];
                tol = expected_stats[i][1][6];
                TSM_ASSERT_DELTA(var_cname, q3_val, exp_val, tol);
            }

            // Write results out to file
            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("tyson_novak_2001_" + solverName + ".dat");
            for (unsigned i = 0; i < times.size(); i++)
            {
                (*file) << times[i];
                for (unsigned j = 0; j < ODE_SIZE; j++)
                {
                    (*file) << ", " << solutions[i][j];
                }
                for (unsigned j = 0; j < NUM_DERIVED_QUANTITIES; j++)
                {
                    (*file) << ", " << derived_quantities[j][i];
                }
                (*file) << "\n"
                        << std::flush;
            }
            file->close();
        }
        catch (Exception& e)
        {
            throw e;
        }
        catch (...)
        {
            exit(EXIT_FAILURE);
        }
    }

public:
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak_2001_ode.arch";
        std::vector<std::string> var_names;

        // Save archive
        {
            TysonNovak2001SbmlOdeSystem ode_system;

            // Set state variables to 0...ODE_SIZE-1
            std::vector<double> state_variables;
            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                state_variables.push_back(static_cast<double>(i));
            }
            ode_system.SetStateVariables(state_variables);

            // Check initial conditions and state variables
            ode_system.SetDefaultInitialCondition(0, 3.141593);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);

            var_names = ode_system.rGetStateVariableNames();
            TSM_ASSERT_DELTA(var_names[0].c_str(), initial_conditions[0], 3.141593, 1e-6);
            for (unsigned i = 1; i < ODE_SIZE; i++)
            {
                TSM_ASSERT_DELTA(var_names[i].c_str(), initial_conditions[i], default_initial_conditions[i], 1e-6);
            }

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        // Load archive
        {
            AbstractOdeSystem* p_ode_system = nullptr;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);

            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                // Check initial conditions
                TSM_ASSERT_DELTA(var_names[i].c_str(), initial_conditions[i], default_initial_conditions[i], 1e-6);

                // Check state variables
                double var = p_ode_system->GetStateVariable(i);
                TSM_ASSERT_DELTA(var_names[i].c_str(), var, static_cast<double>(i), 1e-6);
            }

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        TysonNovak2001SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> derivatives(ODE_SIZE);
        ode_system.EvaluateYDerivatives(time, default_initial_conditions, derivatives);

        // Compare derivatives with values from Tellurium
        std::vector<double> derivatives_expected = {
            3.99580000e-02,  // CycBt
            -2.50100000e-01, // Cdc20a
            9.70803883e-01,  // Cdh1
            2.37500000e-03,  // m
            4.90000000e-03,  // Cdc20t
            1.08707977e-05,  // IEP
            9.99719098e-01,  // CKIt
            2.77174939e-02   // SK
        };

        std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            TSM_ASSERT_DELTA(var_names[i].c_str(), derivatives[i], derivatives_expected[i], 1e-6);
        }

        // Check derived quantities
        TS_ASSERT_EQUALS(ode_system.GetNumberOfDerivedQuantities(), NUM_DERIVED_QUANTITIES);
        TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("CycB"), 0u);
        TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("Trimer"), 1u);
        TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("Mad"), 2u);

        std::vector<double> derived_quantities = ode_system.ComputeDerivedQuantities(0.0, default_initial_conditions);
        TS_ASSERT_EQUALS(derived_quantities.size(), NUM_DERIVED_QUANTITIES);
        TS_ASSERT_DELTA(derived_quantities[0], 0.001, 1e-3); // CycB
        TS_ASSERT_DELTA(derived_quantities[1], 0.001, 1e-3); // Trimer
        TS_ASSERT_DELTA(derived_quantities[2], 1.0, 1e-3);   // Mad
    }

    void TestOdeWithChasteSolver()
    {
        // Solve system using backward Euler solver
        BackwardEulerIvpOdeSolver backward_euler_solver(ODE_SIZE);
        RunOdeWithSolver(backward_euler_solver, "backward_euler");
    }

    void TestOdeWithCvodeSolver()
    {
        // Solve system using CVODE solver
        CvodeAdaptor cvode_solver;
        cvode_solver.CheckForStoppingEvents();
        RunOdeWithSolver(cvode_solver, "cvode");
    }
};

#endif // TEST_TYSON_NOVAK_2001_SBML_ODE_SYSTEM_HPP_
