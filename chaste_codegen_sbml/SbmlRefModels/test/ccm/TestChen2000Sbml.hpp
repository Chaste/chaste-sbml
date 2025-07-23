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

#ifndef TESTCHEN2000SBML_HPP_
#define TESTCHEN2000SBML_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/make_shared.hpp>
#include <boost/serialization/export.hpp>
#include <boost/shared_ptr.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "SbmlTestHelperFunctions.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "Timer.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Chen2000SbmlOdeSystemAndCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace st = sbmltest;

class TestChen2000Sbml : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 13u;

    std::vector<double> default_initial_conditions = {
        0.0078, // Cln2
        0.2342, // Clb2_T
        0.0614, // Clb5_T
        0.1231, // Sic1_T
        0.079,  // Clb2_Sic1
        0.0207, // Clb5_Sic1
        0.8332, // Cdc20_T
        0.6848, // Cdc20
        0.9946, // Hct1
        0.6608, // mass
        0.0,    // ORI
        0.0,    // BUD
        0.0     // SPN
    };

    void RunOdeWithSolver(AbstractIvpOdeSolver &rSolver, const std::string solverName)
    {
        try
        {
            // Solve system using solver
            Chen2000SbmlOdeSystem ode_system;

            const double max_step = 0.01;
            const double sampling_interval = 0.01;

            const double start_time = 0.0;
            const double end_time = 300.0;

            std::vector<double> initial_conditions;
            OdeSolution ode_solution;

            initial_conditions = ode_system.GetInitialConditions();

            Timer::Reset();
            ode_solution = rSolver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
            Timer::Print("Chen 2000 (" + solverName + ")");

            // Compare end solutions with Tellurium values
            std::vector<double> end_solution = ode_solution.rGetSolutions().back();

            std::vector<std::vector<double>> expected_solution = {

                {0.001184, 1e-4},  // Cln2
                {3.619189, 1e-2},  // Clb2_T
                {0.138153, 1e-2},  // Clb5_T
                {0.007854, 1e-4},  // Sic1_T
                {0.007456, 1e-4},  // Clb2_Sic1
                {0.000274, 1e-5},  // Clb5_Sic1
                {2.604753, 1e-2},  // Cdc20_T
                {0.234977, 1e-2},  // Cdc20
                {0.012937, 1e-3},  // Hct1
                {3.737796, 1e-2},  // mass
                {48.574775, 1e-2}, // ORI
                {0.687618, 1e-2},  // BUD
                {1.257020, 1e-2},  // SPN
            };

            std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
            double exp_val;
            double tol;

            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                exp_val = expected_solution[i][0];
                tol = expected_solution[i][1];
                TSM_ASSERT_DELTA(var_names[i].c_str(), end_solution[i], exp_val, tol);
            }

            // Compare solution stats with Tellurium values
            std::vector<double> times = ode_solution.rGetTimes();
            std::vector<std::vector<double>> solutions = ode_solution.rGetSolutions();

            std::vector<std::vector<std::vector<double>>> expected_stats = {
                // {min, max, mean, stddev, q1, q2, q3}, {min_tol, max_tol, ...}
                {{0.000742, 0.519493, 0.050706, 0.118540, 0.000939, 0.005803, 0.017823}, {1e-05, 1e-02, 1e-02, 1e-02, 1e-05, 1e-04, 1e-03}},      // Cln2
                {{0.000677, 3.619189, 1.368818, 1.228197, 0.001046, 1.509475, 2.410650}, {1e-05, 1e-02, 1e-02, 1e-02, 1e-04, 1e-02, 1e-02}},      // Clb2_T
                {{0.023832, 0.267919, 0.096199, 0.049997, 0.060660, 0.090080, 0.118051}, {1e-03, 1e-02, 1e-02, 1e-03, 1e-03, 1e-03, 1e-02}},      // Clb5_T
                {{0.004212, 1.316170, 0.321432, 0.490366, 0.010169, 0.014721, 0.797223}, {1e-04, 1e-02, 1e-02, 1e-02, 1e-03, 1e-03, 1e-02}},      // Sic1_T
                {{0.000076, 0.079031, 0.007810, 0.006538, 0.000921, 0.008694, 0.012629}, {1e-06, 1e-03, 1e-03, 1e-04, 1e-05, 1e-04, 1e-03}},      // Clb2_Sic1
                {{0.000274, 0.165560, 0.017047, 0.030014, 0.000444, 0.000802, 0.029744}, {1e-05, 1e-02, 1e-02, 1e-03, 1e-05, 1e-05, 1e-03}},      // Clb5_Sic1
                {{0.063078, 2.604753, 1.009452, 0.844044, 0.108376, 0.946726, 1.743825}, {1e-03, 1e-02, 1e-02, 1e-02, 1e-02, 1e-02, 1e-02}},      // Cdc20_T
                {{0.005739, 0.689086, 0.132325, 0.102236, 0.058253, 0.121006, 0.177074}, {1e-04, 1e-02, 1e-02, 1e-02, 1e-03, 1e-02, 1e-02}},      // Cdc20
                {{0.006011, 0.999781, 0.315141, 0.447497, 0.013060, 0.013735, 0.993248}, {1e-04, 1e-02, 1e-02, 1e-02, 1e-03, 1e-03, 1e-02}},      // Hct1
                {{0.660800, 3.737796, 1.775750, 0.867305, 1.019076, 1.571603, 2.423702}, {1e-06, 1e-02, 1e-02, 1e-02, 1e-02, 1e-02, 1e-02}},      // mass
                {{0.000000, 48.574775, 18.086541, 16.327771, 0.033456, 17.437463, 32.358107}, {1e-06, 1e-02, 1e-02, 1e-02, 1e-03, 1e-02, 1e-02}}, // ORI
                {{0.000000, 2.518914, 0.633002, 0.613297, 0.094799, 0.526040, 0.666904}, {1e-06, 1e-02, 1e-02, 1e-02, 1e-03, 1e-02, 1e-02}},      // BUD
                {{0.000000, 1.257020, 0.700241, 0.557746, 0.003674, 1.037460, 1.219701}, {1e-06, 1e-02, 1e-02, 1e-02, 1e-04, 1e-02, 1e-02}},      // SPN
            };

            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                const char *var_name = var_names[i].c_str();

                std::vector<double> values;
                for (unsigned j = 0; j < times.size(); j++)
                {
                    values.push_back(solutions[j][i]);
                }

                double min_val = st::min(values);
                exp_val = expected_stats[i][0][0];
                tol = expected_stats[i][1][0];
                TSM_ASSERT_DELTA(var_name, min_val, exp_val, tol);

                double max_val = st::max(values);
                exp_val = expected_stats[i][0][1];
                tol = expected_stats[i][1][1];
                TSM_ASSERT_DELTA(var_name, max_val, exp_val, tol);

                double mean_val = st::mean(values);
                exp_val = expected_stats[i][0][2];
                tol = expected_stats[i][1][2];
                TSM_ASSERT_DELTA(var_name, mean_val, exp_val, tol);

                double std_val = st::stdev(values);
                exp_val = expected_stats[i][0][3];
                tol = expected_stats[i][1][3];
                TSM_ASSERT_DELTA(var_name, std_val, exp_val, tol);

                double q1_val = st::quantile(values, 0.25);
                exp_val = expected_stats[i][0][4];
                tol = expected_stats[i][1][4];
                TSM_ASSERT_DELTA(var_name, q1_val, exp_val, tol);

                double q2_val = st::quantile(values, 0.5);
                exp_val = expected_stats[i][0][5];
                tol = expected_stats[i][1][5];
                TSM_ASSERT_DELTA(var_name, q2_val, exp_val, tol);

                double q3_val = st::quantile(values, 0.75);
                exp_val = expected_stats[i][0][6];
                tol = expected_stats[i][1][6];
                TSM_ASSERT_DELTA(var_name, q3_val, exp_val, tol);
            }

            // The following code provides nice output for gnuplot
            // use the command
            // plot "chen_2000_cvode.dat" u 1:2 etc. for the various species...
            // or
            // plot "chen_2000_cvode.dat" u 1:2, "" u 1:3, "" u 1:4 etc. for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("chen_2000_" + solverName + ".dat");
            for (unsigned i = 0; i < times.size(); i++)
            {
                (*file) << times[i];
                for (unsigned j = 0; j < ODE_SIZE; j++)
                {
                    (*file) << "\t" << solutions[i][j];
                }
                (*file) << "\n"
                        << std::flush;
            }
            file->close();
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

public:
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chen_2000_ode.arch";
        std::vector<std::string> var_names;

        // Save archive
        {
            // Set state variables to 0...ODE_SIZE-1
            std::vector<double> state_variables;
            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                state_variables.push_back(static_cast<double>(i));
            }

            // Check initial conditions
            Chen2000SbmlOdeSystem ode_system(state_variables);
            ode_system.SetDefaultInitialCondition(0, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);

            var_names = ode_system.rGetStateVariableNames();
            TSM_ASSERT_DELTA(var_names[0].c_str(), initial_conditions[0], 3.25, 1e-6);
            for (unsigned i = 1; i < ODE_SIZE; i++)
            {
                TSM_ASSERT_DELTA(var_names[i].c_str(), initial_conditions[i], default_initial_conditions[i], 1e-6);
            }

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem *const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        // Load archive
        {
            AbstractOdeSystem *p_ode_system = nullptr;

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
        Chen2000SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> derivatives(ODE_SIZE);
        ode_system.EvaluateYDerivatives(time, default_initial_conditions, derivatives);

        // Compare derivatives with values from Tellurium
        std::vector<double> derivatives_expected = {
            -0.0006895159121089688, // Cln2
            -0.4556385848470362,    // Clb2_T
            -0.012650686364843589,  // Clb5_T
            0.08924097232728236,    // Sic1_T
            0.0030021479399205686,  // Clb2_Sic1
            0.037098519500713356,   // Clb5_Sic1
            -0.05234400000000001,   // Cdc20_T
            0.025136000000000103,   // Cdc20
            0.02446145765480369,    // Hct1
            0.0038167808000000004,  // mass
            0.20556000000000002,    // ORI
            0.01514524381455681,    // BUD
            0.03495495495495495     // SPN
        };

        std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            TSM_ASSERT_DELTA(var_names[i].c_str(), derivatives[i], derivatives_expected[i], 1e-6);
        }
    }

    void TestOdeWithChasteSolver()
    {
        // Solve system using RK4 solver
        RungeKutta4IvpOdeSolver rk4_solver;
        RunOdeWithSolver(rk4_solver, "rk4");
    }

    void TestOdeWithCvodeSolver()
    {
        // Solve system using CVODE solver
        CvodeAdaptor cvode_solver;
        cvode_solver.CheckForStoppingEvents();
        RunOdeWithSolver(cvode_solver, "cvode");
    }
};

#endif // TESTCHEN2000SBML_HPP_
