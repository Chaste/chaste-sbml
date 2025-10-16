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

#ifndef TEST_CHEN_2004_SBML_ODE_SYSTEM_HPP_
#define TEST_CHEN_2004_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "SbmlTestHelpers.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "Timer.hpp"

#include "Chen2004SbmlOdeSystem.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace sth = sbmltesthelpers;

class TestChen2004SbmlOdeSystem : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 36u;
    const unsigned NUM_DERIVED_QUANTITIES = 35u;

    std::vector<double> default_initial_conditions = {
        0.008473,  // BUD
        0.238404,  // C2
        0.024034,  // C2P
        0.070081,  // C5
        0.006878,  // C5P
        0.468344,  // CDC14
        0.656533,  // CDC15
        0.444296,  // CDC20
        1.472044,  // CDC20i
        0.10758,   // CDC6
        0.015486,  // CDC6P
        0.930499,  // CDH1
        0.0695,    // CDH1i
        0.1469227, // CLB2
        0.0518014, // CLB5
        0.0652511, // CLN2
        0.301313,  // ESP1
        0.236058,  // F2
        0.0273938, // F2P
        7.24e-05,  // F5
        7.91e-05,  // F5P
        0.1015,    // IEP
        1.206019,  // MASS
        0.018645,  // NET1
        0.970271,  // NET1P
        0.000909,  // ORI
        0.025612,  // PDS1
        0.123179,  // PPX
        1.04954,   // RENT
        0.6,       // RENTP
        0.0228776, // SIC1
        0.00641,   // SIC1P
        0.03,      // SPN
        0.95,      // SWI5
        0.02,      // SWI5P
        0.9,       // TEM1GTP
    };

    void RunOdeWithSolver(AbstractIvpOdeSolver& rSolver, const std::string solverName)
    {
        try
        {
            // Solve system using solver
            Chen2004SbmlOdeSystem ode_system;

            const double max_step = 0.01;
            const double sampling_interval = 0.01;

            const double run_length = 100.0;
            double start_time = 0.0;
            double end_time = start_time + run_length;

            std::vector<double> initial_conditions;
            OdeSolution ode_solution;

            std::vector<double> times;
            std::vector<std::vector<double> > solutions;

            std::vector<double> expected_stop_times = {
                35.85,  // DNA synthesis
                84.01,  // Spindle checkpoint
                101.2,  // Cell division
                101.71, // Reset ORI
                136.81, // DNA synthesis
                184.87, // Spindle checkpoint
                202.14, // Cell division
                202.65  // Reset ORI
            };

            initial_conditions = ode_system.GetInitialConditions();

            // Run ODE until it stops, then start again with updated initial conditions
            for (unsigned i = 0; i < 8; i++)
            {
                Timer::Reset();
                ode_solution = rSolver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
                Timer::Print(std::to_string(i) + ". Chen 2004 (" + solverName + ")");

                // ODE should have stopped
                TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), true);

                // Check stopping times
                TS_ASSERT_DELTA(ode_solution.rGetTimes().back(), expected_stop_times[i], 1e-2);

                // Collate solutions and times from all runs
                solutions.insert(solutions.end(), ode_solution.rGetSolutions().begin(), ode_solution.rGetSolutions().end());
                times.insert(times.end(), ode_solution.rGetTimes().begin(), ode_solution.rGetTimes().end());

                // Update initial conditions and time for next run
                start_time = ode_solution.rGetTimes().back();
                end_time = start_time + run_length;

                ode_system.SetStateVariables(ode_solution.rGetSolutions().back());
                ode_system.AdjustParameters(start_time);
                initial_conditions = ode_system.GetStateVariables();
            }

            // Compare end solution with Tellurium values
            std::vector<double> end_solution = solutions.back();

            std::vector<std::vector<double> > expected_solution = {
                { 0.008181, 1e-3 },  // BUD
                { 0.237895, 1e-2 },  // C2
                { 0.024148, 1e-3 },  // C2P
                { 0.069449, 1e-3 },  // C5
                { 0.006869, 1e-4 },  // C5P
                { 0.471489, 1e-2 },  // CDC14
                { 0.656731, 1e-2 },  // CDC15
                { 0.446186, 1e-2 },  // CDC20
                { 1.473426, 1e-2 },  // CDC20i
                { 0.106310, 1e-2 },  // CDC6
                { 0.015478, 1e-3 },  // CDC6P
                { 0.931438, 1e-2 },  // CDH1
                { 0.068561, 1e-2 },  // CDH1i
                { 0.149486, 1e-2 },  // CLB2
                { 0.052216, 1e-3 },  // CLB5
                { 0.060328, 1e-3 },  // CLN2
                { 0.301714, 1e-2 },  // ESP1
                { 0.235779, 1e-2 },  // F2
                { 0.027547, 1e-3 },  // F2P
                { 0.000072, 1e-6 },  // F5
                { 0.000008, 1e-6 },  // F5P
                { 0.103716, 1e-2 },  // IEP
                { 1.207814, 1e-2 },  // MASS
                { 0.018455, 1e-3 },  // NET1
                { 1.252664, 1e-2 },  // NET1P
                { 21.662612, 1e-1 }, // ORI
                { 0.025638, 1e-3 },  // PDS1
                { 0.122829, 1e-2 },  // PPX
                { 1.043842, 1e-2 },  // RENT
                { 0.484669, 1e-2 },  // RENTP
                { 0.022570, 1e-3 },  // SIC1
                { 0.006414, 1e-4 },  // SIC1P
                { 0.030010, 1e-3 },  // SPN
                { 0.955855, 1e-2 },  // SWI5
                { 0.020352, 1e-3 },  // SWI5P
                { 0.905342, 1e-2 }   // TEM1GTP
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
            std::vector<std::vector<std::vector<double> > > expected_stats = {
                // {min, max, mean, stddev, q1, q2, q3}, {min_tol, max_tol, ...}
                { { 0.000000, 2.993992, 1.633851, 0.982415, 0.724748, 1.813629, 2.522914 }, { 1e-6, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // BUD
                { { 0.000160, 0.247435, 0.021261, 0.045354, 0.000882, 0.003010, 0.015465 }, { 1e-5, 1e-2, 1e-3, 1e-3, 1e-5, 1e-4, 1e-3 } },   // C2
                { { 0.000209, 0.074723, 0.017727, 0.025953, 0.000445, 0.001796, 0.036196 }, { 1e-5, 1e-3, 1e-3, 1e-3, 1e-5, 1e-4, 1e-3 } },   // C2P
                { { 0.004284, 0.186872, 0.054194, 0.064878, 0.006780, 0.009087, 0.112282 }, { 1e-4, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4, 1e-2 } },   // C5
                { { 0.001735, 0.061458, 0.020600, 0.012444, 0.013534, 0.018479, 0.024690 }, { 1e-4, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3 } },   // C5P
                { { 0.015026, 0.521218, 0.075377, 0.121140, 0.015455, 0.023294, 0.057560 }, { 1e-3, 1e-2, 1e-3, 1e-2, 1e-3, 1e-3, 1e-3 } },   // CDC14
                { { 0.023409, 0.661787, 0.209909, 0.239566, 0.023411, 0.143793, 0.341942 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-2, 1e-2 } },   // CDC15
                { { 0.000127, 0.516846, 0.099863, 0.173999, 0.000732, 0.006806, 0.085266 }, { 1e-5, 1e-2, 1e-3, 1e-2, 1e-5, 1e-4, 1e-3 } },   // CDC20
                { { 0.017442, 1.953515, 0.602146, 0.731832, 0.024657, 0.099796, 1.483997 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-2 } },   // CDC20i
                { { 0.004119, 1.240448, 0.310671, 0.433823, 0.008326, 0.020079, 0.642472 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-4, 1e-3, 1e-2 } },   // CDC6
                { { 0.008658, 0.100596, 0.045546, 0.030645, 0.022252, 0.030029, 0.075816 }, { 1e-4, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3 } },   // CDC6P
                { { 0.002462, 0.999609, 0.309029, 0.404586, 0.003032, 0.006550, 0.748230 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-4, 1e-4, 1e-2 } },   // CDH1
                { { 0.000390, 0.997538, 0.690971, 0.404586, 0.251770, 0.993450, 0.996967 }, { 1e-5, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // CDH1i
                { { 0.000035, 1.432926, 0.329852, 0.494323, 0.000311, 0.031457, 0.602707 }, { 1e-6, 1e-2, 1e-2, 1e-2, 1e-5, 1e-3, 1e-2 } },   // CLB2
                { { 0.000380, 0.474383, 0.233235, 0.179519, 0.004696, 0.275632, 0.412377 }, { 1e-5, 1e-2, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2 } },   // CLB5
                { { 0.056928, 2.387922, 1.309544, 0.788405, 0.495485, 1.532883, 2.004161 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // CLN2
                { { 0.010685, 0.575134, 0.231903, 0.211488, 0.028825, 0.152143, 0.455184 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-2, 1e-2 } },   // ESP1
                { { 0.000674, 0.236174, 0.018855, 0.038493, 0.003153, 0.005082, 0.018494 }, { 1e-5, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4, 1e-3 } },   // F2
                { { 0.000034, 0.087091, 0.024470, 0.031469, 0.000279, 0.006306, 0.054651 }, { 1e-6, 1e-3, 1e-3, 1e-3, 1e-5, 1e-4, 1e-3 } },   // F2P
                { { 0.000004, 0.001145, 0.000126, 0.000233, 0.000012, 0.000039, 0.000104 }, { 1e-6, 1e-4, 1e-5, 1e-5, 1e-6, 1e-6, 1e-5 } },   // F5
                { { 0.000002, 0.000498, 0.000064, 0.000104, 0.000011, 0.000026, 0.000060 }, { 1e-6, 1e-5, 1e-6, 1e-5, 1e-6, 1e-6, 1e-6 } },   // F5P
                { { 0.000002, 0.305134, 0.056806, 0.097790, 0.000028, 0.002304, 0.065103 }, { 1e-6, 1e-2, 1e-3, 1e-3, 1e-6, 1e-4, 1e-3 } },   // IEP
                { { 1.203171, 2.629232, 1.823300, 0.408710, 1.462739, 1.777846, 2.160835 }, { 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // MASS
                { { 0.015414, 0.755283, 0.477500, 0.267869, 0.217321, 0.498261, 0.733923 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // NET1
                { { 0.058685, 1.297002, 0.365543, 0.351461, 0.066714, 0.309558, 0.611109 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-2, 1e-2 } },   // NET1P
                { { 0.000909, 24.507328, 8.374791, 8.615211, 0.088982, 5.977344, 14.754822 }, { 1e-3, 1e-1, 1e-2, 1e-2, 1e-2, 1e-2, 1e-1 } }, // ORI
                { { 0.008394, 0.958757, 0.228468, 0.294696, 0.013110, 0.060313, 0.360250 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-2 } },   // PDS1
                { { 0.109931, 0.535635, 0.331474, 0.130887, 0.236150, 0.301166, 0.461647 }, { 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // PPX
                { { 0.869767, 1.958538, 1.779831, 0.269933, 1.691375, 1.891688, 1.957719 }, { 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // RENT
                { { 0.026285, 0.615347, 0.150612, 0.169743, 0.027015, 0.093227, 0.218653 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-2 } },   // RENTP
                { { 0.000857, 0.928910, 0.168594, 0.298425, 0.001049, 0.001389, 0.200006 }, { 1e-5, 1e-2, 1e-2, 1e-2, 1e-4, 1e-4, 1e-2 } },   // SIC1
                { { 0.002604, 0.073190, 0.017173, 0.023132, 0.002992, 0.004815, 0.016632 }, { 1e-4, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-3 } },   // SIC1P
                { { 0.000121, 1.302051, 0.367971, 0.449246, 0.024035, 0.095064, 0.712505 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-2 } },   // SPN
                { { 0.085173, 0.956892, 0.362550, 0.270644, 0.106128, 0.265494, 0.604890 }, { 1e-3, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2 } },   // SWI5
                { { 0.000005, 0.228048, 0.040593, 0.072724, 0.000033, 0.000842, 0.041815 }, { 1e-6, 1e-2, 1e-3, 1e-3, 1e-6, 1e-5, 1e-3 } },   // SWI5P
                { { 0.009990, 0.977834, 0.235812, 0.359780, 0.009990, 0.082109, 0.234611 }, { 1e-4, 1e-2, 1e-2, 1e-2, 1e-4, 1e-3, 1e-2 } },   // TEM1GTP
            };

            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                const char* var_name = var_names[i].c_str();

                std::vector<double> values;
                for (unsigned j = 0; j < times.size(); j++)
                {
                    values.push_back(solutions[j][i]);
                }

                double min_val = sth::Min(values);
                exp_val = expected_stats[i][0][0];
                tol = expected_stats[i][1][0];
                TSM_ASSERT_DELTA(var_name, min_val, exp_val, tol);

                double max_val = sth::Max(values);
                exp_val = expected_stats[i][0][1];
                tol = expected_stats[i][1][1];
                TSM_ASSERT_DELTA(var_name, max_val, exp_val, tol);

                double mean_val = sth::Mean(values);
                exp_val = expected_stats[i][0][2];
                tol = expected_stats[i][1][2];
                TSM_ASSERT_DELTA(var_name, mean_val, exp_val, tol);

                double std_val = sth::Stdev(values);
                exp_val = expected_stats[i][0][3];
                tol = expected_stats[i][1][3];
                TSM_ASSERT_DELTA(var_name, std_val, exp_val, tol);

                double q1_val = sth::Quantile(values, 0.25);
                exp_val = expected_stats[i][0][4];
                tol = expected_stats[i][1][4];
                TSM_ASSERT_DELTA(var_name, q1_val, exp_val, tol);

                double q2_val = sth::Quantile(values, 0.5);
                exp_val = expected_stats[i][0][5];
                tol = expected_stats[i][1][5];
                TSM_ASSERT_DELTA(var_name, q2_val, exp_val, tol);

                double q3_val = sth::Quantile(values, 0.75);
                exp_val = expected_stats[i][0][6];
                tol = expected_stats[i][1][6];
                TSM_ASSERT_DELTA(var_name, q3_val, exp_val, tol);
            }

            // Write results out to file
            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("chen_2004_" + solverName + ".csv");
            for (unsigned i = 0; i < times.size(); i++)
            {
                (*file) << times[i];
                for (unsigned j = 0; j < ODE_SIZE; j++)
                {
                    (*file) << ", " << solutions[i][j];
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
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chen_2004_ode.arch";
        std::vector<std::string> var_names;

        // Save archive
        {
            Chen2004SbmlOdeSystem ode_system;

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
        Chen2004SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> derivatives(ODE_SIZE);
        ode_system.EvaluateYDerivatives(time, default_initial_conditions, derivatives);

        // Compare derivatives with values from Tellurium
        std::vector<double> derivatives_expected = {
            0.013784,  // BUD **
            0.038071,  // C2
            -0.021920, // C2P
            0.045274,  // C5
            -0.003342, // C5P
            0.251316,  // CDC14
            -0.018917, // CDC15
            -0.034247, // CDC20
            -0.253250, // CDC20i
            0.123568,  // CDC6
            -0.003733, // CDC6P
            0.219318,  // CDH1
            -0.219318, // CDH1i
            -0.264860, // CLB2
            -0.051354, // CLB5
            -0.006941, // CLN2
            0.058559,  // ESP1
            0.003204,  // F2
            -0.027943, // F2P
            0.000178,  // F5
            -0.000214, // F5P
            -0.062337, // IEP
            0.009288,  // MASS
            -0.110538, // NET1
            0.378489,  // NET1P
            0.225418,  // ORI
            -0.014060, // PDS1
            -0.010198, // PPX
            0.300929,  // RENT
            -0.564034, // RENTP
            0.032294,  // SIC1
            -0.001699, // SIC1P
            0.049406,  // SPN
            -0.021724, // SWI5
            -0.013355, // SWI5P
            -0.130000, // TEM1GTP
        };
        // ** getRatesOfChange() shows 0.0142 in Tellurium for BUD.
        // This is slightly odd, but the difference in results is negligible.
        // Also, Chaste results for BUD are slightly closer to Matlab's results.

        TS_ASSERT_EQUALS(ode_system.GetNumberOfStateVariables(), ODE_SIZE);

        std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            TSM_ASSERT_DELTA(var_names[i].c_str(), derivatives[i], derivatives_expected[i], 1e-6);
        }

        // Check derived quantity indices
        TS_ASSERT_EQUALS(ode_system.GetNumberOfDerivedQuantities(), NUM_DERIVED_QUANTITIES);

        std::vector<std::string> dq_names = {
            "BCK2", "CDC14T", "CDC15i", "CDC6T", "CKIT", "CLB2T", "CLB5T",
            "CLN3", "IE", "MCM1", "NET1T", "PE", "SBF", "SIC1T", "TEM1GDP",
            "D", "mu", "Vdb5", "Vdb2", "Vasbf", "Visbf", "Vkpc1", "Vkpf6",
            "Vacdh", "Vicdh", "Vppnet", "Vkpnet", "Vdppx", "Vdpds",
            "Vaiep", "Vd2c1", "Vd2f6", "Vppc1", "Vppf6", "F"
        };

        for (unsigned i = 0; i < ode_system.GetNumberOfDerivedQuantities(); i++)
        {
            TSM_ASSERT_EQUALS(dq_names[i].c_str(), ode_system.GetDerivedQuantityIndex(dq_names[i]), i);
        }

        // Compare derived quantities with Tellurium values
        std::vector<double> dqs = ode_system.ComputeDerivedQuantities(0.0, default_initial_conditions);
        TS_ASSERT_EQUALS(dqs.size(), ode_system.GetNumberOfDerivedQuantities());

        std::vector<double> dqs_expected = {
            0.065125026,          // BCK2
            2.117884,             // CDC14T
            0.34346699999999997,  // CDC15i
            0.3866693,            // CDC6T
            0.7553539,            // CKIT
            0.6728125,            // CLB2T
            0.1289119,            // CLB5T
            0.06694509131879892,  // CLN3
            0.8985,               // IE
            0.4690076182110798,   // MCM1
            2.638456,             // NET1T
            0.6986870000000001,   // PE
            0.004913738759040801, // SBF
            0.36868460000000003,  // SIC1T
            0.09999999999999998,  // TEM1GDP
            101.21846007568689,   // D
            0.007701635339554948, // mu
            0.08108736,           // Vdb5
            0.441844,             // Vdb2
            0.5908263458114359,   // Vasbf
            1.7753815999999998,   // Visbf
            0.24226958712032795,  // Vkpc1
            0.2663495285233921,   // Vkpf6
            0.3846752,            // Vacdh
            0.05168441222637597,  // Vicdh
            0.419537,             // Vppnet
            0.4871349532762,      // Vkpnet
            0.8946186673169543,   // Vdppx
            0.13607916,           // Vdpds
            0.01469227,           // Vaiep
            0.09724769917563966,  // Vd2c1
            0.11193996917563967,  // Vd2f6
            1.873376,             // Vppc1
            1.873376,             // Vppf6
            0.4586134093959288,   // F
        };

        for (unsigned i = 0; i < dqs.size(); i++)
        {
            TSM_ASSERT_DELTA(dq_names[i].c_str(), dqs[i], dqs_expected[i], 1e-3);
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

#endif // TEST_CHEN_2004_SBML_ODE_SYSTEM_HPP_
