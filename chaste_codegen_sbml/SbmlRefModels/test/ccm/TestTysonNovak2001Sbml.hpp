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

#include <algorithm>
#include <fstream>
#include <iostream>
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
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"
#include "SbmlTestHelperFunctions.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "Timer.hpp"
#include "WildTypeCellMutationState.hpp"

#include "TysonNovak2001SbmlOdeSystemAndCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace st = sbmltest;

class TestTysonNovak2001Sbml : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 8u;

    void RunOdeWithSolver(AbstractIvpOdeSolver &rSolver, const std::string solverName)
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

            std::vector<std::vector<double>> solutions;
            std::vector<double> times;

            // Run ODE until it stops, then start again with updated initial conditions
            for (unsigned i = 0; i < 2; i++)
            {
                initial_conditions = ode_system.GetInitialConditions();

                Timer::Reset();
                ode_solution = rSolver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
                Timer::Print(std::to_string(i) + ". Tyson Novak (" + solverName + ")");

                // ODE should have stopped
                TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), true);

                // Collate solutions and times from all runs
                solutions.insert(solutions.end(), ode_solution.rGetSolutions().begin(), ode_solution.rGetSolutions().end());
                times.insert(times.end(), ode_solution.rGetTimes().begin(), ode_solution.rGetTimes().end());

                // Check first run
                if (i == 0)
                {
                    // Check that the ODE stopped at the right time
                    TS_ASSERT_DELTA(ode_solution.rGetTimes().back(), 103.80, 1e-2);

                    // Compare solutions for first run with Tellurium values
                    for (unsigned j = 0; j < ODE_SIZE; j++)
                    {
                        std::vector<double> values;
                        for (unsigned k = 0; k < solutions.size(); k++)
                        {
                            values.push_back(solutions[k][j]);
                        }
                        double min_val = st::min(values);
                        double max_val = st::max(values);
                        double mean_val = st::mean(values);
                        double std_val = st::stdev(values);
                        double q1_val = st::quantile(values, 0.25);
                        double q2_val = st::quantile(values, 0.5);
                        double q3_val = st::quantile(values, 0.75);

                        // Compare with values from Tellurium
                        if (j == 0) // CycBt
                        {
                            TS_ASSERT_DELTA(min_val, 0.0010, 1e-4);
                            TS_ASSERT_DELTA(max_val, 0.6493, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.1962, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.2210, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.0386, 1e-3);
                            TS_ASSERT_DELTA(q2_val, 0.0393, 1e-3);
                            TS_ASSERT_DELTA(q3_val, 0.3806, 1e-2);
                        }
                        else if (j == 1) // Cdc20a
                        {
                            TS_ASSERT_DELTA(min_val, 0.00000107, 1e-5); // Matlab is -5.51016e-07
                            TS_ASSERT_DELTA(max_val, 0.3022, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.0153, 1e-3);
                            TS_ASSERT_DELTA(std_val, 0.0527, 1e-3);
                            TS_ASSERT_DELTA(q1_val, 0.00000137, 1e-7);
                            TS_ASSERT_DELTA(q2_val, 0.00000174, 1e-7);
                            TS_ASSERT_DELTA(q3_val, 0.000571, 1e-5);
                        }
                        else if (j == 2) // Cdh1
                        {
                            TS_ASSERT_DELTA(min_val, 0.001, 1e-3);
                            TS_ASSERT_DELTA(max_val, 0.9981, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.5790, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.4666, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.00685, 1e-4);
                            TS_ASSERT_DELTA(q2_val, 0.9710, 1e-2);
                            TS_ASSERT_DELTA(q3_val, 0.9953, 1e-2);
                        }
                        else if (j == 3) // m
                        {
                            TS_ASSERT_DELTA(min_val, 0.50, 1e-2);
                            TS_ASSERT_DELTA(max_val, 0.8122, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.6444, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.0900, 1e-3);
                            TS_ASSERT_DELTA(q1_val, 0.5653, 1e-2);
                            TS_ASSERT_DELTA(q2_val, 0.6385, 1e-2);
                            TS_ASSERT_DELTA(q3_val, 0.7205, 1e-2);
                        }
                        else if (j == 4) // Cdc20t
                        {
                            TS_ASSERT_DELTA(min_val, 0.001, 1e-3);
                            TS_ASSERT_DELTA(max_val, 1.5116, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.3198, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.4912, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.0463, 1e-3);
                            TS_ASSERT_DELTA(q2_val, 0.0497, 1e-3);
                            TS_ASSERT_DELTA(q3_val, 0.3572, 1e-2);
                        }
                        else if (j == 5) // IEP
                        {
                            TS_ASSERT_DELTA(min_val, 0.000635, 1e-5);
                            TS_ASSERT_DELTA(max_val, 0.5554, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.1149, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.1879, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.00070, 1e-5);
                            TS_ASSERT_DELTA(q2_val, 0.000930, 1e-5);
                            TS_ASSERT_DELTA(q3_val, 0.18219, 1e-2);
                        }
                        else if (j == 6) // CKIt
                        {
                            TS_ASSERT_DELTA(min_val, 0.001, 1e-3);
                            TS_ASSERT_DELTA(max_val, 0.7162, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.2038, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.1938, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.0382, 1e-3);
                            TS_ASSERT_DELTA(q2_val, 0.1033, 1e-2);
                            TS_ASSERT_DELTA(q3_val, 0.3586, 1e-2);
                        }
                        else if (j == 7) // SK
                        {
                            TS_ASSERT_DELTA(min_val, 0.001, 1e-3);
                            TS_ASSERT_DELTA(max_val, 0.4464, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.0809, 1e-3);
                            TS_ASSERT_DELTA(std_val, 0.1006, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.0160, 1e-3);
                            TS_ASSERT_DELTA(q2_val, 0.0417, 1e-3);
                            TS_ASSERT_DELTA(q3_val, 0.0910, 1e-3);
                        }
                    }

                    // Compare derived quantities for first run with Tellurium values
                    ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);
                    for (const std::string name : {"CycB", "Trimer", "Mad"})
                    {
                        std::vector<double> values = ode_solution.GetAnyVariable(name);
                        double min_val = st::min(values);
                        double max_val = st::max(values);
                        double mean_val = st::mean(values);
                        double std_val = st::stdev(values);
                        double q1_val = st::quantile(values, 0.25);
                        double q2_val = st::quantile(values, 0.5);
                        double q3_val = st::quantile(values, 0.75);

                        if (name == "CycB")
                        {
                            TS_ASSERT_DELTA(min_val, 0.000047, 1e-5);
                            TS_ASSERT_DELTA(max_val, 0.6290, 1e-2);
                            TS_ASSERT_DELTA(mean_val, 0.1572, 1e-2);
                            TS_ASSERT_DELTA(std_val, 0.2270, 1e-2);
                            TS_ASSERT_DELTA(q1_val, 0.000119, 1e-5);
                            TS_ASSERT_DELTA(q2_val, 0.000589, 1e-5);
                            TS_ASSERT_DELTA(q3_val, 0.3423, 1e-2);
                        }
                        else if (name == "Trimer")
                        {
                            TS_ASSERT_DELTA(min_val, 0.000382, 1e-5);
                            TS_ASSERT_DELTA(max_val, 0.09026, 1e-3);
                            TS_ASSERT_DELTA(mean_val, 0.0390, 1e-3);
                            TS_ASSERT_DELTA(std_val, 0.0138, 1e-3);
                            TS_ASSERT_DELTA(q1_val, 0.03488, 1e-3);
                            TS_ASSERT_DELTA(q2_val, 0.03850, 1e-3);
                            TS_ASSERT_DELTA(q3_val, 0.03862, 1e-3);
                        }
                        else if (name == "Mad")
                        {
                            TS_ASSERT_DELTA(min_val, 1.0, 1e-3);
                            TS_ASSERT_DELTA(max_val, 1.0, 1e-3);
                            TS_ASSERT_DELTA(mean_val, 1.0, 1e-3);
                            TS_ASSERT_DELTA(std_val, 0.0, 1e-3);
                            TS_ASSERT_DELTA(q1_val, 1.0, 1e-3);
                            TS_ASSERT_DELTA(q2_val, 1.0, 1e-3);
                            TS_ASSERT_DELTA(q3_val, 1.0, 1e-3);
                        }
                    }
                }

                // Update time for next run
                start_time = ode_solution.rGetTimes().back();
                end_time = start_time + run_length;
            }
            // Check that the last run stopped at the right time
            TS_ASSERT_DELTA(ode_solution.rGetTimes().back(), 250.32, 1e-2);

            // The following code provides nice output for gnuplot
            // use the command
            // plot "tysonnovak_cvode.dat" u 1:2 etc. for the various species...
            // or
            // plot "tysonnovak_cvode.dat" u 1:2, "" u 1:3, "" u 1:4 etc. for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("tysonnovak_" + solverName + ".dat");
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
    void TestCellCycleModel()
    {
        // Setup time
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        const unsigned num_timesteps = 10000;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(300.0, num_timesteps);

        // Create a healthy cell
        auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        auto p_cell_0 = boost::make_shared<Cell>(p_wild_state, new TysonNovak2001SbmlCellCycleModel);
        p_cell_0->SetCellProliferativeType(p_stem_type);

        // Set up the cell cycle model - this should use CVODE by default
        auto p_ccm_0 = static_cast<TysonNovak2001SbmlCellCycleModel *>(p_cell_0->GetCellCycleModel());
        p_ccm_0->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_0->CanCellTerminallyDifferentiate(), false);

        p_cell_0->InitialiseCellCycleModel();
        p_ccm_0->SetDt(0.01);

        // Create another cell with a cell-cycle model that uses a BackwardEulerIvpOdeSolver
        auto solver = CellCycleModelOdeSolver<TysonNovak2001SbmlCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovak2001SbmlCellCycleModel, BackwardEulerIvpOdeSolver>> p_solver(solver);
        p_solver->SetSizeOfOdeSystem(ODE_SIZE);
        p_solver->Initialise();

        auto p_cell_1 = boost::make_shared<Cell>(p_wild_state, new TysonNovak2001SbmlCellCycleModel(p_solver));
        p_cell_1->SetCellProliferativeType(p_stem_type);

        auto p_ccm_1 = static_cast<TysonNovak2001SbmlCellCycleModel *>(p_cell_1->GetCellCycleModel());
        p_ccm_1->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_1->CanCellTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(p_ccm_1->GetOdeSolver()->GetSizeOfOdeSystem(), ODE_SIZE);

        p_cell_1->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_ccm_1->GetDt(), 0.0001); // Timestep for non-adaptive solvers defaults to 0.0001
        p_ccm_1->SetDt(0.01);

        // Test the cell is ready to divide at the right time
        double standard_divide_time = 103.80;
        double tolerance = 0.05;
        for (unsigned i = 0; i < num_timesteps / 2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool division_ready_0 = p_ccm_0->ReadyToDivide();
            bool division_ready_1 = p_ccm_1->ReadyToDivide();

            if (time > standard_divide_time + tolerance)
            {
                TS_ASSERT_EQUALS(division_ready_0, true);
                TS_ASSERT_EQUALS(division_ready_1, true);
            }
            else if (time < standard_divide_time - tolerance)
            {
                TS_ASSERT_EQUALS(division_ready_0, false);
                TS_ASSERT_EQUALS(division_ready_1, false);
            }
        }

        // Check CVODE vs BackwardEuler solution
        std::vector<double> proteins_0 = p_ccm_0->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_0.size(), ODE_SIZE);
        TS_ASSERT_DELTA(proteins_0[0], 0.1789, 1e-4); // CycBt
        TS_ASSERT_DELTA(proteins_0[1], 0.3039, 1e-4); // Cdc20a
        TS_ASSERT_DELTA(proteins_0[2], 0.4455, 1e-4); // Cdh1
        TS_ASSERT_DELTA(proteins_0[3], 0.8125, 1e-4); // m
        TS_ASSERT_DELTA(proteins_0[4], 1.1626, 1e-4); // Cdc20t
        TS_ASSERT_DELTA(proteins_0[5], 0.5465, 1e-4); // IEP
        TS_ASSERT_DELTA(proteins_0[6], 0.0800, 1e-4); // CKIt
        TS_ASSERT_DELTA(proteins_0[7], 0.0816, 1e-4); // SK

        std::vector<double> proteins_1 = p_ccm_1->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_1.size(), ODE_SIZE);
        TS_ASSERT_DELTA(proteins_1[0], proteins_0[0], 1e-2); // CycBt
        TS_ASSERT_DELTA(proteins_1[1], proteins_0[1], 1e-2); // Cdc20a
        TS_ASSERT_DELTA(proteins_1[2], proteins_0[2], 1e-2); // Cdh1
        TS_ASSERT_DELTA(proteins_1[3], proteins_0[3], 1e-2); // m
        TS_ASSERT_DELTA(proteins_1[4], proteins_0[4], 1e-2); // Cdc20t
        TS_ASSERT_DELTA(proteins_1[5], proteins_0[5], 1e-2); // IEP
        TS_ASSERT_DELTA(proteins_1[6], proteins_0[6], 1e-2); // CKIt
        TS_ASSERT_DELTA(proteins_1[7], proteins_0[7], 1e-2); // SK

        // Test for a mutant cell
        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), true);
        p_ccm_0->ResetForDivision();
        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), false);

        auto p_mutation = boost::make_shared<ApcOneHitCellMutationState>();
        auto p_ccm_2 = static_cast<TysonNovak2001SbmlCellCycleModel *>(p_ccm_0->CreateCellCycleModel());
        auto p_cell_2 = boost::make_shared<Cell>(p_mutation, p_ccm_2);
        p_cell_2->SetCellProliferativeType(p_stem_type);

        TS_ASSERT_EQUALS(p_cell_2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_ccm_2->ReadyToDivide(), false);

        // Test the cell is ready to divide at the right time
        standard_divide_time = 243.69;
        for (unsigned i = 0; i < num_timesteps / 2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool division_ready_0 = p_ccm_0->ReadyToDivide();
            bool division_ready_2 = p_ccm_2->ReadyToDivide();

            if (time > standard_divide_time)
            {
                TS_ASSERT_EQUALS(division_ready_0, true);
                TS_ASSERT_EQUALS(division_ready_2, true);
                break;
            }
            else if (time < standard_divide_time)
            {
                TS_ASSERT_EQUALS(division_ready_0, false);
                TS_ASSERT_EQUALS(division_ready_2, false);
            }
        }

        // Check ODE solution
        proteins_0 = p_ccm_0->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_0.size(), ODE_SIZE);
        TS_ASSERT_DELTA(proteins_0[0], 0.1153, 1e-4); // CycBt
        TS_ASSERT_DELTA(proteins_0[1], 0.2764, 1e-4); // Cdc20a
        TS_ASSERT_DELTA(proteins_0[2], 0.0397, 1e-4); // Cdh1
        TS_ASSERT_DELTA(proteins_0[3], 1.5108, 1e-4); // m
        TS_ASSERT_DELTA(proteins_0[4], 0.4814, 1e-4); // Cdc20t
        TS_ASSERT_DELTA(proteins_0[5], 0.5205, 1e-4); // IEP
        TS_ASSERT_DELTA(proteins_0[6], 0.0155, 1e-4); // CKIt
        TS_ASSERT_DELTA(proteins_0[7], 0.9831, 1e-4); // SK

        std::vector<double> proteins_2 = p_ccm_2->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_2.size(), ODE_SIZE);
        TS_ASSERT_DELTA(proteins_2[0], proteins_0[0], 1e-4); // CycBt
        TS_ASSERT_DELTA(proteins_2[1], proteins_0[1], 1e-4); // Cdc20a
        TS_ASSERT_DELTA(proteins_2[2], proteins_0[2], 1e-4); // Cdh1
        TS_ASSERT_DELTA(proteins_2[3], proteins_0[3], 1e-4); // m
        TS_ASSERT_DELTA(proteins_2[4], proteins_0[4], 1e-4); // Cdc20t
        TS_ASSERT_DELTA(proteins_2[5], proteins_0[5], 1e-4); // IEP
        TS_ASSERT_DELTA(proteins_2[6], proteins_0[6], 1e-4); // CKIt
        TS_ASSERT_DELTA(proteins_2[7], proteins_0[7], 1e-4); // SK
    }

    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak_2001_ode.arch";

        {
            std::vector<double> state_variables = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

            TysonNovak2001SbmlOdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);
            TS_ASSERT_DELTA(initial_conditions[0], 0.001, 1e-3); // CycBt
            TS_ASSERT_DELTA(initial_conditions[1], 0.001, 1e-3); // Cdc20a
            TS_ASSERT_DELTA(initial_conditions[2], 3.25, 1e-3);  // Cdh1
            TS_ASSERT_DELTA(initial_conditions[3], 0.5, 1e-3);   // m
            TS_ASSERT_DELTA(initial_conditions[4], 0.001, 1e-3); // Cdc20t
            TS_ASSERT_DELTA(initial_conditions[5], 0.001, 1e-3); // IEP
            TS_ASSERT_DELTA(initial_conditions[6], 0.001, 1e-3); // CKIt
            TS_ASSERT_DELTA(initial_conditions[7], 0.001, 1e-3); // SK

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
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);
            TS_ASSERT_DELTA(initial_conditions[0], 0.001, 1e-3); // CycBt
            TS_ASSERT_DELTA(initial_conditions[1], 0.001, 1e-3); // Cdc20a
            TS_ASSERT_DELTA(initial_conditions[2], 0.001, 1e-3); // Cdh1
            TS_ASSERT_DELTA(initial_conditions[3], 0.5, 1e-3);   // m
            TS_ASSERT_DELTA(initial_conditions[4], 0.001, 1e-3); // Cdc20t
            TS_ASSERT_DELTA(initial_conditions[5], 0.001, 1e-3); // IEP
            TS_ASSERT_DELTA(initial_conditions[6], 0.001, 1e-3); // CKIt
            TS_ASSERT_DELTA(initial_conditions[7], 0.001, 1e-3); // SK

            double var0 = p_ode_system->GetStateVariable(0); // CycBt
            double var1 = p_ode_system->GetStateVariable(1); // Cdc20a
            double var2 = p_ode_system->GetStateVariable(2); // Cdh1
            double var3 = p_ode_system->GetStateVariable(3); // m
            double var4 = p_ode_system->GetStateVariable(4); // Cdc20t
            double var5 = p_ode_system->GetStateVariable(5); // IEP
            double var6 = p_ode_system->GetStateVariable(6); // CKIt
            double var7 = p_ode_system->GetStateVariable(7); // SK

            TS_ASSERT_DELTA(var0, 0.0, 1e-3); // CycBt
            TS_ASSERT_DELTA(var1, 1.0, 1e-3); // Cdc20a
            TS_ASSERT_DELTA(var2, 2.0, 1e-3); // Cdh1
            TS_ASSERT_DELTA(var3, 3.0, 1e-3); // m
            TS_ASSERT_DELTA(var4, 4.0, 1e-3); // Cdc20t
            TS_ASSERT_DELTA(var5, 5.0, 1e-3); // IEP
            TS_ASSERT_DELTA(var6, 6.0, 1e-3); // CKIt
            TS_ASSERT_DELTA(var7, 7.0, 1e-3); // SK

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        TysonNovak2001SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.001); // CycBt
        initial_conditions.push_back(0.001); // Cdc20a
        initial_conditions.push_back(0.001); // Cdh1
        initial_conditions.push_back(0.5);   // m
        initial_conditions.push_back(0.001); // Cdc20t
        initial_conditions.push_back(0.001); // IEP
        initial_conditions.push_back(0.001); // CKIt
        initial_conditions.push_back(0.001); // SK

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Compare derivatives with values from Tellurium
        TS_ASSERT_DELTA(derivs[0], 3.99580000e-02, 1e-3);  // CycBt
        TS_ASSERT_DELTA(derivs[1], -2.50100000e-01, 1e-2); // Cdc20a
        TS_ASSERT_DELTA(derivs[2], 9.70803883e-01, 1e-2);  // Cdh1
        TS_ASSERT_DELTA(derivs[3], 2.37500000e-03, 1e-4);  // m
        TS_ASSERT_DELTA(derivs[4], 4.90000000e-03, 1e-4);  // Cdc20t
        TS_ASSERT_DELTA(derivs[5], 1.08707977e-05, 1e-6);  // IEP
        TS_ASSERT_DELTA(derivs[6], 9.99719098e-01, 1e-2);  // CKIt
        TS_ASSERT_DELTA(derivs[7], 2.77174939e-02, 1e-3);  // SK

        // Check derived quantities
        TS_ASSERT_EQUALS(ode_system.GetNumberOfDerivedQuantities(), 3u);
        TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("CycB"), 0u);
        TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("Trimer"), 1u);
        TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("Mad"), 2u);

        std::vector<double> derived_quantities = ode_system.ComputeDerivedQuantities(0.0, initial_conditions);
        TS_ASSERT_EQUALS(derived_quantities.size(), 3u);
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

#endif // TESTTYSONNOVAK2001SBML_HPP_
