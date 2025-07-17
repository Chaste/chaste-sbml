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

#ifndef TESTCHEN2004SBML_HPP_
#define TESTCHEN2004SBML_HPP_

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

#include "Chen2004SbmlOdeSystemAndCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace st = sbmltest;

class TestChen2004Sbml : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 36u;

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

    void RunOdeWithSolver(AbstractIvpOdeSolver &rSolver, const std::string solverName)
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
            std::vector<std::vector<double>> solutions;
            std::vector<std::vector<double>> derived_quantities(3);

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
                initial_conditions = ode_system.GetStateVariables();
                start_time = ode_solution.rGetTimes().back();
                end_time = start_time + run_length;
            }

            // Compare end solutions with Tellurium values
            std::vector<double> end_solution = solutions.back();
            std::vector<std::vector<double>> expected_end_solutions = {
                {0.008181, 1e-3},  // BUD
                {0.237895, 1e-2},  // C2
                {0.024148, 1e-3},  // C2P
                {0.069449, 1e-3},  // C5
                {0.006869, 1e-4},  // C5P
                {0.471489, 1e-2},  // CDC14
                {0.656731, 1e-2},  // CDC15
                {0.446186, 1e-2},  // CDC20
                {1.473426, 1e-2},  // CDC20i
                {0.106310, 1e-2},  // CDC6
                {0.015478, 1e-3},  // CDC6P
                {0.931438, 1e-2},  // CDH1
                {0.068561, 1e-3},  // CDH1i ***
                {0.149486, 1e-2},  // CLB2
                {0.052216, 1e-3},  // CLB5
                {0.060328, 1e-3},  // CLN2
                {0.301714, 1e-2},  // ESP1
                {0.235779, 1e-2},  // F2
                {0.027547, 1e-3},  // F2P
                {0.000072, 1e-6},  // F5
                {0.000008, 1e-6},  // F5P
                {0.103716, 1e-2},  // IEP
                {1.207814, 1e-2},  // MASS
                {0.018455, 1e-3},  // NET1
                {1.252664, 1e-2},  // NET1P
                {21.662612, 1e-1}, // ORI
                {0.025638, 1e-3},  // PDS1
                {0.122829, 1e-2},  // PPX
                {1.043842, 1e-2},  // RENT
                {0.484669, 1e-2},  // RENTP
                {0.022570, 1e-3},  // SIC1
                {0.006414, 1e-4},  // SIC1P
                {0.030010, 1e-3},  // SPN
                {0.955855, 1e-2},  // SWI5
                {0.020352, 1e-3},  // SWI5P
                {0.905342, 1e-2}   // TEM1GTP
            };
            std::vector<std::string> var_names = ode_system.rGetStateVariableNames();

            for (unsigned i = 0; i < end_solution.size(); i++)
            {
                double expected_value = expected_end_solutions[i][0];
                double tolerance = expected_end_solutions[i][1];
                TSM_ASSERT_DELTA(var_names[i].c_str(), end_solution[i], expected_value, tolerance);
            }

            // The following code provides nice output for gnuplot
            // use the command
            // plot "chen_2004_cvode.dat" u 1:2 etc. for the various species...
            // or
            // plot "chen_2004_cvode.dat" u 1:2, "" u 1:3, "" u 1:4 etc. for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("chen_2004_" + solverName + ".dat");
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
    void DontTestCellCycleModel()
    {
        // Setup time
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        const unsigned num_timesteps = 10000;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(300.0, num_timesteps);

        // Create a healthy cell
        auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        auto p_cell_0 = boost::make_shared<Cell>(p_wild_state, new Chen2004SbmlCellCycleModel);
        p_cell_0->SetCellProliferativeType(p_stem_type);

        // Set up the cell cycle model - this should use CVODE by default
        auto p_ccm_0 = static_cast<Chen2004SbmlCellCycleModel *>(p_cell_0->GetCellCycleModel());
        p_ccm_0->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_0->CanCellTerminallyDifferentiate(), false);

        p_cell_0->InitialiseCellCycleModel();
        p_ccm_0->SetDt(0.01);

        // Create another cell with a cell-cycle model that uses a BackwardEulerIvpOdeSolver
        auto solver = CellCycleModelOdeSolver<Chen2004SbmlCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        boost::shared_ptr<CellCycleModelOdeSolver<Chen2004SbmlCellCycleModel, BackwardEulerIvpOdeSolver>> p_solver(solver);
        p_solver->SetSizeOfOdeSystem(ODE_SIZE);
        p_solver->Initialise();

        auto p_cell_1 = boost::make_shared<Cell>(p_wild_state, new Chen2004SbmlCellCycleModel(p_solver));
        p_cell_1->SetCellProliferativeType(p_stem_type);

        auto p_ccm_1 = static_cast<Chen2004SbmlCellCycleModel *>(p_cell_1->GetCellCycleModel());
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
        TS_ASSERT_DELTA(proteins_0[0], 0.1789, 1e-4); // BUD
        TS_ASSERT_DELTA(proteins_0[1], 0.3039, 1e-4); // C2
        TS_ASSERT_DELTA(proteins_0[2], 0.4455, 1e-4); // C2P
        TS_ASSERT_DELTA(proteins_0[3], 0.8125, 1e-4); // C5
        TS_ASSERT_DELTA(proteins_0[4], 1.1626, 1e-4); // C5P
        TS_ASSERT_DELTA(proteins_0[5], 0.5465, 1e-4); // CDC14
        TS_ASSERT_DELTA(proteins_0[6], 0.0800, 1e-4); // CDC15
        TS_ASSERT_DELTA(proteins_0[7], 0.0816, 1e-4); // CDC15i

        std::vector<double> proteins_1 = p_ccm_1->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_1.size(), ODE_SIZE);
        TS_ASSERT_DELTA(proteins_1[0], proteins_0[0], 1e-2); // BUD
        TS_ASSERT_DELTA(proteins_1[1], proteins_0[1], 1e-2); // C2
        TS_ASSERT_DELTA(proteins_1[2], proteins_0[2], 1e-2); // C2P
        TS_ASSERT_DELTA(proteins_1[3], proteins_0[3], 1e-2); // C5
        TS_ASSERT_DELTA(proteins_1[4], proteins_0[4], 1e-2); // C5P
        TS_ASSERT_DELTA(proteins_1[5], proteins_0[5], 1e-2); // CDC14
        TS_ASSERT_DELTA(proteins_1[6], proteins_0[6], 1e-2); // CDC15
        TS_ASSERT_DELTA(proteins_1[7], proteins_0[7], 1e-2); // CDC15i

        // Test for a mutant cell
        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), true);
        p_ccm_0->ResetForDivision();
        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), false);

        auto p_mutation = boost::make_shared<ApcOneHitCellMutationState>();
        auto p_ccm_2 = static_cast<Chen2004SbmlCellCycleModel *>(p_ccm_0->CreateCellCycleModel());
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
        TS_ASSERT_DELTA(proteins_0[0], 0.1153, 1e-4); // BUD
        TS_ASSERT_DELTA(proteins_0[1], 0.2764, 1e-4); // C2
        TS_ASSERT_DELTA(proteins_0[2], 0.0397, 1e-4); // C2P
        TS_ASSERT_DELTA(proteins_0[3], 1.5108, 1e-4); // C5
        TS_ASSERT_DELTA(proteins_0[4], 0.4814, 1e-4); // C5P
        TS_ASSERT_DELTA(proteins_0[5], 0.5205, 1e-4); // CDC14
        TS_ASSERT_DELTA(proteins_0[6], 0.0155, 1e-4); // CDC15
        TS_ASSERT_DELTA(proteins_0[7], 0.9831, 1e-4); // CDC15i

        std::vector<double> proteins_2 = p_ccm_2->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_2.size(), ODE_SIZE);
        TS_ASSERT_DELTA(proteins_2[0], proteins_0[0], 1e-4); // BUD
        TS_ASSERT_DELTA(proteins_2[1], proteins_0[1], 1e-4); // C2
        TS_ASSERT_DELTA(proteins_2[2], proteins_0[2], 1e-4); // C2P
        TS_ASSERT_DELTA(proteins_2[3], proteins_0[3], 1e-4); // C5
        TS_ASSERT_DELTA(proteins_2[4], proteins_0[4], 1e-4); // C5P
        TS_ASSERT_DELTA(proteins_2[5], proteins_0[5], 1e-4); // CDC14
        TS_ASSERT_DELTA(proteins_2[6], proteins_0[6], 1e-4); // CDC15
        TS_ASSERT_DELTA(proteins_2[7], proteins_0[7], 1e-4); // CDC15i
    }

    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chen_2004_ode.arch";
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
            Chen2004SbmlOdeSystem ode_system(state_variables);
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
        Chen2004SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> derivatives(ODE_SIZE);
        ode_system.EvaluateYDerivatives(time, default_initial_conditions, derivatives);

        // Compare derivatives with values from Tellurium
        std::vector<double> derivatives_expected = {
            0.014292285913187989,    // BUD ***
            0.03807104123616535,     // C2
            -0.02191995883216533,    // C2P
            0.0452736273488603,      // C5
            -0.003342304055100295,   // C5P
            0.251316422776,          // CDC14
            -0.018916645891352057,   // CDC15
            -0.034247066800000024,   // CDC20
            -0.2532503622733521,     // CDC20i
            0.12356839015323366,     // CDC6
            -0.0037330182800774767,  // CDC6P
            0.2193176241948343,      // CDH1
            -0.2193176141948343,     // CDH1i
            -0.2648597095765076,     // CLB2
            -0.051353622636901804,   // CLB5
            -0.006941222654334056,   // CLN2
            0.058558812262920035,    // ESP1
            0.0032043005026250737,   // F2
            -0.02794253663182511,    // F2P
            0.0001780335569909064,   // F5
            -0.00021441434591090637, // F5P
            -0.062337476806313685,   // IEP
            0.009288318550574718,    // MASS
            -0.11053817267683463,    // NET1
            0.37848931545283476,     // NET1P
            0.22541840999999999,     // ORI
            -0.014060356081539327,   // PDS1
            -0.010198232821435102,   // PPX
            0.300929157138497,       // RENT
            -0.5640339799144971,     // RENTP
            0.032294390197855995,    // SIC1
            -0.0016987958956159834,  // SIC1P
            0.04940637021748367,     // SPN
            -0.021724458793113615,   // SWI5
            -0.01335493175,          // SWI5P
            -0.13000000000000003     // TEM1GTP
        };

        std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            TSM_ASSERT_DELTA(var_names[i].c_str(), derivatives[i], derivatives_expected[i], 1e-6);
        }

        // Check derived quantities
        // TS_ASSERT_EQUALS(ode_system.GetNumberOfDerivedQuantities(), 3u);
        // TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("CycB"), 0u);
        // TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("Trimer"), 1u);
        // TS_ASSERT_EQUALS(ode_system.GetDerivedQuantityIndex("Mad"), 2u);

        // std::vector<double> derived_quantities = ode_system.ComputeDerivedQuantities(0.0, initial_conditions);
        // TS_ASSERT_EQUALS(derived_quantities.size(), 3u);
        // TS_ASSERT_DELTA(derived_quantities[0], 0.001, 1e-3); // CycB
        // TS_ASSERT_DELTA(derived_quantities[1], 0.001, 1e-3); // Trimer
        // TS_ASSERT_DELTA(derived_quantities[2], 1.0, 1e-3);   // Mad
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

#endif // TESTCHEN2004SBML_HPP_
