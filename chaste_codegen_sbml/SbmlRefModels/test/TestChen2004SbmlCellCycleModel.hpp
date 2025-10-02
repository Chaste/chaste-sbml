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

#ifndef TESTCHEN2004SBMLCELLCYCLEMODEL_HPP_
#define TESTCHEN2004SBMLCELLCYCLEMODEL_HPP_

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Chen2004SbmlCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestChen2004SbmlCellCycleModel : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 36u;

public:
    void TestCellCycleModel()
    {
        // Setup time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        const double dt = 0.01;
        const double end_time = 220.0;
        const unsigned num_timesteps = static_cast<unsigned>(end_time / dt);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create a healthy cell
        auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        auto p_cell_0 = boost::make_shared<Cell>(p_wild_state, new Chen2004SbmlCellCycleModel);
        p_cell_0->SetCellProliferativeType(p_stem_type);

        // Set up the cell cycle model - this should use CVODE by default
        auto p_ccm_0 = static_cast<Chen2004SbmlCellCycleModel*>(p_cell_0->GetCellCycleModel());
        p_ccm_0->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_0->CanCellTerminallyDifferentiate(), false);

        p_cell_0->InitialiseCellCycleModel();
        p_ccm_0->SetDt(dt);

        // Create another cell with a cell-cycle model that uses a BackwardEulerIvpOdeSolver
        auto solver = CellCycleModelOdeSolver<Chen2004SbmlCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        boost::shared_ptr<CellCycleModelOdeSolver<Chen2004SbmlCellCycleModel, BackwardEulerIvpOdeSolver> > p_solver(solver);
        p_solver->SetSizeOfOdeSystem(ODE_SIZE);
        p_solver->Initialise();

        auto p_cell_1 = boost::make_shared<Cell>(p_wild_state, new Chen2004SbmlCellCycleModel(p_solver));
        p_cell_1->SetCellProliferativeType(p_stem_type);

        auto p_ccm_1 = static_cast<Chen2004SbmlCellCycleModel*>(p_cell_1->GetCellCycleModel());
        p_ccm_1->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_1->CanCellTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(p_ccm_1->GetOdeSolver()->GetSizeOfOdeSystem(), ODE_SIZE);

        p_cell_1->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_ccm_1->GetDt(), 0.0001); // Timestep for non-adaptive solvers defaults to 0.0001
        p_ccm_1->SetDt(dt);

        // Test the cell is ready to divide at the right time
        double standard_divide_time = 101.20;
        double divide_time_tolerance = 0.02;
        for (unsigned i = 0; i < num_timesteps / 2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool division_ready_0 = p_ccm_0->ReadyToDivide();
            bool division_ready_1 = p_ccm_1->ReadyToDivide();

            if (time > standard_divide_time + divide_time_tolerance)
            {
                TS_ASSERT_EQUALS(division_ready_0, true);
                TS_ASSERT_EQUALS(division_ready_1, true);
            }
            else if (time < standard_divide_time - divide_time_tolerance)
            {
                TS_ASSERT_EQUALS(division_ready_0, false);
                TS_ASSERT_EQUALS(division_ready_1, false);
            }
        }

        // Check CVODE vs BackwardEuler solution
        std::vector<double> proteins_0 = p_ccm_0->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_0.size(), ODE_SIZE);

        std::vector<double> proteins_1 = p_ccm_1->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_1.size(), ODE_SIZE);

        std::vector<std::vector<double> > expected_solution = {
            { 1.525341, 1e-2 },  // BUD
            { 0.210599, 1e-2 },  // C2
            { 0.037114, 1e-3 },  // C2P
            { 0.050391, 1e-3 },  // C5
            { 0.008844, 1e-4 },  // C5P
            { 0.514681, 1e-2 },  // CDC14
            { 0.661787, 1e-2 },  // CDC15
            { 0.462518, 1e-2 },  // CDC20
            { 1.506304, 1e-2 },  // CDC20i
            { 0.061446, 1e-3 },  // CDC6
            { 0.017502, 1e-3 },  // CDC6P
            { 0.798235, 1e-2 },  // CDH1
            { 0.201764, 1e-2 },  // CDH1i
            { 0.302153, 1e-2 },  // CLB2
            { 0.075761, 1e-3 },  // CLB5
            { 0.064637, 1e-3 },  // CLN2
            { 0.277167, 1e-2 },  // ESP1
            { 0.218717, 1e-2 },  // F2
            { 0.043554, 1e-3 },  // F2P
            { 0.000052, 1e-6 },  // F5
            { 0.000010, 1e-6 },  // F5P
            { 0.136154, 1e-2 },  // IEP
            { 2.629131, 1e-2 },  // MASS
            { 0.015425, 1e-3 },  // NET1
            { 1.291491, 1e-2 },  // NET1P
            { 22.272502, 1e-1 }, // ORI
            { 0.029937, 1e-3 },  // PDS1
            { 0.128275, 1e-2 },  // PPX
            { 0.870067, 1e-2 },  // RENT
            { 0.615256, 1e-2 },  // RENTP
            { 0.011849, 1e-3 },  // SIC1
            { 0.007308, 1e-4 },  // SIC1P
            { 1.298397, 1e-2 },  // SPN
            { 0.952391, 1e-2 },  // SWI5
            { 0.027953, 1e-3 },  // SWI5P
            { 0.977834, 1e-2 },  // TEM1GTP
        };

        std::vector<std::string> var_names = p_ccm_0->GetOdeSystem()->rGetStateVariableNames();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            const char* var_name = var_names[i].c_str();
            double exp_val = expected_solution[i][0];
            double exp_tol = expected_solution[i][1];

            TSM_ASSERT_DELTA(var_name, proteins_0[i], exp_val, exp_tol);
            TSM_ASSERT_DELTA(var_name, proteins_1[i], exp_val, exp_tol);
        }

        // Test for a mutant cell
        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), true);
        p_ccm_0->ResetForDivision();
        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), false);

        auto p_mutation = boost::make_shared<ApcOneHitCellMutationState>();
        auto p_ccm_2 = static_cast<Chen2004SbmlCellCycleModel*>(p_ccm_0->CreateCellCycleModel());
        auto p_cell_2 = boost::make_shared<Cell>(p_mutation, p_ccm_2);
        p_cell_2->SetCellProliferativeType(p_stem_type);

        TS_ASSERT_EQUALS(p_cell_0->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_cell_2->ReadyToDivide(), false);

        TS_ASSERT_EQUALS(p_ccm_0->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_ccm_2->ReadyToDivide(), false);

        std::vector<double> state_0 = p_ccm_0->GetOdeSystem()->GetStateVariables();
        std::vector<double> state_2 = p_ccm_2->GetOdeSystem()->GetStateVariables();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            TS_ASSERT_DELTA(state_0[i], state_2[i], 1e-6);
        }

        // Test the cell is ready to divide at the right time
        standard_divide_time = 202.14;
        divide_time_tolerance = 0.02;
        for (unsigned i = 0; i < num_timesteps / 2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool division_ready_0 = p_ccm_0->ReadyToDivide();
            bool division_ready_2 = p_ccm_2->ReadyToDivide();

            if (time > standard_divide_time + divide_time_tolerance)
            {
                TS_ASSERT_EQUALS(division_ready_0, true);
                TS_ASSERT_EQUALS(division_ready_2, true);
            }
            else if (time < standard_divide_time - divide_time_tolerance)
            {
                TS_ASSERT_EQUALS(division_ready_0, false);
                TS_ASSERT_EQUALS(division_ready_2, false);
            }
        }

        // Check ODE solution
        proteins_0 = p_ccm_0->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_0.size(), ODE_SIZE);

        std::vector<double> proteins_2 = p_ccm_2->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins_2.size(), ODE_SIZE);

        expected_solution = {
            { 1.510533, 1e-2 },  // BUD
            { 0.210634, 1e-2 },  // C2
            { 0.037258, 1e-3 },  // C2P
            { 0.049837, 1e-3 },  // C5
            { 0.008779, 1e-4 },  // C5P
            { 0.511999, 1e-2 },  // CDC14
            { 0.661786, 1e-2 },  // CDC15
            { 0.461485, 1e-2 },  // CDC20
            { 1.507291, 1e-2 },  // CDC20i
            { 0.061318, 1e-3 },  // CDC6
            { 0.017538, 1e-3 },  // CDC6P
            { 0.795886, 1e-2 },  // CDH1
            { 0.204114, 1e-2 },  // CDH1i
            { 0.302283, 1e-2 },  // CLB2
            { 0.074949, 1e-3 },  // CLB5
            { 0.063794, 1e-3 },  // CLN2
            { 0.277556, 1e-2 },  // ESP1
            { 0.218300, 1e-2 },  // F2
            { 0.043643, 1e-3 },  // F2P
            { 0.000051, 1e-6 },  // F5
            { 0.000010, 1e-6 },  // F5P
            { 0.135086, 1e-2 },  // IEP
            { 2.623296, 1e-2 },  // MASS
            { 0.015549, 1e-3 },  // NET1
            { 1.296074, 1e-2 },  // NET1P
            { 22.172016, 1e-2 }, // ORI
            { 0.029879, 1e-3 },  // PDS1
            { 0.128173, 1e-2 },  // PPX
            { 0.872748, 1e-2 },  // RENT
            { 0.615252, 1e-2 },  // RENTP
            { 0.011854, 1e-3 },  // SIC1
            { 0.007335, 1e-4 },  // SIC1P
            { 1.299038, 1e-2 },  // SPN
            { 0.952581, 1e-2 },  // SWI5
            { 0.028176, 1e-3 },  // SWI5P
            { 0.977834, 1e-2 },  // TEM1GTP
        };

        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            const char* var_name = var_names[i].c_str();
            double exp_val = expected_solution[i][0];
            double exp_tol = expected_solution[i][1];

            TSM_ASSERT_DELTA(var_name, proteins_0[i], exp_val, exp_tol);
            TSM_ASSERT_DELTA(var_name, proteins_2[i], exp_val, exp_tol);
        }
    }
};

#endif // TESTCHEN2004SBMLCELLCYCLEMODEL_HPP_
