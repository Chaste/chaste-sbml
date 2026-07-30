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

#ifndef TEST_TYSON_NOVAK_2001_SBML_CELL_CYCLE_MODEL_HPP_
#define TEST_TYSON_NOVAK_2001_SBML_CELL_CYCLE_MODEL_HPP_

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "OutputFileHandler.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "TysonNovak2001SbmlCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace
{
// Native time units per hour: the model's time is in minutes, but Chaste integrates in hours.
constexpr double TIMESCALE_MULTIPLIER = 60.0;
} // namespace

class TestTysonNovak2001SbmlCellCycleModel : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 8u;

public:
    void TestCellCycleModel()
    {
        // Setup time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        // The model is in minutes but Chaste integrates in hours (the ODE system scales derivatives
        // by 60), so divide the end time and timestep by 60 to integrate over the same native span.
        // The simulation time (hours) is converted back to minutes below to compare against the
        // native divide times.
        const unsigned num_timesteps = 10000;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(300.0 / TIMESCALE_MULTIPLIER, num_timesteps);

        // Create a healthy cell
        auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        auto p_cell_0 = boost::make_shared<Cell>(p_wild_state, new TysonNovak2001SbmlCellCycleModel);
        p_cell_0->SetCellProliferativeType(p_stem_type);

        // Set up the cell cycle model - this should use CVODE by default
        auto p_ccm_0 = static_cast<TysonNovak2001SbmlCellCycleModel*>(p_cell_0->GetCellCycleModel());
        p_ccm_0->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_0->CanCellTerminallyDifferentiate(), false);

        p_cell_0->InitialiseCellCycleModel();
        p_ccm_0->SetDt(0.01 / TIMESCALE_MULTIPLIER);

        // Create another cell with a cell-cycle model that uses a BackwardEulerIvpOdeSolver
        auto solver = CellCycleModelOdeSolver<TysonNovak2001SbmlCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovak2001SbmlCellCycleModel, BackwardEulerIvpOdeSolver> > p_solver(solver);
        p_solver->SetSizeOfOdeSystem(ODE_SIZE);
        p_solver->Initialise();

        auto p_cell_1 = boost::make_shared<Cell>(p_wild_state, new TysonNovak2001SbmlCellCycleModel(p_solver));
        p_cell_1->SetCellProliferativeType(p_stem_type);

        auto p_ccm_1 = static_cast<TysonNovak2001SbmlCellCycleModel*>(p_cell_1->GetCellCycleModel());
        p_ccm_1->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm_1->CanCellTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(p_ccm_1->GetOdeSolver()->GetSizeOfOdeSystem(), ODE_SIZE);

        p_cell_1->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_ccm_1->GetDt(), 0.0001); // Timestep for non-adaptive solvers defaults to 0.0001
        p_ccm_1->SetDt(0.01 / TIMESCALE_MULTIPLIER);

        // Test the cell is ready to divide at the right time
        double standard_divide_time = 103.80;
        double tolerance = 0.05;
        for (unsigned i = 0; i < num_timesteps / 2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime() * TIMESCALE_MULTIPLIER; // hours -> native minutes

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
        TS_ASSERT_DELTA(proteins_0[0], 0.1791, 1e-4); // CycBt
        TS_ASSERT_DELTA(proteins_0[1], 0.3039, 1e-4); // Cdc20a
        TS_ASSERT_DELTA(proteins_0[2], 0.4435, 1e-4); // Cdh1
        TS_ASSERT_DELTA(proteins_0[3], 0.8125, 1e-4); // m
        TS_ASSERT_DELTA(proteins_0[4], 1.1629, 1e-4); // Cdc20t
        TS_ASSERT_DELTA(proteins_0[5], 0.5465, 1e-4); // IEP
        TS_ASSERT_DELTA(proteins_0[6], 0.0800, 1e-4); // CKIt
        TS_ASSERT_DELTA(proteins_0[7], 0.0811, 1e-4); // SK

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
        auto p_ccm_2 = static_cast<TysonNovak2001SbmlCellCycleModel*>(p_ccm_0->CreateCellCycleModel());
        auto p_cell_2 = boost::make_shared<Cell>(p_mutation, p_ccm_2);
        p_cell_2->SetCellProliferativeType(p_stem_type);

        TS_ASSERT_EQUALS(p_cell_2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_ccm_2->ReadyToDivide(), false);

        // Test the cell is ready to divide at the right time
        standard_divide_time = 250.32;
        for (unsigned i = 0; i < num_timesteps / 2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime() * TIMESCALE_MULTIPLIER; // hours -> native minutes

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
        TS_ASSERT_DELTA(proteins_0[0], 0.1797, 1e-4); // CycBt
        TS_ASSERT_DELTA(proteins_0[1], 0.2977, 1e-4); // Cdc20a
        TS_ASSERT_DELTA(proteins_0[2], 0.4312, 1e-4); // Cdh1
        TS_ASSERT_DELTA(proteins_0[3], 0.8096, 1e-4); // m
        TS_ASSERT_DELTA(proteins_0[4], 1.1312, 1e-4); // Cdc20t
        TS_ASSERT_DELTA(proteins_0[5], 0.5448, 1e-4); // IEP
        TS_ASSERT_DELTA(proteins_0[6], 0.0806, 1e-4); // CKIt
        TS_ASSERT_DELTA(proteins_0[7], 0.0791, 1e-4); // SK

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

    void TestBaseClassMethods()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        auto p_cell = boost::make_shared<Cell>(p_wild_state, new TysonNovak2001SbmlCellCycleModel);
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

        auto p_ccm = static_cast<TysonNovak2001SbmlCellCycleModel*>(p_cell->GetCellCycleModel());

        // Base-class accessors (defaults defined on AbstractSbmlCellCycleModel)
        TS_ASSERT_DELTA(p_ccm->GetAverageTransitCellCycleTime(), 1.25, 1e-9);
        TS_ASSERT_DELTA(p_ccm->GetAverageStemCellCycleTime(), 1.25, 1e-9);
        TS_ASSERT_THROWS_NOTHING(p_ccm->GetStateVariable("CycBt"));

        // Base-class parameter output
        OutputFileHandler handler("TestTysonNovakCcOutputParameters", false);
        out_stream parameter_file = handler.OpenOutputFile("tyson_novak_cc_results.parameters");
        TS_ASSERT_THROWS_NOTHING(p_ccm->OutputCellCycleModelParameters(parameter_file));
        parameter_file->close();
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak_cc.arch";

        {
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            AbstractCellCycleModel* const p_model = new TysonNovak2001SbmlCellCycleModel;
            p_model->SetDimension(3);
            p_model->SetBirthTime(-1.0);

            auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
            auto p_stem_type = boost::make_shared<StemCellProliferativeType>();
            CellPtr p_cell(new Cell(p_wild_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->InitialiseCellCycleModel();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_model;

            SimulationTime::Destroy();
        }

        {
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_model2;

            TS_ASSERT_EQUALS(p_model2->GetDimension(), 3u);
            delete p_model2;
        }
    }
};

#endif // TEST_TYSON_NOVAK_2001_SBML_CELL_CYCLE_MODEL_HPP_
