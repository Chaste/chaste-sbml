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

#ifndef TEST_GARDNER_1998_SBML_CELL_CYCLE_MODEL_HPP_
#define TEST_GARDNER_1998_SBML_CELL_CYCLE_MODEL_HPP_

#include <fstream>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Gardner1998SbmlCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace
{
// Native time units per hour: the model's time is in seconds, but Chaste integrates in hours.
constexpr double TIMESCALE_MULTIPLIER = 3600.0;
} // namespace

class TestGardner1998SbmlCellCycleModel : public AbstractCellBasedTestSuite
{
public:
    void TestCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(Gardner1998SbmlCellCycleModel cell_cycle_model);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        // The model is in seconds but Chaste integrates in hours (the ODE system scales derivatives
        // by 3600), so divide the timestep and end time by 3600 to integrate over the same span.
        const double dt = 0.01 / TIMESCALE_MULTIPLIER;
        const double end_time = 20.0 / TIMESCALE_MULTIPLIER;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, static_cast<unsigned>(end_time / dt));

        auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        auto p_cell = boost::make_shared<Cell>(p_wild_state, new Gardner1998SbmlCellCycleModel);
        p_cell->SetCellProliferativeType(p_stem_type);

        auto p_ccm = static_cast<Gardner1998SbmlCellCycleModel*>(p_cell->GetCellCycleModel());
        p_ccm->SetBirthTime(p_simulation_time->GetTime());
        TS_ASSERT_EQUALS(p_ccm->CanCellTerminallyDifferentiate(), false);

        p_cell->InitialiseCellCycleModel();
        p_ccm->SetDt(dt);

        // Base-class accessors
        TS_ASSERT_DELTA(p_ccm->GetAverageTransitCellCycleTime(), 1.25, 1e-9);
        TS_ASSERT_DELTA(p_ccm->GetAverageStemCellCycleTime(), 1.25, 1e-9);
        TS_ASSERT_THROWS_NOTHING(p_ccm->GetStateVariable("C"));

        // Advance a few steps so the ODE system is integrated
        for (unsigned i = 0; i < 10; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_THROWS_NOTHING(p_ccm->ReadyToDivide());
        }

        // Copy via CreateCellCycleModel (exercises the protected copy-constructor)
        Gardner1998SbmlCellCycleModel* p_ccm2 = static_cast<Gardner1998SbmlCellCycleModel*>(p_ccm->CreateCellCycleModel());
        TS_ASSERT(p_ccm2 != nullptr);
        delete p_ccm2;

        // Parameter output
        OutputFileHandler handler("TestGardner1998CcOutputParameters", false);
        out_stream parameter_file = handler.OpenOutputFile("gardner_1998_cc_results.parameters");
        TS_ASSERT_THROWS_NOTHING(p_ccm->OutputCellCycleModelParameters(parameter_file));
        parameter_file->close();
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "gardner_1998_cc.arch";

        {
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            AbstractCellCycleModel* const p_model = new Gardner1998SbmlCellCycleModel;
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
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.0, 1e-12);

            delete p_model2;
        }
    }
};

#endif // TEST_GARDNER_1998_SBML_CELL_CYCLE_MODEL_HPP_
