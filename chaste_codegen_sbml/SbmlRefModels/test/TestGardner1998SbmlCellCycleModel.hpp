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

#ifndef TESTGARDNER1998SBMLCELLCYCLEMODEL_HPP_
#define TESTGARDNER1998SBMLCELLCYCLEMODEL_HPP_

#include <iostream>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Gardner1998SbmlCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestGardner1998SbmlCellCycleModel : public AbstractCellBasedTestSuite
{
    // TODO: Add tests
public:
    void TestCellCycleModelCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(Gardner1998SbmlCellCycleModel cell_cycle_model);
    }

    void TestCellCycleArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "gardner_1998_srn.arch";

        double C0, X0, M0, Y0, Z0;

        // Save archive
        {
            // Setup time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            const double end_time = 10.0;
            const unsigned num_timesteps = 100;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            // Create a healthy cell so we can initialise the cell-cycle model's ODE system
            auto p_wild_state = boost::make_shared<WildTypeCellMutationState>();
            auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

            auto p_cell = boost::make_shared<Cell>(p_wild_state, new Gardner1998SbmlCellCycleModel);
            p_cell->SetCellProliferativeType(p_stem_type);

            // Initialise the cell cycle model
            auto p_cc_model = static_cast<Gardner1998SbmlCellCycleModel*>(p_cell->GetCellCycleModel());
            p_cc_model->SetBirthTime(p_simulation_time->GetTime());
            TS_ASSERT_EQUALS(p_cc_model->CanCellTerminallyDifferentiate(), false);

            p_cell->InitialiseCellCycleModel();
            p_cc_model->SetDt(0.01);

            // Update the cell-cycle model so the state variables are different from initial conditions
            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_cc_model->ReadyToDivide();
            }

            C0 = p_cc_model->GetStateVariable("C");
            X0 = p_cc_model->GetStateVariable("X");
            M0 = p_cc_model->GetStateVariable("M");
            Y0 = p_cc_model->GetStateVariable("Y");
            Z0 = p_cc_model->GetStateVariable("Z");

            // Archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* p_cc_arch = p_cc_model;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cc_arch; // TODO: Fix archiving error thrown by this line

            // Deletion of the cell-cycle model is handled by the cell destructor
            SimulationTime::Destroy();
        }

        // Load archive
        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_cc_model;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cc_model;

            double C1 = dynamic_cast<Gardner1998SbmlCellCycleModel*>(p_cc_model)->GetStateVariable("C");
            double X1 = dynamic_cast<Gardner1998SbmlCellCycleModel*>(p_cc_model)->GetStateVariable("X");
            double M1 = dynamic_cast<Gardner1998SbmlCellCycleModel*>(p_cc_model)->GetStateVariable("M");
            double Y1 = dynamic_cast<Gardner1998SbmlCellCycleModel*>(p_cc_model)->GetStateVariable("Y");
            double Z1 = dynamic_cast<Gardner1998SbmlCellCycleModel*>(p_cc_model)->GetStateVariable("Z");
            TS_ASSERT_DELTA(C1, C0, 1e-6);
            TS_ASSERT_DELTA(M1, M0, 1e-6);
            TS_ASSERT_DELTA(X1, X0, 1e-6);
            TS_ASSERT_DELTA(Y1, Y0, 1e-6);
            TS_ASSERT_DELTA(Z1, Z0, 1e-6);

            // Destroy model
            delete p_cc_model;
        }
    }
};

#endif // TESTGARDNER1998SBMLCELLCYCLEMODEL_HPP_
