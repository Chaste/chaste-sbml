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

#ifndef TEST_TAN_2014_SBML_SRN_MODEL_HPP_
#define TEST_TAN_2014_SBML_SRN_MODEL_HPP_

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractSrnModel.hpp"
#include "Cell.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Tan2014SbmlOdeSystem.hpp"
#include "Tan2014SbmlSrnModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestTan2014SbmlSrnModel : public AbstractCellBasedTestSuite
{
public:
    void TestSrnCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(Tan2014SbmlSrnModel srn_model);

        Tan2014SbmlSrnModel* p_srn_model = new Tan2014SbmlSrnModel();

        // Set non-default initial conditions (one per state variable)
        std::vector<double> starter_conditions = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
        p_srn_model->SetInitialConditions(starter_conditions);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        // Model is in minutes; Chaste integrates in hours, so divide the end time by 60.
        double end_time = 10.0 / 60.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_srn_model->SimulateToCurrentTime();
        }

        TS_ASSERT_THROWS_NOTHING(p_srn_model->GetStateVariable("bcat_cm"));
    }

    void TestSrnCreateCopy()
    {
        Tan2014SbmlSrnModel* p_model = new Tan2014SbmlSrnModel;

        // Set ODE system with known state variables
        std::vector<double> state_variables = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        Tan2014SbmlOdeSystem* p_ode_system = new Tan2014SbmlOdeSystem;
        p_ode_system->SetStateVariables(state_variables);
        p_model->SetOdeSystem(p_ode_system);

        // Create a copy
        Tan2014SbmlSrnModel* p_model2 = static_cast<Tan2014SbmlSrnModel*>(p_model->CreateSrnModel());

        // The copy carries the same state
        TS_ASSERT_EQUALS(dynamic_cast<AbstractSbmlSrnModel*>(p_model2)->GetStateVariable("bcat_cm"), 1.0);
        TS_ASSERT_EQUALS(dynamic_cast<AbstractSbmlSrnModel*>(p_model2)->GetStateVariable("complex_nu"), 6.0);

        delete p_model;
        delete p_model2;
    }

    void TestSrnModelOutputParameters()
    {
        OutputFileHandler output_file_handler("TestTan2014SrnOutputParameters", false);

        Tan2014SbmlSrnModel srn_model;
        TS_ASSERT_EQUALS(srn_model.GetIdentifier(), "Tan2014SbmlSrnModel");

        out_stream parameter_file = output_file_handler.OpenOutputFile("tan_2014_srn_results.parameters");
        TS_ASSERT_THROWS_NOTHING(srn_model.OutputSrnModelParameters(parameter_file));
        parameter_file->close();
    }

    void TestSrnArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tan_2014_srn.arch";

        double var0;

        // Save
        {
            // Model is in minutes; Chaste integrates in hours, so divide the end time by 60.
            double end_time = 10.0 / 60.0;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, 100);

            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel;
            AbstractOdeSrnModel* p_srn_model = new Tan2014SbmlSrnModel;

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_cell->SetBirthTime(0.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            var0 = dynamic_cast<AbstractSbmlSrnModel*>(p_srn_model)->GetStateVariable("bcat_cm");

            output_arch << p_srn_model;
            SimulationTime::Destroy();
        }

        // Load
        {
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractSrnModel* p_srn_model;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_srn_model;

            double var1 = dynamic_cast<Tan2014SbmlSrnModel*>(p_srn_model)->GetStateVariable("bcat_cm");
            TS_ASSERT_DELTA(var1, var0, 1e-6);

            delete p_srn_model;
        }
    }
};

#endif // TEST_TAN_2014_SBML_SRN_MODEL_HPP_
