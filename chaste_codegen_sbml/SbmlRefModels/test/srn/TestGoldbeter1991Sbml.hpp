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

#ifndef TESTGOLDBETER1991SBML_HPP_
#define TESTGOLDBETER1991SBML_HPP_

#include <fstream>
#include <iostream>
#include <memory>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "Cell.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "Timer.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Goldbeter1991SbmlOdeSystemAndSrnModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestGoldbeter1991Sbml : public AbstractCellBasedTestSuite
{
public:
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "goldbeter_1991_ode.arch";

        {
            std::vector<double> state_variables;
            state_variables.push_back(3.0);
            state_variables.push_back(4.0);
            state_variables.push_back(5.0);

            Goldbeter1991SbmlOdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 3.2500, 1e-4);

            double var1 = ode_system.GetStateVariable(0);
            double var2 = ode_system.GetStateVariable(1);
            double var3 = ode_system.GetStateVariable(2);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.01, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 0.01, 1e-4);

            double var1 = p_ode_system->GetStateVariable(0);
            double var2 = p_ode_system->GetStateVariable(1);
            double var3 = p_ode_system->GetStateVariable(2);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        Goldbeter1991SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.01);
        initial_conditions.push_back(0.01);
        initial_conditions.push_back(0.01);

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct
        TS_ASSERT_DELTA(derivs[0], 0.0240, 1e-4);
        TS_ASSERT_DELTA(derivs[1], -0.9414, 1e-4);
        TS_ASSERT_DELTA(derivs[2], -0.3233, 1e-4);
    }

    void TestOdeWithChasteSolver()
    {
        try
        {
            Goldbeter1991SbmlOdeSystem ode_system;

            // Solve system using RK4 solver

            double dt = 0.0001;

            // RK4 solver solution worked out
            RungeKutta4IvpOdeSolver rk4_solver;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            Timer::Reset();
            OdeSolution solutions = rk4_solver.Solve(&ode_system, state_variables, 0.0, 100.0, dt, dt);
            Timer::Print("1. Goldbeter RK4");

            unsigned end = solutions.rGetSolutions().size() - 1;

            //  Decent results - checked with numpy # [ 0.54706214  0.29369527  0.00678837]
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.5470, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.2936, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 0.0067, 1e-3);

            // The following code provides nice output for gnuplot
            // use the command
            // plot "goldbeter_1991.dat" u 1:2
            // or
            // plot "goldbeter_1991.dat" u 1:3 etc. for the various species...
            // or
            // plot "goldbeter_1991.dat" u 1:2, "" u 1:3, "" u 1:4 ... for all species

            // OutputFileHandler handler("");
            // out_stream file = handler.OpenOutputFile("goldbeter_1991.dat");
            // for (unsigned i = 0; i < solutions.rGetSolutions().size(); i++)
            // {
            //   (*file) << solutions.rGetTimes()[i];
            //   for (unsigned j = 0; j < solutions.rGetSolutions()[i].size(); j++)
            //   {
            //     (*file) << "\t" << solutions.rGetSolutions()[i][j];
            //   }
            //   (*file) << "\n"
            //           << std::flush;
            // }
            // file->close();
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

    void TestOdeWithCvodeSolver()
    {
        try
        {
            Goldbeter1991SbmlOdeSystem ode_system;

            double end_time = 100;
            double h_value = 0.01;

            CvodeAdaptor solver;
            OdeSolution solutions;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            Timer::Reset();
            solutions = solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 0.1);
            Timer::Print("1. Goldbeter CVODE");

            unsigned end = solutions.rGetSolutions().size() - 1;

            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.5470, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.2936, 1e-3);
            TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 0.0067, 1e-3);
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

    void TestSrnArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "goldbeter_1991_srn.arch";

        double C0, M0, X0;

        // Create an output archive
        {
            double end_time = 10.0;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, 100);

            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractSrnModel* p_srn_model = new Goldbeter1991SbmlSrnModel;

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            // We must create a cell to be able to initialise the cell SRN model's ODE system
            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_cell->SetBirthTime(0.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Now update the SRN so the state variables are different from ICS
            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            C0 = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_srn_model)->GetStateVariable("C");
            M0 = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_srn_model)->GetStateVariable("M");
            X0 = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_srn_model)->GetStateVariable("X");

            output_arch << p_srn_model;

            // Note that here, deletion of the cell-cycle model and SRN is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractSrnModel* p_srn_model;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_srn_model;

            double C1 = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_srn_model)->GetStateVariable("C");
            double M1 = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_srn_model)->GetStateVariable("M");
            double X1 = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_srn_model)->GetStateVariable("X");
            TS_ASSERT_DELTA(C1, C0, 1e-6);
            TS_ASSERT_DELTA(M1, M0, 1e-6);
            TS_ASSERT_DELTA(X1, X0, 1e-6);

            // Destroy model
            delete p_srn_model;
        }
    }

    void TestSrnCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(Goldbeter1991SbmlSrnModel srn_model);

        Goldbeter1991SbmlSrnModel* p_srn_model = new Goldbeter1991SbmlSrnModel();

        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.6);
        starter_conditions.push_back(0.7);
        p_srn_model->SetInitialConditions(starter_conditions);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->GetCellData()->SetItem("mean delta", 1.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        double end_time = 10.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_srn_model->SimulateToCurrentTime();
        }
    }

    void TestSrnCreateCopy()
    {
        // Test with Goldbeter1991SbmlSrnModel
        Goldbeter1991SbmlSrnModel* p_model = new Goldbeter1991SbmlSrnModel;

        // Set ODE system
        std::vector<double> state_variables;
        state_variables.push_back(2.0);
        state_variables.push_back(3.0);
        state_variables.push_back(4.0);

        p_model->SetOdeSystem(new Goldbeter1991SbmlOdeSystem(state_variables));

        // Create a copy
        Goldbeter1991SbmlSrnModel* p_model2 = static_cast<Goldbeter1991SbmlSrnModel*>(p_model->CreateSrnModel());

        // Check correct initializations
        double C = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_model2)->GetStateVariable("C");
        double M = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_model2)->GetStateVariable("M");
        double X = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_model2)->GetStateVariable("X");
        TS_ASSERT_EQUALS(C, 2.0);
        TS_ASSERT_EQUALS(M, 3.0);
        TS_ASSERT_EQUALS(X, 4.0);

        // Destroy models
        delete p_model;
        delete p_model2;
    }

    void TestSrnModelOutputParameters()
    {
        std::string output_directory = "TestSrnModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        Goldbeter1991SbmlSrnModel srn_model;

        TS_ASSERT_EQUALS(srn_model.GetIdentifier(), "SbmlSrnWrapperModel-Goldbeter1991SbmlOdeSystem-3");

        out_stream parameter_file = output_file_handler.OpenOutputFile("goldbeter_1991_srn_results.parameters");
        srn_model.OutputSrnModelParameters(parameter_file);
        parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("goldbeter_1991_srn_results.parameters");
            FileFinder reference("projects/SbmlRefModels/test/data/TestGoldbeter1991SBML/goldbeter_1991_srn_results.parameters",
                                 RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestSteadyStateSimulation()
    {
        // Keep running until we reach steady state
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        // run until 100, with dt=0.001
        double end_time = 100;
        double dt = 0.001;
        unsigned num_steps = (unsigned)end_time / dt;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps + 1);

        // Create a cell-cycle model
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SbmlSrnModel* p_srn_model = new Goldbeter1991SbmlSrnModel();

        CellPtr p_tn_cell(new Cell(p_healthy_state, p_cell_model, p_srn_model, false, CellPropertyCollection()));
        p_tn_cell->SetCellProliferativeType(p_transit_type);
        p_tn_cell->InitialiseCellCycleModel();
        p_tn_cell->InitialiseSrnModel();

        // Run the cell simulation until end_time
        while (!p_simulation_time->IsFinished())
        {
            p_simulation_time->IncrementTimeOneStep();
            if (p_tn_cell->ReadyToDivide())
            {
                p_tn_cell->Divide();
            }
        }

        // Direct access to state variables
        double C = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_tn_cell->GetSrnModel())->GetStateVariable("C");
        double M = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_tn_cell->GetSrnModel())->GetStateVariable("M");
        double X = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_tn_cell->GetSrnModel())->GetStateVariable("X");

        TS_ASSERT_DELTA(C, 0.5470, 1e-2);
        TS_ASSERT_DELTA(M, 0.2936, 1e-2);
        TS_ASSERT_DELTA(X, 0.0067, 1e-3);

        // Indirect access to state vector
        C = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_tn_cell->GetSrnModel())->GetProteinConcentrations()[0];
        M = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_tn_cell->GetSrnModel())->GetProteinConcentrations()[1];
        X = dynamic_cast<Goldbeter1991SbmlSrnModel*>(p_tn_cell->GetSrnModel())->GetProteinConcentrations()[2];

        TS_ASSERT_DELTA(C, 0.5470, 1e-2);
        TS_ASSERT_DELTA(M, 0.2936, 1e-2);
        TS_ASSERT_DELTA(X, 0.0067, 1e-3);
    }
};

#endif // TESTGOLDBETER1991SBML_HPP_
