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

#ifndef TEST_GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_
#define TEST_GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_

#include <iostream>
#include <memory>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "OutputFileHandler.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "SbmlTestHelpers.hpp"
#include "SmartPointers.hpp"
#include "Timer.hpp"

#include "Goldbeter1991SbmlOdeSystem.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace sth = sbmltesthelpers;

class TestGoldbeter1991SbmlOdeSystem : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 3u;

    std::vector<double> default_initial_conditions = {
        0.01, // C
        0.01, // M
        0.01, // X
    };

    void RunOdeWithSolver(AbstractIvpOdeSolver& rSolver, const std::string solverName)
    {
        try
        {
            // Solve system using solver
            Goldbeter1991SbmlOdeSystem ode_system;

            double start_time = 0.0;
            double end_time = 100.0;
            double max_step = 0.01;
            double sampling_interval = 0.1;

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            OdeSolution solutions;

            Timer::Reset();
            solutions = rSolver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
            Timer::Print("Goldbeter 1991 (" + solverName + ")");

            // No events, so ODE shouldn't have stopped
            TS_ASSERT_EQUALS(rSolver.StoppingEventOccurred(), false);

            // Check final solution with Tellurium values
            std::vector<double> end_solution = solutions.rGetSolutions().back();
            TS_ASSERT_DELTA(end_solution[0], 0.547050, 1e-3);
            TS_ASSERT_DELTA(end_solution[1], 0.293640, 1e-3);
            TS_ASSERT_DELTA(end_solution[2], 0.006785, 1e-3);

            // Exports results to csv
            sth::export_csv("goldbeter_1991_" + solverName + ".csv", solutions, ode_system);
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
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "goldbeter_1991_ode.arch";

        // Save archive
        {
            Goldbeter1991SbmlOdeSystem ode_system;

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

            std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
            TSM_ASSERT_DELTA(var_names[0].c_str(), initial_conditions[0], 3.141593, 1e-6);
            for (unsigned i = 1; i < ODE_SIZE; i++)
            {
                const char* var_cname = var_names[i].c_str();

                // Check initial condition
                TSM_ASSERT_DELTA(var_cname, initial_conditions[i], default_initial_conditions[i], 1e-6);

                // Check state variable
                double var = ode_system.GetStateVariable(i);
                TSM_ASSERT_DELTA(var_cname, var, static_cast<double>(i), 1e-6);
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

            std::vector<std::string> var_names = p_ode_system->rGetStateVariableNames();
            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                const char* var_cname = var_names[i].c_str();

                // Check initial condition
                TSM_ASSERT_DELTA(var_cname, initial_conditions[i], default_initial_conditions[i], 1e-6);

                // Check state variable
                double var = p_ode_system->GetStateVariable(i);
                TSM_ASSERT_DELTA(var_cname, var, static_cast<double>(i), 1e-6);
            }

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        Goldbeter1991SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> derivs(ODE_SIZE);
        ode_system.EvaluateYDerivatives(time, default_initial_conditions, derivs);

        // Compare derivatives with values from Tellurium
        std::vector<double> derivs_expected = {
            0.02406667,  // C
            -0.94147207, // M
            -0.32338358, // X
        };

        std::vector<std::string> var_names = ode_system.rGetStateVariableNames();
        for (unsigned i = 0; i < ODE_SIZE; i++)
        {
            TSM_ASSERT_DELTA(var_names[i].c_str(), derivs[i], derivs_expected[i], 1e-6);
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

#endif // TEST_GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_
