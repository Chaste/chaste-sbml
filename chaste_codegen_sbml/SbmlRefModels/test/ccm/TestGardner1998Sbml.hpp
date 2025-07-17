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

#ifndef TESTGARDNER1998SBML_HPP_
#define TESTGARDNER1998SBML_HPP_

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

#include "Gardner1998SbmlOdeSystemAndCellCycleModel.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace st = sbmltest;

class TestGardner1998Sbml : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 5u;

    void RunOdeWithSolver(AbstractIvpOdeSolver &rSolver, const std::string solverName)
    {
        try
        {
            Gardner1998SbmlOdeSystem ode_system;

            double dt = 0.01;
            double end_time = 200.0;

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();

            Timer::Reset();
            OdeSolution ode_solution = rSolver.Solve(&ode_system, initial_conditions, 0.0, end_time, dt, dt);
            Timer::Print("Gardner 1998 (" + solverName + ")");

            // Compare end solutions with Tellurium values
            std::vector<double> end_solution = ode_solution.rGetSolutions().back();
            TS_ASSERT_DELTA(end_solution[0], 0.199671, 1e-3); // C
            TS_ASSERT_DELTA(end_solution[1], 0.144840, 1e-2); // X
            TS_ASSERT_DELTA(end_solution[2], 0.332992, 1e-2); // M
            TS_ASSERT_DELTA(end_solution[3], 3.960409, 1e-2); // Y
            TS_ASSERT_DELTA(end_solution[4], 0.405661, 1e-2); // Z

            // The following code provides nice output for gnuplot
            // use the command
            // plot "gardner_1998_cvode.dat" u 1:2
            // or
            // plot "gardner_1998_cvode.dat" u 1:3 etc. for the various species...
            // or
            // plot "gardner_1998_cvode.dat" u 1:2, "" u 1:3, "" u 1:4 ... for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("gardner_1998_" + solverName + ".dat");
            std::vector<double> times = ode_solution.rGetTimes();
            std::vector<std::vector<double>> solutions = ode_solution.rGetSolutions();
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
    void TestOdeArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "gardner_1998_ode.arch";

        {
            std::vector<double> state_variables = {0.0, 1.0, 2.0, 3.0, 4.0};

            Gardner1998SbmlOdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);
            TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-6);  // C
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6);  // X
            TS_ASSERT_DELTA(initial_conditions[2], 3.25, 1e-6); // M
            TS_ASSERT_DELTA(initial_conditions[3], 1.0, 1e-6);  // Y
            TS_ASSERT_DELTA(initial_conditions[4], 1.0, 1e-6);  // Z

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
            TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-6); // C
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6); // X
            TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-6); // M
            TS_ASSERT_DELTA(initial_conditions[3], 1.0, 1e-6); // Y
            TS_ASSERT_DELTA(initial_conditions[4], 1.0, 1e-6); // Z

            double var0 = p_ode_system->GetStateVariable(0); // C
            double var1 = p_ode_system->GetStateVariable(1); // X
            double var2 = p_ode_system->GetStateVariable(2); // M
            double var3 = p_ode_system->GetStateVariable(3); // Y
            double var4 = p_ode_system->GetStateVariable(4); // Z

            TS_ASSERT_DELTA(var0, 0.0, 1e-6); // C
            TS_ASSERT_DELTA(var1, 1.0, 1e-6); // X
            TS_ASSERT_DELTA(var2, 2.0, 1e-6); // M
            TS_ASSERT_DELTA(var3, 3.0, 1e-6); // Y
            TS_ASSERT_DELTA(var4, 4.0, 1e-6); // Z

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        Gardner1998SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.0); // C
        initial_conditions.push_back(0.0); // X
        initial_conditions.push_back(0.0); // M
        initial_conditions.push_back(1.0); // Y
        initial_conditions.push_back(1.0); // Z

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Compare derivatives with values from Tellurium
        TS_ASSERT_DELTA(derivs[0], 0.155, 1e-3);  // C
        TS_ASSERT_DELTA(derivs[1], 0.0, 1e-6);    // X
        TS_ASSERT_DELTA(derivs[2], 0.0, 1e-6);    // M
        TS_ASSERT_DELTA(derivs[3], 0.202, 1e-3);  // Y
        TS_ASSERT_DELTA(derivs[4], -0.057, 1e-3); // Z
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

#endif // TESTGARDNER1998SBML_HPP_
