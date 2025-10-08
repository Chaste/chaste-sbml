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

#ifndef TEST_VAN_LEEUWEN_2007_NON_DIM_SBML_ODE_SYSTEM_HPP_
#define TEST_VAN_LEEUWEN_2007_NON_DIM_SBML_ODE_SYSTEM_HPP_

#include <iostream>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "CvodeAdaptor.hpp"
#include "Debug.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "Timer.hpp"

#include "VanLeeuwen2007NonDimSbmlOdeSystem.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestVanLeeuwen2007NonDimSbmlOdeSystem : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 11u;

    void RunOdeWithSolver(AbstractIvpOdeSolver& rSolver, const std::string solverName)
    {
        try
        {
            VanLeeuwen2007NonDimSbmlOdeSystem ode_system;

            double dt = 0.01;
            double end_time = 1000.0;

            std::vector<double> state_variables = ode_system.GetInitialConditions();

            Timer::Reset();
            OdeSolution ode_solution = rSolver.Solve(&ode_system, state_variables, 0.0, end_time, dt, dt);
            Timer::Print("VanLeeuwen 2007 Non-dim (" + solverName + ")");

            // Compare end solutions with Tellurium values
            std::vector<double> end_solution = ode_solution.rGetSolutions().back();
            TS_ASSERT_DELTA(end_solution[0], 0.266667, 1e-2);   // X
            TS_ASSERT_DELTA(end_solution[1], 2.666669, 1e-2);   // D
            TS_ASSERT_DELTA(end_solution[2], 76.190309, 1e-2);  // C_o
            TS_ASSERT_DELTA(end_solution[3], 13.473635, 1e-2);  // C_u
            TS_ASSERT_DELTA(end_solution[4], 0.0, 1e-6);        // C_c
            TS_ASSERT_DELTA(end_solution[5], 299.964754, 1e-2); // A
            TS_ASSERT_DELTA(end_solution[6], 544.152296, 1e-2); // C_A
            TS_ASSERT_DELTA(end_solution[7], 749.997274, 1e-2); // T
            TS_ASSERT_DELTA(end_solution[8], 76.190010, 1e-2);  // C_oT
            TS_ASSERT_DELTA(end_solution[9], 0.0, 1e-6);        // C_cT
            TS_ASSERT_DELTA(end_solution[10], 36.180879, 1e-2); // Y

            // The following code provides nice output for gnuplot
            // use the command
            // plot "vanleeuwen_2007_nondim_cvode.dat" u 1:2
            // or
            // plot "vanleeuwen_2007_nondim_cvode.dat" u 1:3 etc. for the various species...
            // or
            // plot "vanleeuwen_2007_nondim_cvode.dat" u 1:2, "" u 1:3, "" u 1:4 ... for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("vanleeuwen_2007_nondim_" + solverName + ".dat");
            std::vector<double> times = ode_solution.rGetTimes();
            std::vector<std::vector<double> > solutions = ode_solution.rGetSolutions();
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
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "van_leeuwen_2007_nondim_ode.arch";

        // Save archive
        {
            VanLeeuwen2007NonDimSbmlOdeSystem ode_system;

            // Set state variables to 0...ODE_SIZE-1
            std::vector<double> state_variables;
            state_variables.push_back(0.0);  // X
            state_variables.push_back(1.0);  // D
            state_variables.push_back(2.0);  // C_o
            state_variables.push_back(3.0);  // C_u
            state_variables.push_back(4.0);  // C_c
            state_variables.push_back(5.0);  // A
            state_variables.push_back(6.0);  // C_A
            state_variables.push_back(7.0);  // T
            state_variables.push_back(8.0);  // C_oT
            state_variables.push_back(9.0);  // C_cT
            state_variables.push_back(10.0); // Y
            ode_system.SetStateVariables(state_variables);

            // Check initial conditions and state variables
            ode_system.SetDefaultInitialCondition(0, 3.141593);
            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);
            TS_ASSERT_DELTA(initial_conditions[0], 3.141593, 1e-3);  // X
            TS_ASSERT_DELTA(initial_conditions[1], 2.68, 1e-3);   // D
            TS_ASSERT_DELTA(initial_conditions[2], 76.2, 1e-3); // C_o
            TS_ASSERT_DELTA(initial_conditions[3], 13.5, 1e-3);   // C_u
            TS_ASSERT_DELTA(initial_conditions[4], 0.0, 1e-3);    // C_c
            TS_ASSERT_DELTA(initial_conditions[5], 300.0, 1e-3);  // A
            TS_ASSERT_DELTA(initial_conditions[6], 544.2, 1e-3);  // C_A
            TS_ASSERT_DELTA(initial_conditions[7], 750.0, 1e-3);  // T
            TS_ASSERT_DELTA(initial_conditions[8], 76.2, 1e-3);   // C_oT
            TS_ASSERT_DELTA(initial_conditions[9], 0.0, 1e-3);    // C_cT
            TS_ASSERT_DELTA(initial_conditions[10], 36.0, 1e-3);  // Y

            double var0 = ode_system.GetStateVariable(0);   // X
            double var1 = ode_system.GetStateVariable(1);   // D
            double var2 = ode_system.GetStateVariable(2);   // C_o
            double var3 = ode_system.GetStateVariable(3);   // C_u
            double var4 = ode_system.GetStateVariable(4);   // C_c
            double var5 = ode_system.GetStateVariable(5);   // A
            double var6 = ode_system.GetStateVariable(6);   // C_A
            double var7 = ode_system.GetStateVariable(7);   // T
            double var8 = ode_system.GetStateVariable(8);   // C_oT
            double var9 = ode_system.GetStateVariable(9);   // C_cT
            double var10 = ode_system.GetStateVariable(10); // Y

            TS_ASSERT_DELTA(var0, 0.0, 1e-3);   // X
            TS_ASSERT_DELTA(var1, 1.0, 1e-3);   // D
            TS_ASSERT_DELTA(var2, 2.0, 1e-3);   // C_o
            TS_ASSERT_DELTA(var3, 3.0, 1e-3);   // C_u
            TS_ASSERT_DELTA(var4, 4.0, 1e-3);   // C_c
            TS_ASSERT_DELTA(var5, 5.0, 1e-3);   // A
            TS_ASSERT_DELTA(var6, 6.0, 1e-3);   // C_A
            TS_ASSERT_DELTA(var7, 7.0, 1e-3);   // T
            TS_ASSERT_DELTA(var8, 8.0, 1e-3);   // C_oT
            TS_ASSERT_DELTA(var9, 9.0, 1e-3);   // C_cT
            TS_ASSERT_DELTA(var10, 10.0, 1e-3); // Y

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        // Load archive
        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);
            TS_ASSERT_DELTA(initial_conditions[0], 0.268, 1e-3); // X
            TS_ASSERT_DELTA(initial_conditions[1], 2.68, 1e-3);  // D
            TS_ASSERT_DELTA(initial_conditions[2], 76.2, 1e-3);  // C_o
            TS_ASSERT_DELTA(initial_conditions[3], 13.5, 1e-3);  // C_u
            TS_ASSERT_DELTA(initial_conditions[4], 0.0, 1e-3);   // C_c
            TS_ASSERT_DELTA(initial_conditions[5], 300.0, 1e-3); // A
            TS_ASSERT_DELTA(initial_conditions[6], 544.2, 1e-3); // C_A
            TS_ASSERT_DELTA(initial_conditions[7], 750.0, 1e-3); // T
            TS_ASSERT_DELTA(initial_conditions[8], 76.2, 1e-3);  // C_oT
            TS_ASSERT_DELTA(initial_conditions[9], 0.0, 1e-3);   // C_cT
            TS_ASSERT_DELTA(initial_conditions[10], 36.0, 1e-3); // Y

            double var0 = p_ode_system->GetStateVariable(0);   // X
            double var1 = p_ode_system->GetStateVariable(1);   // D
            double var2 = p_ode_system->GetStateVariable(2);   // C_o
            double var3 = p_ode_system->GetStateVariable(3);   // C_u
            double var4 = p_ode_system->GetStateVariable(4);   // C_c
            double var5 = p_ode_system->GetStateVariable(5);   // A
            double var6 = p_ode_system->GetStateVariable(6);   // C_A
            double var7 = p_ode_system->GetStateVariable(7);   // T
            double var8 = p_ode_system->GetStateVariable(8);   // C_oT
            double var9 = p_ode_system->GetStateVariable(9);   // C_cT
            double var10 = p_ode_system->GetStateVariable(10); // Y

            TS_ASSERT_DELTA(var0, 0.0, 1e-3);   // X
            TS_ASSERT_DELTA(var1, 1.0, 1e-3);   // D
            TS_ASSERT_DELTA(var2, 2.0, 1e-3);   // C_o
            TS_ASSERT_DELTA(var3, 3.0, 1e-3);   // C_u
            TS_ASSERT_DELTA(var4, 4.0, 1e-3);   // C_c
            TS_ASSERT_DELTA(var5, 5.0, 1e-3);   // A
            TS_ASSERT_DELTA(var6, 6.0, 1e-3);   // C_A
            TS_ASSERT_DELTA(var7, 7.0, 1e-3);   // T
            TS_ASSERT_DELTA(var8, 8.0, 1e-3);   // C_oT
            TS_ASSERT_DELTA(var9, 9.0, 1e-3);   // C_cT
            TS_ASSERT_DELTA(var10, 10.0, 1e-3); // Y

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        VanLeeuwen2007NonDimSbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.268); // X
        initial_conditions.push_back(2.68);  // D
        initial_conditions.push_back(76.2);  // C_o
        initial_conditions.push_back(13.5);  // C_u
        initial_conditions.push_back(0.0);   // C_c
        initial_conditions.push_back(300.0); // A
        initial_conditions.push_back(544.2); // C_A
        initial_conditions.push_back(750.0); // T
        initial_conditions.push_back(76.2);  // C_oT
        initial_conditions.push_back(0.0);   // C_cT
        initial_conditions.push_back(36.0);  // Y

        std::vector<double> derivs(initial_conditions.size());
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Compare derivatives with values from Tellurium (wnt-level=0, gamma1=1, gamma2=1)
        TS_ASSERT_EQUALS(ode_system.GetParameter("wnt_level"), 0.0);
        TS_ASSERT_EQUALS(ode_system.GetParameter("gamma1"), 1.0);
        TS_ASSERT_EQUALS(ode_system.GetParameter("gamma2"), 1.0);
        TS_ASSERT_DELTA(derivs[0], -2.66666667e-04, 1e-6); // X
        TS_ASSERT_DELTA(derivs[1], -1.04083409e-16, 1e-6); // D
        TS_ASSERT_DELTA(derivs[2], -4.44116711e-02, 1e-4); // C_o
        TS_ASSERT_DELTA(derivs[3], 2.81167109e-03, 1e-5);  // C_u
        TS_ASSERT_DELTA(derivs[4], 0.0, 1e-6);             // C_c
        TS_ASSERT_DELTA(derivs[5], -4.00000000e-02, 1e-4); // A
        TS_ASSERT_DELTA(derivs[6], 4.00000000e-02, 1e-6);  // C_A
        TS_ASSERT_DELTA(derivs[7], 1.99174011e-13, 1e-6);  // T
        TS_ASSERT_DELTA(derivs[8], -1.98951966e-13, 1e-6); // C_oT
        TS_ASSERT_DELTA(derivs[9], 0.0, 1e-5);             // C_cT
        TS_ASSERT_DELTA(derivs[10], 3.44118767e-04, 1e-6); // Y

        // Set wnt_level=1 and check derivatives again
        ode_system.SetParameter("wnt_level", 1.0);
        TS_ASSERT_EQUALS(ode_system.GetParameter("wnt_level"), 1.0);
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        TS_ASSERT_DELTA(derivs[0], -5.38666667e-02, 1e-4); // X
        TS_ASSERT_DELTA(derivs[1], -3.57333333e-02, 1e-4); // D
        TS_ASSERT_DELTA(derivs[2], -4.44116711e-02, 1e-4); // C_o
        TS_ASSERT_DELTA(derivs[3], 2.81167109e-03, 1e-5);  // C_u
        TS_ASSERT_DELTA(derivs[4], 0.0, 1e-6);             // C_c
        TS_ASSERT_DELTA(derivs[5], -4.00000000e-02, 1e-4); // A
        TS_ASSERT_DELTA(derivs[6], 4.00000000e-02, 1e-4);  // C_A
        TS_ASSERT_DELTA(derivs[7], 1.99174011e-13, 1e-6);  // T
        TS_ASSERT_DELTA(derivs[8], -1.98951966e-13, 1e-6); // C_oT
        TS_ASSERT_DELTA(derivs[9], 0.0, 1e-6);             // C_cT
        TS_ASSERT_DELTA(derivs[10], 3.44118767e-04, 1e-6); // Y
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
        CvodeAdaptor solver;
        RunOdeWithSolver(solver, "cvode");
    }
};

#endif // TEST_VAN_LEEUWEN_2007_NON_DIM_SBML_ODE_SYSTEM_HPP_
