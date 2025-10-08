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

#ifndef TEST_TAN_2014_SBML_ODE_SYSTEM_HPP_
#define TEST_TAN_2014_SBML_ODE_SYSTEM_HPP_

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
#include "SbmlTestHelpers.hpp"
#include "Timer.hpp"

#include "Tan2014SbmlOdeSystem.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace sth = sbmltesthelpers;

class TestTan2014SbmlOdeSystem : public AbstractCellBasedTestSuite
{
private:
    const unsigned ODE_SIZE = 6u;

    void RunOdeWithSolver(AbstractIvpOdeSolver& rSolver, const std::string solverName)
    {
        try
        {
            Tan2014SbmlOdeSystem ode_system;

            double start_time = 0.0;
            double end_time = 5000.0;
            double max_step = 0.01;
            double sampling_interval = 0.01;

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();

            Timer::Reset();
            OdeSolution ode_solution = rSolver.Solve(&ode_system, initial_conditions, start_time, end_time, max_step, sampling_interval);
            Timer::Print("Tan 2014 (" + solverName + ")");

            // Compare end solutions with Tellurium values
            std::vector<double> end_solution = ode_solution.rGetSolutions().back();
            TS_ASSERT_DELTA(end_solution[0], 80.051972, 1e-3);  // bcat_cm
            TS_ASSERT_DELTA(end_solution[1], 447.597749, 1e-3); // ligand_cm
            TS_ASSERT_DELTA(end_solution[2], 552.402251, 1e-3); // complex_cm
            TS_ASSERT_DELTA(end_solution[3], 62.047741, 1e-3);  // bcat_nu
            TS_ASSERT_DELTA(end_solution[4], 360.013140, 1e-3); // ligand_nu
            TS_ASSERT_DELTA(end_solution[5], 639.986860, 1e-3); // complex_nu

            ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);
            std::vector<double> drag = ode_solution.GetAnyVariable("drag");
            TS_ASSERT_DELTA(drag.back(), 1.0, 1e-3);

            // Compare solution stats with Tellurium values
            for (unsigned i = 0; i < ODE_SIZE; i++)
            {
                std::vector<double> values;
                for (unsigned j = 0; j < ode_solution.rGetSolutions().size(); j++)
                {
                    values.push_back(ode_solution.rGetSolutions()[j][i]);
                }
                double min_val = sth::min(values);
                double max_val = sth::max(values);
                double mean_val = sth::mean(values);
                double std_val = sth::stdev(values);
                double q1_val = sth::quantile(values, 0.25);
                double q2_val = sth::quantile(values, 0.5);
                double q3_val = sth::quantile(values, 0.75);

                if (i == 0) // bcat_cm
                {
                    TS_ASSERT_DELTA(min_val, 45.24753509322006, 1e-2);
                    TS_ASSERT_DELTA(max_val, 80.05197189552682, 1e-2);
                    TS_ASSERT_DELTA(mean_val, 76.79374377296425, 1e-2);
                    TS_ASSERT_DELTA(std_val, 5.332948929206765, 1e-2);
                    TS_ASSERT_DELTA(q1_val, 76.30072625724083, 1e-2);
                    TS_ASSERT_DELTA(q2_val, 79.10049811306001, 1e-2);
                    TS_ASSERT_DELTA(q3_val, 79.8527692249068, 1e-2);
                }
                else if (i == 1) // ligand_cm
                {
                    TS_ASSERT_DELTA(min_val, 447.5977489158307, 1e-2);
                    TS_ASSERT_DELTA(max_val, 581.1170004363729, 1e-2);
                    TS_ASSERT_DELTA(mean_val, 477.46482882069574, 1e-2);
                    TS_ASSERT_DELTA(std_val, 37.585157908239694, 1e-2);
                    TS_ASSERT_DELTA(q1_val, 449.98139359136957, 1e-2);
                    TS_ASSERT_DELTA(q2_val, 459.0489354564536, 1e-2);
                    TS_ASSERT_DELTA(q3_val, 492.86092359129896, 1e-2);
                }
                else if (i == 2) // complex_cm
                {
                    TS_ASSERT_DELTA(min_val, 418.8829995636274, 1e-2);
                    TS_ASSERT_DELTA(max_val, 552.4022510841701, 1e-2);
                    TS_ASSERT_DELTA(mean_val, 522.5351711793052, 1e-2);
                    TS_ASSERT_DELTA(std_val, 37.5851579082397, 1e-2);
                    TS_ASSERT_DELTA(q1_val, 507.13907640870184, 1e-2);
                    TS_ASSERT_DELTA(q2_val, 540.9510645435473, 1e-2);
                    TS_ASSERT_DELTA(q3_val, 550.0186064086312, 1e-2);
                }
                else if (i == 3) // bcat_nu
                {
                    TS_ASSERT_DELTA(min_val, 32.6, 1e-2);
                    TS_ASSERT_DELTA(max_val, 62.04774090925595, 1e-2);
                    TS_ASSERT_DELTA(mean_val, 59.52187943417476, 1e-2);
                    TS_ASSERT_DELTA(std_val, 4.134577756185565, 1e-2);
                    TS_ASSERT_DELTA(q1_val, 59.13996823289651, 1e-2);
                    TS_ASSERT_DELTA(q2_val, 61.310212705453466, 1e-2);
                    TS_ASSERT_DELTA(q3_val, 61.89333036475111, 1e-2);
                }
                else if (i == 4) // ligand_nu
                {
                    TS_ASSERT_DELTA(min_val, 360.0131396545373, 1e-2);
                    TS_ASSERT_DELTA(max_val, 516.8, 1e-2);
                    TS_ASSERT_DELTA(mean_val, 374.1294962358405, 1e-2);
                    TS_ASSERT_DELTA(std_val, 27.472630884273887, 1e-2);
                    TS_ASSERT_DELTA(q1_val, 360.65896071465426, 1e-2);
                    TS_ASSERT_DELTA(q2_val, 363.11746747180143, 1e-2);
                    TS_ASSERT_DELTA(q3_val, 372.7214970204209, 1e-2);
                }
                else if (i == 5) // complex_nu
                {
                    TS_ASSERT_DELTA(min_val, 483.19999999999993, 1e-2);
                    TS_ASSERT_DELTA(max_val, 639.9868603454628, 1e-2);
                    TS_ASSERT_DELTA(mean_val, 625.8705037641597, 1e-2);
                    TS_ASSERT_DELTA(std_val, 27.472630884273933, 1e-2);
                    TS_ASSERT_DELTA(q1_val, 627.2785029795791, 1e-2);
                    TS_ASSERT_DELTA(q2_val, 636.8825325281987, 1e-2);
                    TS_ASSERT_DELTA(q3_val, 639.3410392853458, 1e-2);
                }
            }

            // Compare derived quantity stats with Tellurium values
            {
                double min_val = sth::min(drag);
                double max_val = sth::max(drag);
                double mean_val = sth::mean(drag);
                double std_val = sth::stdev(drag);
                double q1_val = sth::quantile(drag, 0.25);
                double q2_val = sth::quantile(drag, 0.5);
                double q3_val = sth::quantile(drag, 0.75);

                TS_ASSERT_DELTA(min_val, 1.0, 1e-3);
                TS_ASSERT_DELTA(max_val, 1.0, 1e-3);
                TS_ASSERT_DELTA(mean_val, 1.0, 1e-3);
                TS_ASSERT_DELTA(std_val, 0.0, 1e-3);
                TS_ASSERT_DELTA(q1_val, 1.0, 1e-3);
                TS_ASSERT_DELTA(q2_val, 1.0, 1e-3);
                TS_ASSERT_DELTA(q3_val, 1.0, 1e-3);
            }

            // The following code provides nice output for gnuplot
            // use the command
            // plot "tan_2014_cvode.dat" u 1:2
            // or
            // plot "tan_2014_cvode.dat" u 1:3 etc. for the various species...
            // or
            // plot "tan_2014_cvode.dat" u 1:2, "" u 1:3, "" u 1:4 ... for all species

            OutputFileHandler handler("");
            out_stream file = handler.OpenOutputFile("tan_2014_" + solverName + ".dat");
            std::vector<double> times = ode_solution.rGetTimes();
            std::vector<std::vector<double> > solutions = ode_solution.rGetSolutions();
            for (unsigned i = 0; i < solutions.size(); i++)
            {
                (*file) << times[i];
                for (unsigned j = 0; j < solutions[i].size(); j++)
                {
                    (*file) << "\t" << solutions[i][j];
                }
                (*file) << "\t" << drag[i] << "\n"
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
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tan_2014_ode.arch";

        // Save archive
        {
            Tan2014SbmlOdeSystem ode_system;

            // Set state variables to 0...ODE_SIZE-1
            std::vector<double> state_variables;
            for (unsigned i = 0; i < ODE_SIZE; ++i)
            {
                state_variables.push_back(static_cast<double>(i));
            }
            ode_system.SetStateVariables(state_variables);

            // Check initial conditions and state variables
            ode_system.SetDefaultInitialCondition(0, 3.141593);
            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), ODE_SIZE);
            TS_ASSERT_DELTA(initial_conditions[0], 3.141593, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[1], 581.10, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[2], 418.90, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[3], 32.60, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[4], 516.80, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[5], 483.20, 1e-3);

            double var0 = ode_system.GetStateVariable(0);
            double var1 = ode_system.GetStateVariable(1);
            double var2 = ode_system.GetStateVariable(2);
            double var3 = ode_system.GetStateVariable(3);
            double var4 = ode_system.GetStateVariable(4);
            double var5 = ode_system.GetStateVariable(5);

            TS_ASSERT_DELTA(var0, 0.0, 1e-3);
            TS_ASSERT_DELTA(var1, 1.0, 1e-3);
            TS_ASSERT_DELTA(var2, 2.0, 1e-3);
            TS_ASSERT_DELTA(var3, 3.0, 1e-3);
            TS_ASSERT_DELTA(var4, 4.0, 1e-3);
            TS_ASSERT_DELTA(var5, 5.0, 1e-3);

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
            TS_ASSERT_DELTA(initial_conditions[0], 46.60, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[1], 581.10, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[2], 418.90, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[3], 32.60, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[4], 516.80, 1e-3);
            TS_ASSERT_DELTA(initial_conditions[5], 483.20, 1e-3);

            double var0 = p_ode_system->GetStateVariable(0);
            double var1 = p_ode_system->GetStateVariable(1);
            double var2 = p_ode_system->GetStateVariable(2);
            double var3 = p_ode_system->GetStateVariable(3);
            double var4 = p_ode_system->GetStateVariable(4);
            double var5 = p_ode_system->GetStateVariable(5);

            TS_ASSERT_DELTA(var0, 0.0, 1e-3);
            TS_ASSERT_DELTA(var1, 1.0, 1e-3);
            TS_ASSERT_DELTA(var2, 2.0, 1e-3);
            TS_ASSERT_DELTA(var3, 3.0, 1e-3);
            TS_ASSERT_DELTA(var4, 4.0, 1e-3);
            TS_ASSERT_DELTA(var5, 5.0, 1e-3);

            // Tidy up
            delete p_ode_system;
        }
    }

    void TestOdeEquation()
    {
        Tan2014SbmlOdeSystem ode_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(46.6);  // bcat_cm
        initial_conditions.push_back(581.1); // ligand_cm
        initial_conditions.push_back(418.9); // complex_cm
        initial_conditions.push_back(32.6);  // bcat_nu
        initial_conditions.push_back(516.8); // ligand_nu
        initial_conditions.push_back(483.2); // complex_nu

        std::vector<double> derivs(initial_conditions.size(), 0.0);
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Compare derivatives with values from Tellurium
        // (Tellurium's `getRatesOfChange()` is not scaled by compartment)
        const double CytosolMembrane = 1.16;
        const double nucleus = 0.65;
        TS_ASSERT_EQUALS(ode_system.GetParameter("wnt_level"), 0.0);
        TS_ASSERT_DELTA(derivs[0], -1.97469879e+02 / CytosolMembrane, 1e-3); // bcat_cm
        TS_ASSERT_DELTA(derivs[1], 2.73412000e-04 / CytosolMembrane, 1e-6);  // ligand_cm
        TS_ASSERT_DELTA(derivs[2], -2.73412000e-04 / CytosolMembrane, 1e-6); // complex_cm
        TS_ASSERT_DELTA(derivs[3], 1.98105040e+02 / nucleus, 1e-3);          // bcat_nu
        TS_ASSERT_DELTA(derivs[4], 1.04000000e-03 / nucleus, 1e-5);          // ligand_nu
        TS_ASSERT_DELTA(derivs[5], -1.04000000e-03 / nucleus, 1e-5);         // complex_nu

        // Set wnt_level=1 and check derivatives again
        ode_system.SetParameter("wnt_level", 1.0);
        TS_ASSERT_EQUALS(ode_system.GetParameter("wnt_level"), 1.0);
        ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        TS_ASSERT_DELTA(derivs[0], -1.97029323e+02 / CytosolMembrane, 1e-3); // bcat_cm
        TS_ASSERT_DELTA(derivs[1], 2.73412000e-04 / CytosolMembrane, 1e-6);  // ligand_cm
        TS_ASSERT_DELTA(derivs[2], -2.73412000e-04 / CytosolMembrane, 1e-6); // complex_cm
        TS_ASSERT_DELTA(derivs[3], 1.98105040e+02 / nucleus, 1e-3);          // bcat_nu
        TS_ASSERT_DELTA(derivs[4], 1.04000000e-03 / nucleus, 1e-5);          // ligand_nu
        TS_ASSERT_DELTA(derivs[5], -1.04000000e-03 / nucleus, 1e-5);         // complex_nu
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

#endif // TEST_TAN_2014_SBML_ODE_SYSTEM_HPP_
