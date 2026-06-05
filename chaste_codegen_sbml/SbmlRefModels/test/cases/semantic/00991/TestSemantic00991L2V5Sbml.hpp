#ifndef TEST_SEMANTIC_00991_L2_V5_SBML_HPP_
#define TEST_SEMANTIC_00991_L2_V5_SBML_HPP_

#include <cmath>
#include <string>
#include <vector>

#include <cxxtest/TestSuite.h>

#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
#include "SbmlTestHelpers.hpp"

#include "Semantic00991L2V5SbmlOdeSystem.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace sth = sbmltesthelpers;

class Test : public CxxTest::TestSuite
{
public:
    void TestOdeSystem()
    {
        try
        {
            Semantic00991L2V5SbmlOdeSystem ode_system;
            OdeSolution ode_solution;

            CvodeAdaptor solver;
            solver.CheckForStoppingEvents();

            // Settings
            double start = 0;
            double duration = 2;
            double steps = 20;

            double end = start + duration;
            double sampling = duration / steps;
            double timestep = sampling / 10.0;

            // Solve
            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            ode_solution = solver.Solve(&ode_system, initial_conditions, start, end, timestep, sampling);

            while (solver.StoppingEventOccurred() && ode_solution.rGetTimes().back() < end)
            {
                // Get state at stopping time
                std::vector<double> state_at_event = ode_solution.rGetSolutions().back();
                double time_at_event = ode_solution.rGetTimes().back();

                // Update ODE system
                ode_system.SetStateVariables(state_at_event);
                ode_system.AdjustParameters(time_at_event);

                // Continue solve
                start = time_at_event;
                end = start + (duration - (time_at_event - 0));
                initial_conditions = ode_system.GetStateVariables();

                // Force CVODE to reinitialize so BDF history from before the event
                // does not pollute the first integration step after the restart.
                solver.SetForceReset(true);
                OdeSolution next_solution = solver.Solve(&ode_system, initial_conditions, start, end, timestep, sampling);
                solver.SetForceReset(false);

                // Append new solution to existing solution
                sth::AppendOdeSolution(&ode_solution, &next_solution);
            }

            ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);

            // Expected results
            std::vector<std::string> expected_result_columns = { "time", "X" };
            std::set<std::string> expected_amounts{ "X" };

            std::vector<std::vector<double> > expected_result_data = {
                { 0, 1 },
                { 0.1, 1.1 },
                { 0.2, 1.2 },
                { 0.3, 1.3 },
                { 0.4, 1.4 },
                { 0.5, 1.5 },
                { 0.6, 1.6 },
                { 0.7, 1.7 },
                { 0.8, 1.8 },
                { 0.9, 1.9 },
                { 1, 2 },
                { 1.1, 2.2 },
                { 1.2, 2.4 },
                { 1.3, 2.6 },
                { 1.4, 2.8 },
                { 1.5, 3 },
                { 1.6, 3.2 },
                { 1.7, 3.4 },
                { 1.8, 3.6 },
                { 1.9, 3.8 },
                { 2, 4 }
            };

            // Check variable values
            double tol_absolute = 0.0001 * 10.0; // TODO: review tolerance values
            double tol_relative = 0.0001 * 10.0;

            for (unsigned j = 1; j < expected_result_columns.size(); j++)
            {
                std::string var_name = expected_result_columns[j];
                TSM_ASSERT_EQUALS(var_name.c_str(), ode_system.HasAnyVariable(var_name), true);

                if (expected_amounts.find(var_name) != expected_amounts.end())
                {
                    // Use amount variable
                    if (ode_system.HasAnyVariable("amt__" + var_name))
                    {
                        var_name = "amt__" + var_name;
                    }
                }

                std::vector<double> values = ode_solution.GetAnyVariable(var_name);
                TS_ASSERT_EQUALS(values.size(), expected_result_data.size());
                for (unsigned i = 0; i < expected_result_data.size(); i++)
                {
                    double delta = std::abs(expected_result_data[i][j] - values[i]);
                    double tol = tol_absolute + tol_relative * std::abs(expected_result_data[i][j]);
                    tol = std::max(tol, 1e-6); // Set minimum tolerance to avoid false failures
                    std::string msg(sth::ToString(values[i]) + " vs " + sth::ToString(expected_result_data[i][j])
                                    + " at " + sth::ToString(ode_solution.rGetTimes()[i], 3) + " for " + var_name);
                    TSM_ASSERT_LESS_THAN_EQUALS(msg.c_str(), delta, tol);
                }
            }

            // Exports results to csv
            sth::ExportCsv("Semantic00991L2V5SbmlOdeSystem.csv", ode_solution, ode_system);
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
};

#endif // TEST_SEMANTIC_00991_L2_V5_SBML_HPP_