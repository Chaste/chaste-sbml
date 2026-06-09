#ifndef {{ test_header_guard }}
#define {{ test_header_guard }}

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <cxxtest/TestSuite.h>

#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
#include "SbmlTestOdeSolution.hpp"
#include "SbmlTestHelpers.hpp"

#include "{{ ode_class_name }}.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

namespace sth = sbmltesthelpers;

class Test{{ model_name }} : public CxxTest::TestSuite
{
public:
    void TestOdeSystem()
    {
        try
        {
            {{ ode_class_name }} ode_system;
            SbmlTestOdeSolution ode_solution;

            CvodeAdaptor solver;
            solver.CheckForStoppingEvents();

            // Settings
            double start = {{ test_settings["start"] }};
            double duration = {{ test_settings["duration"] }};
            double steps = {{ test_settings["steps"] }};

            double end = start + duration;
            double sampling = duration / steps;
            double timestep = sampling / 10.0;

            // Solve. Each segment is appended via AppendSegment, which also records the
            // current parameters for every step so that parameters changed by events are
            // time-resolved (plain OdeSolution keeps only a single parameter snapshot).
            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            OdeSolution segment = solver.Solve(&ode_system, initial_conditions, start, end, timestep, sampling);
            ode_solution.AppendSegment(segment, &ode_system);

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
                ode_solution.AppendSegment(next_solution, &ode_system);
            }

            ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);

            // Expected results
            std::vector<std::string> expected_result_columns = { {{ test_result_columns }} };
            std::set<std::string> expected_amounts{ {{ test_amounts }} };

            std::vector<std::vector<double> > expected_result_data = {
                {{ test_result_data }}
            };

            // Check variable values
            double tol_absolute = {{ test_settings["absolute"] }} * 10.0; // TODO: review tolerance values
            double tol_relative = {{ test_settings["relative"] }} * 10.0;

            for (unsigned j = 1; j < expected_result_columns.size(); j++)
            {
                std::string var_name = expected_result_columns[j];
                TSM_ASSERT_EQUALS(var_name.c_str(), ode_system.HasAnyVariable(var_name), true);

                if (expected_amounts.find(var_name) != expected_amounts.end())
                {
                    // Use amount variable
                    if (ode_system.HasAnyVariable("{{ AMOUNT_PREFIX }}{{ PREFIX_SEP }}" + var_name))
                    {
                        var_name = "{{ AMOUNT_PREFIX }}{{ PREFIX_SEP }}" + var_name;
                    }
                }

                // Parameters changed by events vary in time, but OdeSolution stores only a
                // single parameter snapshot, so read those from the per-step record instead.
                const std::vector<std::string>& param_names = ode_system.rGetParameterNames();
                bool is_parameter = std::find(param_names.begin(), param_names.end(), var_name) != param_names.end();
                std::vector<double> values = is_parameter
                    ? ode_solution.GetParameterSeries(var_name, &ode_system)
                    : ode_solution.GetAnyVariable(var_name);
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

            // Exports results to csv (with time-resolved parameters)
            sth::ExportCsv("{{ ode_class_name }}.csv", ode_solution, ode_system, &ode_solution.rGetParametersPerStep());
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

#endif // {{ test_header_guard }}
