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

            double sampling = duration / steps;
            double timestep = sampling / 10.0;

            // Solve to each sample-grid point in turn. An event stops the integration early and
            // is applied, but only grid points are recorded - so the solution stays aligned with
            // the expected grid even when an event fires between grid points. Consecutive Solve
            // calls that start where the previous one ended reuse CVODE's state (no reset), so
            // continuity is preserved within a segment; CVODE is reset only after an event.
            // RecordPoint also stores the parameters at each grid point, so event-modified
            // parameters are time-resolved.
            const double grid_tol = sampling * 1e-6;
            const unsigned num_steps = static_cast<unsigned>(steps + 0.5);

            std::vector<double> y = ode_system.GetInitialConditions();
            ode_solution.RecordPoint(start, y, &ode_system);

            double current = start;
            for (unsigned k = 1; k <= num_steps; k++)
            {
                double target = start + k * sampling;
                bool event_at_target = false;

                while (current < target - grid_tol)
                {
                    OdeSolution segment = solver.Solve(&ode_system, y, current, target, timestep, target - current);
                    solver.SetForceReset(false); // chained Solves continue without re-initialising
                    y = segment.rGetSolutions().back();
                    double reached = segment.rGetTimes().back();

                    if (solver.StoppingEventOccurred() && reached < target - grid_tol)
                    {
                        // Event before the grid point: apply it and keep integrating to the grid point.
                        ode_system.SetStateVariables(y);
                        ode_system.AdjustParameters(reached);
                        y = ode_system.GetStateVariables();
                        solver.SetForceReset(true); // event changed the state/parameters: reset CVODE
                        current = reached;
                    }
                    else
                    {
                        event_at_target = solver.StoppingEventOccurred();
                        current = target;
                    }
                }

                ode_solution.RecordPoint(target, y, &ode_system);

                if (event_at_target)
                {
                    // Event exactly on the grid point: the pre-event value is recorded above, then apply it.
                    ode_system.SetStateVariables(y);
                    ode_system.AdjustParameters(target);
                    y = ode_system.GetStateVariables();
                    solver.SetForceReset(true);
                }
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
