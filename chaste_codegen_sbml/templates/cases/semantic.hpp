#ifndef {{ test_header_guard }}
#define {{ test_header_guard }}

#include <cmath>
#include <string>
#include <vector>

#include <cxxtest/TestSuite.h>

#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
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
            OdeSolution ode_solution;

            CvodeAdaptor solver;
            solver.CheckForStoppingEvents();

            // Settings
            double start = {{ test_settings["start"] }};
            double duration = {{ test_settings["duration"] }};
            double steps = {{ test_settings["steps"] }};
            double tol_absolute = {{ test_settings["absolute"] }} * 10.0; // TODO: review tolerance values
            double tol_relative = {{ test_settings["relative"] }} * 10.0;

            double end = start + duration;
            double timestep = duration / steps;
            double sampling = timestep;

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

                OdeSolution next_solution = solver.Solve(&ode_system, initial_conditions, start, end, timestep, sampling);

                // Append new solution to existing solution
                sth::AppendOdeSolution(&ode_solution, &next_solution);
            }

            ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);

            // Expected results
            std::vector<std::string> expected_result_columns = { {{ test_result_columns }} };
            std::set<std::string> expected_amounts{ {{ test_amounts }} };
            std::set<std::string> expected_concentrations{ {{ test_concentrations }} };

            std::vector<std::vector<double> > expected_result_data = {
                {{ test_result_data }}
            };

            // Check sampling times
            TS_ASSERT_EQUALS(ode_solution.rGetTimes().size(), expected_result_data.size());
            for (unsigned i = 0; i < ode_solution.rGetTimes().size(); i++)
            {
                TS_ASSERT_DELTA(ode_solution.rGetTimes()[i], expected_result_data[i][0], 1e-6);
            }

            // Check variable values
            for (unsigned j = 1; j < expected_result_columns.size(); j++)
            {
                std::string var_name = expected_result_columns[j];
                TSM_ASSERT_EQUALS(var_name.c_str(), ode_system.HasAnyVariable(var_name), true);

                if (expected_amounts.find(var_name) != expected_amounts.end())
                {
                    // Use amount variable
                    if (ode_system.HasAnyVariable("amount__" + var_name))
                    {
                        var_name = "amount__" + var_name;
                    }
                }
                else
                {
                    TSM_ASSERT(var_name + " is neither amount nor concentration",
                               expected_concentrations.find(var_name) != expected_concentrations.end());
                }

                std::vector<double> values = ode_solution.GetAnyVariable(var_name);
                TS_ASSERT_EQUALS(values.size(), expected_result_data.size());
                for (unsigned i = 0; i < values.size(); i++)
                {
                    double delta = std::abs(expected_result_data[i][j] - values[i]);
                    double tol = tol_absolute + tol_relative * std::abs(expected_result_data[i][j]);
                    std::string msg(sth::ToString(values[i]) + " vs " + sth::ToString(expected_result_data[i][j])
                                    + " at " + sth::ToString(ode_solution.rGetTimes()[i], 3) + " for " + var_name);
                    TSM_ASSERT_LESS_THAN_EQUALS(msg.c_str(), delta, tol);
                }
            }

            // Exports results to csv
            // sth::ExportCsv("{{ ode_class_name }}.csv", ode_solution, ode_system);
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
