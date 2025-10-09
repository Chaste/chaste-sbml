#ifndef {{ test_header_guard }}
#define {{ test_header_guard }}

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
            CvodeAdaptor solver;
            OdeSolution ode_solution;

            // Settings
            double start = {{ test_settings["start"] }};
            double duration = {{ test_settings["duration"] }};
            double steps = {{ test_settings["steps"] }};

            double end = start + duration;
            double timestep = duration / steps;
            double sampling = timestep;

            // Solve
            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            ode_solution = solver.Solve(&ode_system, initial_conditions, start, end, timestep, sampling);
            ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);

            // Expected results
            std::vector<std::string> expected_result_columns = { {{ test_result_columns }} };

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
                const std::string& var_name = expected_result_columns[j];
                TSM_ASSERT_EQUALS(var_name.c_str(), ode_system.HasAnyVariable(var_name), true);

                std::vector<double> values = ode_solution.GetAnyVariable(var_name);
                TS_ASSERT_EQUALS(values.size(), expected_result_data.size());
                for (unsigned i = 0; i < values.size(); i++)
                {
                    TSM_ASSERT_DELTA(var_name.c_str(), values[i], expected_result_data[i][j], 1e-4);
                }
            }

            // Exports results to csv
            // sth::export_csv("{{ ode_class_name }}.csv", ode_system.rGetStateVariableNames(), ode_solution);
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
