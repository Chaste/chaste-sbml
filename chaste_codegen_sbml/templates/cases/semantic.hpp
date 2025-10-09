#ifndef {{ test_header_guard }}
#define {{ test_header_guard }}

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
            OdeSolution solutions;

            double start = {{ test_settings["start"] }};
            double duration = {{ test_settings["duration"] }};
            double steps = {{ test_settings["steps"] }};

            double end = start + duration;
            double timestep = duration / steps;
            double sampling = timestep;

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            solutions = solver.Solve(&ode_system, initial_conditions, start, end, timestep, sampling);

            // Check results against expected values
            std::vector<std::vector<double> > expected_results = {
                {{ test_results }}
            };

            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i++)
            {
                TS_ASSERT_DELTA(solutions.rGetTimes()[i], expected_results[i][0], 1e-4);
                TS_ASSERT_DELTA(solutions.rGetSolutions()[i][0], expected_results[i][1], 1e-4);
                TS_ASSERT_DELTA(solutions.rGetSolutions()[i][1], expected_results[i][2], 1e-4);
            }

            // Exports results to csv
            // sth::export_csv("{{ ode_class_name }}.csv", ode_system.rGetStateVariableNames(), solutions);
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
