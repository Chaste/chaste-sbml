#ifndef TEST_SEMANTIC_00928_L3_V2_SBML_HPP_
#define TEST_SEMANTIC_00928_L3_V2_SBML_HPP_

#include <cmath>
#include <string>
#include <vector>

#include <cxxtest/TestSuite.h>

#include "CvodeAdaptor.hpp"
#include "OdeSolution.hpp"
#include "SbmlTestHelpers.hpp"

#include "Semantic00928L3V2SbmlOdeSystem.hpp"

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
            Semantic00928L3V2SbmlOdeSystem ode_system;
            OdeSolution ode_solution;

            CvodeAdaptor solver;
            solver.CheckForStoppingEvents();

            // Settings
            double start = 0;
            double duration = 5;
            double steps = 50;

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
            std::vector<std::string> expected_result_columns = { "time", "S1", "S2" };
            std::set<std::string> expected_amounts{ "S1", "S2" };

            std::vector<std::vector<double> > expected_result_data = {
                { 0, 0.00015, 0 },
                { 0.1, 0.0001357256127053939, 1.427438729460607e-005 },
                { 0.2, 0.0001228096129616973, 2.719038703830272e-005 },
                { 0.3, 0.0001111227331022577, 3.887726689774233e-005 },
                { 0.4, 0.0001005480069053459, 4.945199309465411e-005 },
                { 0.5, 9.097959895689501e-005, 5.902040104310499e-005 },
                { 0.6, 8.232174541410396e-005, 6.767825458589604e-005 },
                { 0.7, 7.448779556871142e-005, 7.551220443128858e-005 },
                { 0.8, 6.739934461758323e-005, 8.260065538241677e-005 },
                { 0.9, 6.098544896108986e-005, 8.901455103891014e-005 },
                { 1, 5.518191617571635e-005, 9.481808382428365e-005 },
                { 1.1, 4.993066255471193e-005, 0.0001000693374452881 },
                { 1.2, 4.517913178683031e-005, 0.0001048208682131697 },
                { 1.3, 4.087976895510189e-005, 0.0001091202310448981 },
                { 1.4, 3.698954459124097e-005, 0.000113010455408759 },
                { 1.5, 3.346952402226447e-005, 0.0001165304759777355 },
                { 1.6, 3.028447769919831e-005, 0.0001197155223008017 },
                { 1.7, 2.740252860791019e-005, 0.0001225974713920898 },
                { 1.8, 2.479483323323798e-005, 0.000125205166766762 },
                { 1.9, 2.243529288339526e-005, 0.0001275647071166047 },
                { 2, 2.03002924854919e-005, 0.0001296997075145081 },
                { 2.1, 1.836846423794729e-005, 0.0001316315357620527 },
                { 2.2, 1.662047375435008e-005, 0.0001333795262456499 },
                { 2.3, 1.503882655842056e-005, 0.0001349611734415794 },
                { 2.4, 1.360769299341187e-005, 0.0001363923070065881 },
                { 2.5, 1.231274979358482e-005, 0.0001376872502064152 },
                { 2.6, 1.114103673215008e-005, 0.0001388589632678499 },
                { 2.7, 1.008082691096246e-005, 0.0001399191730890375 },
                { 2.8, 9.121509393782692e-006, 0.0001408784906062173 },
                { 2.9, 8.253483008461082e-006, 0.0001417465169915389 },
                { 3, 7.468060255179591e-006, 0.0001425319397448204 },
                { 3.1, 6.757380359033671e-006, 0.0001432426196409663 },
                { 3.2, 6.114330596754931e-006, 0.0001438856694032451 },
                { 3.3, 5.532475110186e-006, 0.000144467524889814 },
                { 3.4, 5.00599049404891e-006, 0.0001449940095059511 },
                { 3.5, 4.529607513347775e-006, 0.0001454703924866522 },
                { 3.6, 4.098558367093884e-006, 0.0001459014416329061 },
                { 3.7, 3.708528970550908e-006, 0.0001462914710294491 },
                { 3.8, 3.355615778424839e-006, 0.0001466443842215752 },
                { 3.9, 3.036286716870657e-006, 0.0001469637132831294 },
                { 4, 2.747345833310127e-006, 0.0001472526541666899 },
                { 4.1, 2.485901310264186e-006, 0.0001475140986897358 },
                { 4.2, 2.249336523071655e-006, 0.0001477506634769283 },
                { 4.3, 2.03528385183014e-006, 0.0001479647161481699 },
                { 4.4, 1.841600985460266e-006, 0.0001481583990145397 },
                { 4.5, 1.666349480736346e-006, 0.0001483336505192637 },
                { 4.6, 1.507775361695036e-006, 0.000148492224638305 },
                { 4.7, 1.364291565254372e-006, 0.0001486357084347456 },
                { 4.8, 1.234462057353003e-006, 0.000148765537942647 },
                { 4.9, 1.116987460638651e-006, 0.0001488830125393614 },
                { 5, 1.01069204986282e-006, 0.0001489893079501372 }
            };

            // Check variable values
            double tol_absolute = 1.000000e-007 * 10.0; // TODO: review tolerance values
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
            sth::ExportCsv("Semantic00928L3V2SbmlOdeSystem.csv", ode_solution, ode_system);
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

#endif // TEST_SEMANTIC_00928_L3_V2_SBML_HPP_