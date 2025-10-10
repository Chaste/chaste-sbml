
#include <fstream>
#include <stdexcept>

#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"

#include "SbmlTestHelpers.hpp"

void sbmltesthelpers::export_csv(const std::string& filename,
                                 OdeSolution& ode_solution,
                                 AbstractOdeSystem& ode_system)
{
    OutputFileHandler handler("");
    out_stream file = handler.OpenOutputFile(filename);

    // Times
    const std::vector<double>& time_data = ode_solution.rGetTimes();

    // State variables
    const std::vector<std::string>& svar_names = ode_system.rGetStateVariableNames();
    const std::vector<std::vector<double> >& svar_data = ode_solution.rGetSolutions();

    // Derived quantities
    ode_solution.CalculateDerivedQuantitiesAndParameters(&ode_system);
    const std::vector<std::string>& dq_names = ode_system.rGetDerivedQuantityNames();
    const std::vector<std::vector<double> >& dq_data = ode_solution.rGetDerivedQuantities(&ode_system);

    // Sanity checks
    if (time_data.empty())
    {
        throw std::invalid_argument("OdeSolution contains no time points.");
    }

    if (svar_data.empty() && dq_data.empty())
    {
        throw std::invalid_argument("OdeSolution contains no state variables or derived quantities.");
    }

    if (!svar_data.empty() && (svar_data.size() != time_data.size()))
    {
        throw std::length_error("Number of state variable data rows do not match time points.");
    }

    if (!dq_data.empty() && (dq_data.size() != time_data.size()))
    {
        throw std::length_error("Number of derived quantity data rows do not match time points.");
    }

    if ((svar_data.empty() && !svar_names.empty()) || (!svar_data.empty() && svar_data[0].size() != svar_names.size()))
    {
        throw std::length_error("Number of state variable names do not match data.");
    }

    if ((dq_data.empty() && !dq_names.empty()) || (!dq_data.empty() && dq_data[0].size() != dq_names.size()))
    {
        throw std::length_error("Number of derived quantity names do not match data.");
    }

    // Write column headings
    (*file) << "time";
    if (!svar_names.empty())
    {
        for (unsigned i = 0; i < svar_names.size(); i++)
        {
            (*file) << "," << svar_names[i];
        }
    }
    if (!dq_names.empty())
    {
        for (unsigned i = 0; i < dq_names.size(); i++)
        {
            (*file) << "," << dq_names[i];
        }
    }
    (*file) << std::endl
            << std::flush;

    // Write data
    for (unsigned i = 0; i < time_data.size(); i++)
    {
        (*file) << time_data[i];
        if (!svar_data.empty())
        {
            for (unsigned j = 0; j < svar_data[i].size(); j++)
            {
                (*file) << "," << svar_data[i][j];
            }
        }
        if (!dq_data.empty())
        {
            for (unsigned j = 0; j < dq_data[i].size(); j++)
            {
                (*file) << "," << dq_data[i][j];
            }
        }
        (*file) << std::endl
                << std::flush;
    }
    file->close();
}
