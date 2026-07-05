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

#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"

#include "SbmlTestHelpers.hpp"

void sbmltesthelpers::ExportCsv(const std::string& rFilename,
                                OdeSolution& rOdeSolution,
                                AbstractOdeSystem& rOdeSystem,
                                const std::vector<std::vector<double> >* pParamsPerStep)
{
    OutputFileHandler handler("");
    out_stream file = handler.OpenOutputFile(rFilename);

    // Times
    const std::vector<double>& time_data = rOdeSolution.rGetTimes();

    // State variables
    const std::vector<std::string>& svar_names = rOdeSystem.rGetStateVariableNames();
    const std::vector<std::vector<double> >& svar_data = rOdeSolution.rGetSolutions();

    // Derived quantities
    rOdeSolution.CalculateDerivedQuantitiesAndParameters(&rOdeSystem);
    const std::vector<std::string>& dq_names = rOdeSystem.rGetDerivedQuantityNames();
    const std::vector<std::vector<double> >& dq_data = rOdeSolution.rGetDerivedQuantities(&rOdeSystem);

    // Parameters
    const std::vector<std::string>& param_names = rOdeSystem.rGetParameterNames();
    const std::vector<double>& param_data = rOdeSolution.rGetParameters(&rOdeSystem);

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

    if ((param_data.empty() && !param_names.empty()) || (!param_data.empty() && param_data.size() != param_names.size()))
    {
        throw std::length_error("Number of parameter names do not match data.");
    }

    // Write column headings
    (*file) << "time";
    for (const auto& name : svar_names)
    {
        (*file) << "," << name;
    }
    for (const auto& name : dq_names)
    {
        (*file) << "," << name;
    }
    for (const auto& name : param_names)
    {
        (*file) << "," << name;
    }
    (*file) << '\n';

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
        if (pParamsPerStep != nullptr && i < pParamsPerStep->size())
        {
            // Time-resolved parameters (e.g. a parameter changed by an event).
            for (unsigned j = 0; j < (*pParamsPerStep)[i].size(); j++)
            {
                (*file) << "," << (*pParamsPerStep)[i][j];
            }
        }
        else if (!param_data.empty())
        {
            for (unsigned j = 0; j < param_data.size(); j++)
            {
                (*file) << "," << param_data[j];
            }
        }
        (*file) << '\n';
    }
    (*file) << std::flush;
    file->close();
}

std::string sbmltesthelpers::ToString(double value, unsigned precision)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}
