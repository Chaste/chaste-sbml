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

#include "SbmlTestOdeSolution.hpp"

void SbmlTestOdeSolution::AppendSegment(OdeSolution& rSegment, AbstractOdeSystem* pSystem)
{
    // Snapshot the system's current parameters (constant within this segment).
    std::vector<double> params(pSystem->GetNumberOfParameters());
    for (unsigned i = 0; i < params.size(); ++i)
    {
        params[i] = pSystem->GetParameter(i);
    }

    std::vector<double>& r_times = rGetTimes();
    std::vector<std::vector<double> >& r_solutions = rGetSolutions();

    if (r_times.empty())
    {
        // First segment: seed the solution and record one parameter set per step.
        SetOdeSystemInformation(pSystem->GetSystemInformation());
        r_times = rSegment.rGetTimes();
        r_solutions = rSegment.rGetSolutions();
        mParametersPerStep.assign(r_times.size(), params);
    }
    else
    {
        // Later segment: drop the duplicated restart point shared with the previous
        // segment's end, then append the rest.
        r_times.pop_back();
        r_solutions.pop_back();
        mParametersPerStep.pop_back();

        const std::vector<double>& r_seg_times = rSegment.rGetTimes();
        const std::vector<std::vector<double> >& r_seg_solutions = rSegment.rGetSolutions();

        r_times.insert(r_times.end(), r_seg_times.begin() + 1, r_seg_times.end());
        r_solutions.insert(r_solutions.end(), r_seg_solutions.begin() + 1, r_seg_solutions.end());
        mParametersPerStep.insert(mParametersPerStep.end(), r_seg_times.size() - 1, params);
    }

    SetNumberOfTimeSteps(r_times.size());
}

std::vector<double> SbmlTestOdeSolution::GetParameterSeries(const std::string& rName, AbstractOdeSystem* pSystem) const
{
    unsigned index = pSystem->GetParameterIndex(rName);
    std::vector<double> series(mParametersPerStep.size());
    for (unsigned i = 0; i < series.size(); ++i)
    {
        series[i] = mParametersPerStep[i][index];
    }
    return series;
}
