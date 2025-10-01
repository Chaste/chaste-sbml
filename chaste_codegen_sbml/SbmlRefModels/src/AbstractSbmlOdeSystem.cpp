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

#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

#include "AbstractSbmlOdeSystem.hpp"

AbstractSbmlOdeSystem::AbstractSbmlOdeSystem(unsigned numberOfStateVariables, unsigned numberOfParameters, unsigned numberOfEvents)
        : AbstractOdeSystem(numberOfStateVariables),
          mNumberOfParameters(numberOfParameters),
          mNumberOfEvents(numberOfEvents)
{
    if (mNumberOfEvents > 0)
    {
        mEventType.resize(mNumberOfEvents, SbmlEventType::UNKNOWN);

        mEventSatisfied.resize(mNumberOfEvents, true); // Prevent events from triggering at the start
        mEventTriggered.resize(mNumberOfEvents, false);

        if (mNumberOfStateVariables > 0)
        {
            mEventAdjustedStateVars.resize(mNumberOfStateVariables, false);
            mEventAdjustedStateValues.resize(mNumberOfStateVariables, 0.0);
        }

        if (mNumberOfParameters > 0)
        {
            mEventAdjustedParameters.resize(mNumberOfParameters, false);
            mEventAdjustedParameterValues.resize(mNumberOfParameters, 0.0);
        }
    }
}

AbstractSbmlOdeSystem::~AbstractSbmlOdeSystem()
{
}

void AbstractSbmlOdeSystem::AdjustParameters(double time)
{
    for (unsigned i = 0; i < mEventAdjustedParameters.size(); ++i)
    {
        if (mEventAdjustedParameters[i])
        {
            SetParameter(i, mEventAdjustedParameterValues[i]);
        }
    }

    for (unsigned i = 0; i < mEventAdjustedStateVars.size(); ++i)
    {
        if (mEventAdjustedStateVars[i])
        {
            SetStateVariable(i, mEventAdjustedStateValues[i]);
            mEventAdjustedStateVars[i] = false;
        }
    }
}

double AbstractSbmlOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    return ProcessModelEvents(time, rY);
}

bool AbstractSbmlOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    return ProcessModelEvents(time, rY) == 0.0;
}

bool AbstractSbmlOdeSystem::HasEventOccurred(SbmlEventType eventType)
{
    for (unsigned i = 0; i < mEventTriggered.size(); ++i)
    {
        if (mEventTriggered[i] && mEventType[i] == eventType)
        {
            return true;
        }
    }
    return false;
}

void AbstractSbmlOdeSystem::ResetEventsOccurred()
{
    std::fill(mEventTriggered.begin(), mEventTriggered.end(), false);
}
