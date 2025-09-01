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

#include "AbstractOdeSrnModel.hpp"
#include "AbstractSbmlOdeSystem.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

#include "AbstractSbmlSrnModel.hpp"

AbstractSbmlSrnModel::AbstractSbmlSrnModel(unsigned stateSize, boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractOdeSrnModel(stateSize, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        // Default to CVODE where available
        mpOdeSolver = CellCycleModelOdeSolver<AbstractSbmlSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000); // Safe default
        // CVODE needs to be instructed to check for stopping events
        mpOdeSolver->CheckForStoppingEvents();
#else
        // Default to Chaste Runge-Kutta solver where CVODE is not available
        mpOdeSolver = CellCycleModelOdeSolver<AbstractSbmlSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        this->SetDt(0.0001); // Safe default
#endif // CHASTE_CVODE
    }

    assert(mpOdeSolver->IsSetUp());
}

double AbstractSbmlSrnModel::GetStateVariable(const std::string& rName)
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetStateVariable(rName);
}

void AbstractSbmlSrnModel::Initialise(AbstractSbmlOdeSystem* pOdeSystem)
{
    AbstractOdeSrnModel::Initialise(pOdeSystem);

    // Initialise cell data
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    /* Custom behaviour: store the state variables as cell data and set any parameters
     * using cell data, so that we can visualise different concentrations in Paraview.
     */

    std::vector<std::string> stateVariableNames = mpOdeSystem->rGetStateVariableNames();

    for (unsigned i = 0; i < stateVariableNames.size(); i++)
    {
        std::string stateName = stateVariableNames[i];
        double stateValue = mpOdeSystem->rGetStateVariables()[i];

        // Set current state variable value as cell data
        mpCell->GetCellData()->SetItem(stateName, stateValue);
    }
}

void AbstractSbmlSrnModel::SimulateToCurrentTime()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    /* Custom behaviour: store the state variables and state parameters as
     * cell data so that we can visualise different concentrations in Paraview.
     */

    // Set state parameters as cell data
    std::vector<std::string> parameterNames = mpOdeSystem->rGetParameterNames();
    for (unsigned i = 0; i < parameterNames.size(); i++)
    {
        std::string parameterName = parameterNames[i];
        double parameterValue = mpOdeSystem->GetParameter(i);
        mpCell->GetCellData()->SetItem(parameterName, parameterValue);
    }

    // Set state variables as cell data
    std::vector<std::string> stateVariableNames = mpOdeSystem->rGetStateVariableNames();

    for (unsigned i = 0; i < stateVariableNames.size(); i++)
    {
        std::string stateName = stateVariableNames[i];
        double stateValue = mpOdeSystem->GetStateVariable(i);
        mpCell->GetCellData()->SetItem(stateName, stateValue);
    }

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}
