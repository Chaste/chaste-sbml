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

#ifndef SBMLCELLCYCLEWRAPPERMODEL_CPP_
#define SBMLCELLCYCLEWRAPPERMODEL_CPP_

//#include "UblasIncludes.hpp"
#include "SbmlCellCycleWrapperModel.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AbstractOdeBasedCellCycleModel.hpp"
#include "CvodeAdaptor.hpp"
//#include "Exception.hpp"
#include "Debug.hpp"


template<typename SBMLODE, unsigned SIZE>
SbmlCellCycleWrapperModel<SBMLODE, SIZE>::SbmlCellCycleWrapperModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedCellCycleModel(SIZE, pOdeSolver)
{

	// TODO add CVODE solver. See VanLeeuwen CCM in trunk

    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {

    	#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<SbmlCellCycleWrapperModel<SBMLODE, SIZE>, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
		#else
        // TODO: Set solver type and dt from CLI
        mpOdeSolver = CellCycleModelOdeSolver<SbmlCellCycleWrapperModel<SBMLODE, SIZE>, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.0001); // This is small enough for both examples to converge
		#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

//New method for copy constructor
template<typename SBMLODE, unsigned SIZE>
SbmlCellCycleWrapperModel<SBMLODE, SIZE>::SbmlCellCycleWrapperModel(const SbmlCellCycleWrapperModel& rModel)
    : AbstractOdeBasedCellCycleModel(rModel)
{
    /*
     * Set each member variable of the new Cell Cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new Cell Cycle model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new Cell Cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new Cell Cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    assert(rModel.GetOdeSystem());
    SetOdeSystem(new SBMLODE(rModel.GetOdeSystem()->rGetStateVariables()));
}

template<typename SBMLODE, unsigned SIZE>
AbstractCellCycleModel* SbmlCellCycleWrapperModel<SBMLODE, SIZE>::CreateCellCycleModel()
{
	return new SbmlCellCycleWrapperModel(*this);
}

template<typename SBMLODE, unsigned SIZE>
void SbmlCellCycleWrapperModel<SBMLODE, SIZE>::Initialise()
{
	AbstractOdeBasedCellCycleModel::Initialise(new SBMLODE);

	//Initialise cell data
	assert(mpOdeSystem != NULL);
	assert(mpCell != NULL);

    /* Store the state variables and state parameters as
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

}


template<typename SBMLODE, unsigned SIZE>
void SbmlCellCycleWrapperModel<SBMLODE, SIZE>::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output.

    // Call direct parent class
    AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

template<typename SBMLODE, unsigned SIZE>
double SbmlCellCycleWrapperModel<SBMLODE, SIZE>::GetStateVariable(const std::string& rName)
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetStateVariable(rName);
}

#endif /* SBMLCELLCYCLEWRAPPERMODEL_CPP_ */
