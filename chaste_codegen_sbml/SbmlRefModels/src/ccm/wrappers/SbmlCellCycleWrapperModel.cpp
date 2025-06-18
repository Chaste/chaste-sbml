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
#include "AbstractOdeBasedCellCycleModel.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
//#include "Exception.hpp"
// #include "Debug.hpp"

#include "SbmlCellCycleWrapperModel.hpp"


template<typename SBMLODE, unsigned SIZE>
SbmlCellCycleWrapperModel<SBMLODE, SIZE>::SbmlCellCycleWrapperModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedCellCycleModel(SIZE, pOdeSolver)
{
    if (!mpOdeSolver)
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<SbmlCellCycleWrapperModel<SBMLODE, SIZE>, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
        mpOdeSolver->SetTolerances(1e-6, 1e-8);
#else
        mpOdeSolver = CellCycleModelOdeSolver<SbmlCellCycleWrapperModel<SBMLODE, SIZE>, BackwardEulerIvpOdeSolver>::Instance();
        mpOdeSolver->SetSizeOfOdeSystem(SIZE);
        mpOdeSolver->Initialise();
        SetDt(0.1/90.0);
#endif //CHASTE_CVODE
    }
}

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

// template<typename SBMLODE, unsigned SIZE>
// void SbmlCellCycleWrapperModel<SBMLODE, SIZE>::Initialise()
// {
//     assert(mpOdeSystem == nullptr);
//     mpOdeSystem = new SBMLODE;
//     mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

//     AbstractOdeBasedCellCycleModel::Initialise();
// }

template<typename SBMLODE, unsigned SIZE>
void SbmlCellCycleWrapperModel<SBMLODE, SIZE>::Initialise()
{    
    assert(mpOdeSystem == nullptr);
    mpOdeSystem = new SBMLODE;

    AbstractOdeBasedCellCycleModel::Initialise();

	//Initialise cell data

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
void SbmlCellCycleWrapperModel<SBMLODE, SIZE>::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem != nullptr);

    /**
     * This model needs the protein concentrations and phase resetting to G0/G1.
     *
     * In theory, the solution to the Tyson-Novak equations should exhibit stable
     * oscillations, and we only need to halve the mass of the cell each period.
     *
     * However, the backward Euler solver used to solve the equations
     * currently returns a solution that diverges after long times, so
     * we must reset the initial conditions each period.
     *
     * When running with CVODE however we can use the halving the mass of the cell method.
     */
// #ifdef CHASTE_CVODE
//     mpOdeSystem->rGetStateVariables()[5] = 0.5*mpOdeSystem->rGetStateVariables()[5];
// #else
//     mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
// #endif //CHASTE_CVODE
}

template<typename SBMLODE, unsigned SIZE>
void SbmlCellCycleWrapperModel<SBMLODE, SIZE>::InitialiseDaughterCell()
{
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
}

template<typename SBMLODE, unsigned SIZE>
AbstractCellCycleModel* SbmlCellCycleWrapperModel<SBMLODE, SIZE>::CreateCellCycleModel()
{
	return new SbmlCellCycleWrapperModel(*this);
}


// AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
// {
//     return new TysonNovakCellCycleModel(*this);
// }

template<typename SBMLODE, unsigned SIZE>
double SbmlCellCycleWrapperModel<SBMLODE, SIZE>::GetAverageTransitCellCycleTime()
{
    return 1.25;
}

template<typename SBMLODE, unsigned SIZE>
double SbmlCellCycleWrapperModel<SBMLODE, SIZE>::GetAverageStemCellCycleTime()
{
    return 1.25;
}

template<typename SBMLODE, unsigned SIZE>
bool SbmlCellCycleWrapperModel<SBMLODE, SIZE>::CanCellTerminallyDifferentiate()
{
    return false;
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
