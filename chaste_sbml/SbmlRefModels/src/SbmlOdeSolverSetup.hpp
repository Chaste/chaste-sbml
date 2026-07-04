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

#ifndef SBML_ODE_SOLVER_SETUP_HPP_
#define SBML_ODE_SOLVER_SETUP_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCellCycleModelOdeSolver.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/**
 * Default ODE solver settings and setup shared by the SBML cell-cycle and SRN model base classes.
 */

namespace sbmlodesolversetup
{
/** Maximum number of internal steps CVODE may take per Solve() call. */
inline constexpr long int CVODE_MAX_STEPS = 10000;

/** CVODE relative tolerance. */
inline constexpr double CVODE_REL_TOL = 1e-6;

/** CVODE absolute tolerance. */
inline constexpr double CVODE_ABS_TOL = 1e-9;

/** Fixed timestep for the Runge-Kutta fallback solver used when CVODE is unavailable. */
inline constexpr double RK4_TIME_STEP = 0.0001;

/**
 * Give a model the default SBML ODE solver if it does not already have one: CVODE where available
 * (with stopping-event checking enabled), otherwise the fixed-step Runge-Kutta fallback. A solver
 * already supplied to the model is left untouched.
 *
 * @param pModel the model being constructed; used to key the CellCycleModelOdeSolver singleton
 *               (via the template parameter) and to set the fallback timestep
 * @param rpOdeSolver the model's solver pointer, set in place only when currently null
 */
template <typename Model>
void SetUpDefaultOdeSolver([[maybe_unused]] Model* pModel,
                           boost::shared_ptr<AbstractCellCycleModelOdeSolver>& rpOdeSolver)
{
    if (rpOdeSolver)
    {
        return; // A solver was supplied explicitly; leave it as-is.
    }

#ifdef CHASTE_CVODE
    // Default to CVODE where available
    rpOdeSolver = CellCycleModelOdeSolver<Model, CvodeAdaptor>::Instance();
    rpOdeSolver->Initialise();
    rpOdeSolver->SetMaxSteps(CVODE_MAX_STEPS);
    rpOdeSolver->SetTolerances(CVODE_REL_TOL, CVODE_ABS_TOL);
    // CVODE needs to be instructed to check for stopping events
    rpOdeSolver->CheckForStoppingEvents();
#else
    // Default to Chaste Runge-Kutta solver where CVODE is not available
    rpOdeSolver = CellCycleModelOdeSolver<Model, RungeKutta4IvpOdeSolver>::Instance();
    rpOdeSolver->Initialise();
    pModel->SetDt(RK4_TIME_STEP);
#endif // CHASTE_CVODE
}
} // namespace sbmlodesolversetup

#endif // SBML_ODE_SOLVER_SETUP_HPP_
