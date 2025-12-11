#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Tan2014SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Tan2014SbmlOdeSystem::Tan2014SbmlOdeSystem()
        : AbstractSbmlOdeSystem(6, 2, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>);

    // VARIABLE PARAMETERS
    wnt_level = 0.0;
    gamma = 1.0;

    // STATE VARIABLES
    bcat_cm = 46.6;
    ligand_cm = 581.1;
    complex_cm = 418.9;
    bcat_nu = 32.6;
    ligand_nu = 516.8;
    complex_nu = 483.2;

    // DERIVED QUANTITIES
    compartment = 1.0;
    CytosolMembrane = 1.16;
    nucleus = 0.65;
    drag = 1.0;

    // INITIAL ASSIGNMENTS
    RunInitialAssignments(0.0);

    // ODE SYSTEM INFORMATION
    SetDefaultInitialCondition(0, bcat_cm);
    SetDefaultInitialCondition(1, ligand_cm);
    SetDefaultInitialCondition(2, complex_cm);
    SetDefaultInitialCondition(3, bcat_nu);
    SetDefaultInitialCondition(4, ligand_nu);
    SetDefaultInitialCondition(5, complex_nu);

    mStateVariables.push_back(bcat_cm);
    mStateVariables.push_back(ligand_cm);
    mStateVariables.push_back(complex_cm);
    mStateVariables.push_back(bcat_nu);
    mStateVariables.push_back(ligand_nu);
    mStateVariables.push_back(complex_nu);

    mParameters.push_back(wnt_level);
    mParameters.push_back(gamma);

    // REACTIONS
    RunReactions(0.0);

    // EVENTS
}

Tan2014SbmlOdeSystem::~Tan2014SbmlOdeSystem()
{
}

std::vector<double> Tan2014SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;

    RunModelRules(time, rY);

    dqs.push_back(compartment);
    dqs.push_back(CytosolMembrane);
    dqs.push_back(nucleus);
    dqs.push_back(drag);

    // AMOUNTS
    double amt__bcat_cm = bcat_cm * CytosolMembrane;
    double amt__ligand_cm = ligand_cm * CytosolMembrane;
    double amt__complex_cm = complex_cm * CytosolMembrane;
    double amt__bcat_nu = bcat_nu * nucleus;
    double amt__ligand_nu = ligand_nu * nucleus;
    double amt__complex_nu = complex_nu * nucleus;
    double amt__drag = drag * CytosolMembrane;

    dqs.push_back(amt__bcat_cm);
    dqs.push_back(amt__ligand_cm);
    dqs.push_back(amt__complex_cm);
    dqs.push_back(amt__bcat_nu);
    dqs.push_back(amt__ligand_nu);
    dqs.push_back(amt__complex_nu);
    dqs.push_back(amt__drag);
    return dqs;
}

void Tan2014SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // d[bcat_cm]/dt
    rDY[1] = -kC / CytosolMembrane;                                                                     // d[ligand_cm]/dt
    rDY[2] = kC / CytosolMembrane;                                                                      // d[complex_cm]/dt
    rDY[3] = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus;                                    // d[bcat_nu]/dt
    rDY[4] = -kN / nucleus;                                                                             // d[ligand_nu]/dt
    rDY[5] = kN / nucleus;                                                                              // d[complex_nu]/dt

    // TODO: Scale time appropriately
}

double Tan2014SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void Tan2014SbmlOdeSystem::RunAssignmentRules(double time)
{
    drag = sm::piecewise((complex_cm - 700.0) / 10.0, sm::gt((complex_cm - 700.0) / 10.0, 1.0), 1.0);
}

// INITIAL ASSIGNMENTS
void Tan2014SbmlOdeSystem::RunInitialAssignments(double time)
{
    drag = sm::piecewise((complex_cm - 700.0) / 10.0, sm::gt((complex_cm - 700.0) / 10.0, 1.0), 1.0); //
}

// REACTIONS
void Tan2014SbmlOdeSystem::RunReactions(double time)
{
    Bsynthesis = Bsyn * CytosolMembrane;

    kDegradation = CytosolMembrane * kdegradation * gamma * bcat_cm * (1.0 - 0.5 * wnt_level);

    kC = CytosolMembrane * (kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm);

    kN = nucleus * (kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu);

    kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);

    K_c_active = K_c_active_k * bcat_cm;

    K_n_active = K_n_active_k * bcat_nu;
}

// VARIABLE PARAMETERS
void Tan2014SbmlOdeSystem::UpdateParameters(double time)
{
    wnt_level = GetParameter(0);
    gamma = GetParameter(1);
}

// STATE VARIABLES
void Tan2014SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    bcat_cm = rStateVariables[0];
    ligand_cm = rStateVariables[1];
    complex_cm = rStateVariables[2];
    bcat_nu = rStateVariables[3];
    ligand_nu = rStateVariables[4];
    complex_nu = rStateVariables[5];
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("bcat_cm");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(46.6);

    this->mVariableNames.push_back("ligand_cm");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(581.1);

    this->mVariableNames.push_back("complex_cm");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(418.9);

    this->mVariableNames.push_back("bcat_nu");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(32.6);

    this->mVariableNames.push_back("ligand_nu");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(516.8);

    this->mVariableNames.push_back("complex_nu");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(483.2);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("compartment");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CytosolMembrane");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("nucleus");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__bcat_cm");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__ligand_cm");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__complex_cm");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__bcat_nu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__ligand_nu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__complex_nu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)