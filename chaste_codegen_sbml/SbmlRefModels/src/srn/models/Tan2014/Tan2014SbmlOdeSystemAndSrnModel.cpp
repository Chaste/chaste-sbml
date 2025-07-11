#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Tan2014SbmlOdeSystemAndSrnModel.hpp"

namespace sm = sbmlmath;

Tan2014SbmlOdeSystem::Tan2014SbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(6)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>);

    // STATE VARIABLES
    bcat_cm = 46.6;
    ligand_cm = 581.1;
    complex_cm = 418.9;
    bcat_nu = 32.6;
    ligand_nu = 516.8;
    complex_nu = 483.2;

    SetDefaultInitialCondition(0, bcat_cm);
    SetDefaultInitialCondition(1, ligand_cm);
    SetDefaultInitialCondition(2, complex_cm);
    SetDefaultInitialCondition(3, bcat_nu);
    SetDefaultInitialCondition(4, ligand_nu);
    SetDefaultInitialCondition(5, complex_nu);

    if (stateVariables.size() == 6)
    {
        bcat_cm = stateVariables[0];
        ligand_cm = stateVariables[1];
        complex_cm = stateVariables[2];
        bcat_nu = stateVariables[3];
        ligand_nu = stateVariables[4];
        complex_nu = stateVariables[5];
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("Tan2014SbmlOdeSystem: Expected 6 state variables, got " + std::to_string(stateVariables.size()));
    }

    mStateVariables.push_back(bcat_cm);
    mStateVariables.push_back(ligand_cm);
    mStateVariables.push_back(complex_cm);
    mStateVariables.push_back(bcat_nu);
    mStateVariables.push_back(ligand_nu);
    mStateVariables.push_back(complex_nu);

    // DERIVED QUANTITIES
    drag = 0.0;

    // VARIABLE PARAMETERS
    compartment = 1.0;
    CytosolMembrane = 1.16;
    nucleus = 0.65;
    wnt_level = 0.0;
    gamma = 1.0;

    mParameters.push_back(compartment);
    mParameters.push_back(CytosolMembrane);
    mParameters.push_back(nucleus);
    mParameters.push_back(wnt_level);
    mParameters.push_back(gamma);

    // RULE-BASED PARAMETERS

    // REACTIONS
    Bsynthesis = 0.0;
    kDegradation = 0.0;
    kC = 0.0;
    kN = 0.0;
    kdiffusion = 0.0;
    K_c_active = 0.0;
    K_n_active = 0.0;


    ProcessRules(0.0, mStateVariables);
}

Tan2014SbmlOdeSystem::~Tan2014SbmlOdeSystem()
{
}

void Tan2014SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    ProcessRules(time, rY);

    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // d[bcat_cm]/dt
    rDY[1] = (-kC) / CytosolMembrane; // d[ligand_cm]/dt
    rDY[2] = (kC) / CytosolMembrane; // d[complex_cm]/dt
    rDY[3] = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus; // d[bcat_nu]/dt
    rDY[4] = (-kN) / nucleus; // d[ligand_nu]/dt
    rDY[5] = (kN) / nucleus; // d[complex_nu]/dt

    // Scale time appropriately
}

std::vector<double> Tan2014SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double> &rY)
{
    ProcessRules(time, rY);

    std::vector<double> dqs;
    dqs.push_back(drag);
    return dqs;
}

void Tan2014SbmlOdeSystem::ProcessRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
    bcat_cm = rY[0];
    ligand_cm = rY[1];
    complex_cm = rY[2];
    bcat_nu = rY[3];
    ligand_nu = rY[4];
    complex_nu = rY[5];

    // VARIABLE PARAMETERS
    compartment = GetParameter(0);
    CytosolMembrane = GetParameter(1);
    nucleus = GetParameter(2);
    wnt_level = GetParameter(3);
    gamma = GetParameter(4);

    // RULES
    drag = sm::piecewise((complex_cm - 700.0) / 10.0, sm::gt((complex_cm - 700.0) / 10.0, 1.0), 1.0);

    // REACTIONS
    Bsynthesis = Bsyn * CytosolMembrane;

    kDegradation = CytosolMembrane * kdegradation * gamma * bcat_cm * (1.0 - 0.5 * wnt_level);

    kC = CytosolMembrane * (kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm);

    kN = nucleus * (kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu);

    kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);

    K_c_active = K_c_active_k * bcat_cm;

    K_n_active = K_n_active_k * bcat_nu;

}




// FUNCTIONS

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
    this->mDerivedQuantityNames.push_back("drag");
    this->mDerivedQuantityUnits.push_back("non-dim");


    // PARAMETERS
    this->mParameterNames.push_back("compartment");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("CytosolMembrane");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("nucleus");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014SbmlOdeSystem, 6> Tan2014SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014SbmlOdeSystem, 6)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SbmlSrnModel)