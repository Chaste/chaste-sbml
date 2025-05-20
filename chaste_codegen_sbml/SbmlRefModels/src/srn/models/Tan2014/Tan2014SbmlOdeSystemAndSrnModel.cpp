#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Tan2014SbmlOdeSystemAndSrnModel.hpp"

Tan2014SbmlOdeSystem::Tan2014SbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(6)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 46.6); // bcat_cm
    SetDefaultInitialCondition(1, 581.1); // ligand_cm
    SetDefaultInitialCondition(2, 418.9); // complex_cm
    SetDefaultInitialCondition(3, 32.6); // bcat_nu
    SetDefaultInitialCondition(4, 516.8); // ligand_nu
    SetDefaultInitialCondition(5, 483.2); // complex_nu

    mParameters.push_back(0.0); // drag

    mParameters.push_back(0.0); // wnt_level
    mParameters.push_back(1.0); // gamma
    mParameters.push_back(1.0); // ComplexTransitThreshold

    if (stateVariables.size() > 0)
    {
        SetStateVariables(stateVariables);
    }
}

Tan2014SbmlOdeSystem::~Tan2014SbmlOdeSystem()
{
}


void Tan2014SbmlOdeSystem::Init()
{
    // COMPARTMENTS:
    compartment = 1.0;
    CytosolMembrane = 1.16;
    nucleus = 0.65;

    // PARAMETERS:
    wnt_level = 0.0;
    kdegradation = 0.0163;
    Bsyn = 1.306;
    K_n_active_k = 17.16;
    kC_k1 = 1e-05;
    kC_k2 = 0.000647;
    kN_k1 = 0.0001;
    kN_k2 = 0.00349;
    kdiffusion_k = 39.13;
    K_c_active_k = 4.5;
    gamma = 1.0;
    ComplexTransitThreshold = 1.0;

    // RULES:

}

void Tan2014SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    // SPECIES STATE VARIABLES:
    double bcat_cm = rY[0];
    double ligand_cm = rY[1];
    double complex_cm = rY[2];
    double bcat_nu = rY[3];
    double ligand_nu = rY[4];
    double complex_nu = rY[5];

    // SPECIES STATE PARAMETERS:
    double drag = GetParameter("drag");

    // RULES:
    drag = sbmlmath::sm_piecewise((complex_cm - 700) / 10, sbmlmath::sm_gt((complex_cm - 700) / 10, 1), 1);

    // UPDATE STATE PARAMETERS:
    SetParameter("drag", drag);

    // REACTIONS:

    double Bsynthesis = Bsyn * CytosolMembrane;

    double kDegradation = CytosolMembrane * kdegradation * gamma * bcat_cm * (1 - 0.5 * wnt_level);

    double kC = CytosolMembrane * (kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm);

    double kN = nucleus * (kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu);

    double kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);

    double K_c_active = K_c_active_k * bcat_cm;

    double K_n_active = K_n_active_k * bcat_nu;

    // ODES:
    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // dbcat_cm/dt
    rDY[1] = (-kC) / CytosolMembrane; // dligand_cm/dt
    rDY[2] = (kC) / CytosolMembrane; // dcomplex_cm/dt
    rDY[3] = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus; // dbcat_nu/dt
    rDY[4] = (-kN) / nucleus; // dligand_nu/dt
    rDY[5] = (kN) / nucleus; // dcomplex_nu/dt

    // Scale time appropriately
}


template <>
void CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES:
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


    // STATE PARAMETERS:
    this->mParameterNames.push_back("drag");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ComplexTransitThreshold");
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