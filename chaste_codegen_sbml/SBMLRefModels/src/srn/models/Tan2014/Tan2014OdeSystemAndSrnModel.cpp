#include "Tan2014OdeSystemAndSrnModel.hpp"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
Tan2014OdeSystem::Tan2014OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(7)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 46.6); // bcat_cm
    SetDefaultInitialCondition(1, 581.1); // ligand_cm
    SetDefaultInitialCondition(2, 418.9); // complex_cm
    SetDefaultInitialCondition(3, 32.6); // bcat_nu
    SetDefaultInitialCondition(4, 516.8); // ligand_nu
    SetDefaultInitialCondition(5, 483.2); // complex_nu
    SetDefaultInitialCondition(6, 1.0); // drag

    this->mParameters.push_back(0.0); // wnt_level
    this->mParameters.push_back(1.0); // gamma
    this->mParameters.push_back(1.0); // ComplexTransitThreshold

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Tan2014OdeSystem::~Tan2014OdeSystem()
{
}


void Tan2014OdeSystem::Init()
{
    /* Initialise model compartments. */
    compartment = 1.0;
    CytosolMembrane = 1.16;
    nucleus = 0.65;

    /* Initialise model parameters. */
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
    gamma = 0.0;
    ComplexTransitThreshold = 1.0;

}

void Tan2014OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    /* Define state variables */
    double bcat_cm = rY[0]; // bcat_cm
    double ligand_cm = rY[1]; // ligand_cm
    double complex_cm = rY[2]; // complex_cm
    double bcat_nu = rY[3]; // bcat_nu
    double ligand_nu = rY[4]; // ligand_nu
    double complex_nu = rY[5]; // complex_nu
    double drag = rY[6]; // drag

    /* Define state parameters. */

    double wnt_level = this->mParameters[0]; // wnt_level
    double gamma = this->mParameters[1]; // gamma
    double ComplexTransitThreshold = this->mParameters[2]; // ComplexTransitThreshold

    /* Define algebraic rules. */
    drag = fmax((complex_cm - 700) / 10, 1);

    /* Define the reactions in this model. */
    double Bsynthesis = Bsyn * CytosolMembrane;

    double kDegradation = CytosolMembrane * kdegradation * gamma * bcat_cm * (1 - 0.5 * wnt_level);

    double kC = CytosolMembrane * (kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm);

    double kN = nucleus * (kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu);

    double kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);

    double K_c_active = K_c_active_k * bcat_cm;

    double K_n_active = K_n_active_k * bcat_nu;


    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // dbcat_cm/dt
    rDY[1] = (-kC) / CytosolMembrane; // dligand_cm/dt
    rDY[2] = (kC) / CytosolMembrane; // dcomplex_cm/dt
    rDY[3] = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus; // dbcat_nu/dt
    rDY[4] = (-kN) / nucleus; // dligand_nu/dt
    rDY[5] = (kN) / nucleus; // dcomplex_nu/dt
    rDY[6] = (drag - rY[6]) / CytosolMembrane; // ddrag/dt

    /* Account for the differences in timescales. */
}


template <>
void CellwiseOdeSystemInformation<Tan2014OdeSystem>::Initialise()
{
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

    this->mVariableNames.push_back("drag");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ComplexTransitThreshold");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014OdeSystem, 7> Tan2014SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Tan2014OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014OdeSystem, 7)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SrnModel)