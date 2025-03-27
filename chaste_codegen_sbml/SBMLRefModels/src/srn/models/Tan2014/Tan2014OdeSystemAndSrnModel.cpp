#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Tan2014OdeSystemAndSrnModel.hpp"

/* SBML ODE System */
Tan2014OdeSystem::Tan2014OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(7)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014OdeSystem>);

    Init();

<<<<<<< Updated upstream
    SetDefaultInitialCondition(0, 46.6);
    SetDefaultInitialCondition(1, 581.1);
    SetDefaultInitialCondition(2, 418.9);
    SetDefaultInitialCondition(3, 32.6);
    SetDefaultInitialCondition(4, 516.8);
    SetDefaultInitialCondition(5, 483.2);
    SetDefaultInitialCondition(6, 1.0);
=======
    SetDefaultInitialCondition(0, 46.6); // 
    SetDefaultInitialCondition(1, 581.1); // 
    SetDefaultInitialCondition(2, 418.9); // 
    SetDefaultInitialCondition(3, 32.6); // 
    SetDefaultInitialCondition(4, 516.8); // 
    SetDefaultInitialCondition(5, 483.2); // 
    SetDefaultInitialCondition(6, 1.0); // drag
>>>>>>> Stashed changes

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
    /* Initialise the parameters. */
    compartment = 1.0;
    CytosolMembrane = 1.16;
    nucleus = 0.65;
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
<<<<<<< Updated upstream
    double bcat_cm = rY[0];
    double ligand_cm = rY[1];
    double complex_cm = rY[2];
    double bcat_nu = rY[3];
    double ligand_nu = rY[4];
    double complex_nu = rY[5];
=======
    double bcat_cm = rY[0]; // 
    double ligand_cm = rY[1]; // 
    double complex_cm = rY[2]; // 
    double bcat_nu = rY[3]; // 
    double ligand_nu = rY[4]; // 
    double complex_nu = rY[5]; // 
>>>>>>> Stashed changes
    double drag = rY[6]; // drag

    /* Define state parameters */
    double wnt_level = this->mParameters[0];               // wnt_level
    double gamma = this->mParameters[1];                   // gamma
    double ComplexTransitThreshold = this->mParameters[2]; // ComplexTransitThreshold

    /* Define algebraic rules. */
    drag = sbmlmath::sm_max((this->GetStateVariable("complex_cm") - 700) / 10, 1);

    /* Define the reactions in this model. */
    double Bsynthesis = Bsyn * CytosolMembrane;

    double kDegradation = CytosolMembrane * kdegradation * this->GetParameter("gamma") * this->GetStateVariable("bcat_cm") * (1 - 0.5 * this->GetParameter("wnt_level"));

    double kC = CytosolMembrane * (kC_k1 * this->GetStateVariable("bcat_cm") * this->GetStateVariable("ligand_cm") - kC_k2 * this->GetStateVariable("complex_cm"));

    double kN = nucleus * (kN_k1 * this->GetStateVariable("bcat_nu") * this->GetStateVariable("ligand_nu") - kN_k2 * this->GetStateVariable("complex_nu"));

    double kdiffusion = kdiffusion_k * (this->GetStateVariable("bcat_cm") - this->GetStateVariable("bcat_nu"));

    double K_c_active = K_c_active_k * this->GetStateVariable("bcat_cm");

    double K_n_active = K_n_active_k * this->GetStateVariable("bcat_nu");

<<<<<<< Updated upstream
    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // dbcat_cm/dt
    rDY[1] = (-kC) / CytosolMembrane;                                                                   // dligand_cm/dt
    rDY[2] = (kC) / CytosolMembrane;                                                                    // dcomplex_cm/dt
    rDY[3] = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus;                                    // dbcat_nu/dt
    rDY[4] = (-kN) / nucleus;                                                                           // dligand_nu/dt
    rDY[5] = (kN) / nucleus;                                                                            // dcomplex_nu/dt
    rDY[6] = (drag - rY[6]) / CytosolMembrane;                                                          // ddrag/dt
=======

    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // d/dt
    rDY[1] = (-kC) / CytosolMembrane; // d/dt
    rDY[2] = (kC) / CytosolMembrane; // d/dt
    rDY[3] = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus; // d/dt
    rDY[4] = (-kN) / nucleus; // d/dt
    rDY[5] = (kN) / nucleus; // d/dt
    rDY[6] = (drag - rY[6]) / CytosolMembrane; // ddrag/dt
>>>>>>> Stashed changes

    /* Account for the differences in timescales. */
    // rDY[0] *= 60.0;
    // rDY[1] *= 60.0;
    // rDY[2] *= 60.0;
    // rDY[3] *= 60.0;
    // rDY[4] *= 60.0;
    // rDY[5] *= 60.0;
}

template <>
void CellwiseOdeSystemInformation<Tan2014OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(46.6);

    this->mVariableNames.push_back("");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(581.1);

    this->mVariableNames.push_back("");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(418.9);

    this->mVariableNames.push_back("");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(32.6);

    this->mVariableNames.push_back("");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(516.8);

    this->mVariableNames.push_back("");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(483.2);

    this->mVariableNames.push_back("drag");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);

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
