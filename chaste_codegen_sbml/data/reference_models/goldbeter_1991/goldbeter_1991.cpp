#include "goldbeter_1991.hpp"
#include "CellwiseOdeSystemInformation.hpp"
/* SBML ODE System */
goldbeter_1991OdeSystem::goldbeter_1991OdeSystem (std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<goldbeter_1991OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.01);
    SetDefaultInitialCondition(1, 0.01);
    SetDefaultInitialCondition(2, 0.01);


    if (stateVariables != std::vector<double>())    {
        SetStateVariables(stateVariables);
    }
}

goldbeter_1991OdeSystem::~goldbeter_1991OdeSystem()
{
}

void goldbeter_1991OdeSystem::Init()
 {
    /* Initialise the parameters. */
    cell = 1.0;
    V1 = 0.0;
    V3 = 0.0;
    VM1 = 3.0;
    VM3 = 1.0;
    Kc = 0.5;
}

void goldbeter_1991OdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    /* Define state variables */
    double C = rY[0]; // Cyclin
    double M = rY[1]; // cdc_2_kinase
    double X = rY[2]; // Cyclin_Protease

    /* Define state parameters */

     /* Define algebraic rules. */
    V1 = C * VM1 * pow(C + Kc, -1);
    V3 = M * VM3;

    /* Define the reactions in this model. */
    // creation of cyclin
    double vi = 0.025;
    double reaction1 = cell * vi;

    // default degradation of cyclin
    double kd = 0.01;
    double reaction2 = C * cell * kd;

    // cdc2 kinase triggered degration of cyclin
    double vd = 0.25;
    double Kd = 0.02;
    double reaction3 = C * cell * vd * X * pow(C + Kd, -1);

    // activation of cdc2 kinase
    double K1 = 0.005;
    double reaction4 = cell * (1 + -1 * M) * V1 * pow(K1 + -1 * M + 1, -1);

    // deactivation of cdc2 kinase
    double V2 = 1.5;
    double K2 = 0.005;
    double reaction5 = cell * M * V2 * pow(K2 + M, -1);

    // activation of cyclin protease
    double K3 = 0.005;
    double reaction6 = cell * V3 * (1 + -1 * X) * pow(K3 + -1 * X + 1, -1);

    // deactivation of cyclin protease
    double K4 = 0.005;
    double V4 = 0.5;
    double reaction7 = cell * V4 * X * pow(K4 + X, -1);


    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // dCyclin/dt
    rDY[1] = (reaction4 - reaction5) / cell; // dcdc_2_kinase/dt
    rDY[2] = (reaction6 - reaction7) / cell; // dCyclin Protease/dt

    /* Account for the differences in timescales. */
    rDY[0] *= 3600.0;
    rDY[1] *= 3600.0;
    rDY[2] *= 3600.0;

}

template<>
void CellwiseOdeSystemInformation<goldbeter_1991OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("cdc_2_kinase");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("Cyclin Protease");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<goldbeter_1991OdeSystem,5> goldbeter_1991;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(goldbeter_1991OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, goldbeter_1991OdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(goldbeter_1991)

