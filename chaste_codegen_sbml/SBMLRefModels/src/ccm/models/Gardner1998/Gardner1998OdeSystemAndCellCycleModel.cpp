#include "Gardner1998OdeSystemAndCellCycleModel.hpp"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
Gardner1998OdeSystem::Gardner1998OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(5)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Gardner1998OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.0); // cyclin
    SetDefaultInitialCondition(1, 0.0); // protease
    SetDefaultInitialCondition(2, 0.0); // cdc2k
    SetDefaultInitialCondition(3, 1.0); // cyclin inhibitor
    SetDefaultInitialCondition(4, 1.0); // complex inhibitor-cyclin


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Gardner1998OdeSystem::~Gardner1998OdeSystem()
{
}


void Gardner1998OdeSystem::Init()
{
    /* Initialise model compartments. */
    Cell = 1.0;

    /* Initialise model parameters. */
    V1 = 0.0;
    K6 = 0.3;
    V1p = 0.75;
    V3 = 0.0;
    V3p = 0.3;

}

void Gardner1998OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    /* Define state variables */
    double C = rY[0]; // cyclin
    double X = rY[1]; // protease
    double M = rY[2]; // cdc2k
    double Y = rY[3]; // cyclin inhibitor
    double Z = rY[4]; // complex inhibitor-cyclin

    /* Define state parameters. */


    /* Define algebraic rules. */
    V1 = C * V1p * pow(C + K6, -1);
    V3 = M * V3p;

    /* Define the reactions in this model. */
    // creation of cyclin
    double vi = 0.1; // vi
    double reaction1 = vi;

    // cdc2 kinase triggered degration of cyclin
    double k1 = 0.5; // k1
    double K5 = 0.02; // K5
    double reaction2 = C * k1 * X * pow(C + K5, -1);

    // default degradation of cyclin
    double kd = 0.02; // kd
    double reaction3 = C * kd;

    // activation of cdc2 kinase
    double K1 = 0.1; // K1
    double reaction4 = (1 + -1 * M) * V1 * pow(K1 + -1 * M + 1, -1);

    // deactivation of cdc2 kinase
    double V2 = 0.25; // V2
    double K2 = 0.1; // K2
    double reaction5 = M * V2 * pow(K2 + M, -1);

    // activation of cyclin protease
    double K3 = 0.2; // K3
    double reaction6 = V3 * (1 + -1 * X) * pow(K3 + -1 * X + 1, -1);

    // deactivation of cyclin protease
    double K4 = 0.1; // K4
    double V4 = 0.1; // V4
    double reaction7 = V4 * X * pow(K4 + X, -1);

    // reaction8
    double a1 = 0.05; // a1
    double reaction8 = a1 * C * Y;

    // reaction9
    double a2 = 0.05; // a2
    double reaction9 = a2 * Z;

    // desinhibition of cyclin
    double alpha = 0.1; // alpha
    double d1 = 0.05; // d1
    double reaction10 = alpha * d1 * Z;

    // degradation of inhibited cyclin
    double kd = 0.02; // kd
    double alpha = 0.1; // alpha
    double reaction11 = alpha * kd * Z;

    // creation of cyclin inhibitor
    double vs = 0.2; // vs
    double reaction12 = vs;

    // degradation of cyclin inhibitor
    double d1 = 0.05; // d1
    double reaction13 = d1 * Y;


    rDY[0] = (reaction1 - reaction2 - reaction3 - reaction8 + reaction9 + reaction10) / Cell; // dcyclin/dt
    rDY[1] = (reaction6 - reaction7) / Cell; // dprotease/dt
    rDY[2] = (reaction4 - reaction5) / Cell; // dcdc2k/dt
    rDY[3] = (-reaction8 + reaction9 + reaction11 + reaction12 - reaction13) / Cell; // dcyclin inhibitor/dt
    rDY[4] = (reaction8 - reaction9 - reaction10 - reaction11) / Cell; // dcomplex inhibitor-cyclin/dt

    /* Account for the differences in timescales. */
}


template <>
void CellwiseOdeSystemInformation<Gardner1998OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("protease");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("cdc2k");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("cyclin inhibitor");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("complex inhibitor-cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Gardner1998OdeSystem, 5> Gardner1998CellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Gardner1998OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Gardner1998OdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998CellCycleModel)