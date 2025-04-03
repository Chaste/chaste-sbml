#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Gardner1998OdeSystemAndCellCycleModel.hpp"

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
    /* Define algebraic rules. */
    V1 = rY[0] * V1p * std::pow(rY[0] + K6, -1);
    V3 = rY[2] * V3p;

    /* Define the reactions in this model. */
    // creation of cyclin
    double reaction1 = 0.0;
    {
        double vi = 0.1;
        reaction1 = vi;
    }

    // cdc2 kinase triggered degration of cyclin
    double reaction2 = 0.0;
    {
        double k1 = 0.5;
        double K5 = 0.02;
        reaction2 = rY[0] * k1 * rY[1] * std::pow(rY[0] + K5, -1);
    }

    // default degradation of cyclin
    double reaction3 = 0.0;
    {
        double kd = 0.02;
        reaction3 = rY[0] * kd;
    }

    // activation of cdc2 kinase
    double reaction4 = 0.0;
    {
        double K1 = 0.1;
        reaction4 = (1 + -1 * rY[2]) * V1 * std::pow(K1 + -1 * rY[2] + 1, -1);
    }

    // deactivation of cdc2 kinase
    double reaction5 = 0.0;
    {
        double V2 = 0.25;
        double K2 = 0.1;
        reaction5 = rY[2] * V2 * std::pow(K2 + rY[2], -1);
    }

    // activation of cyclin protease
    double reaction6 = 0.0;
    {
        double K3 = 0.2;
        reaction6 = V3 * (1 + -1 * rY[1]) * std::pow(K3 + -1 * rY[1] + 1, -1);
    }

    // deactivation of cyclin protease
    double reaction7 = 0.0;
    {
        double K4 = 0.1;
        double V4 = 0.1;
        reaction7 = V4 * rY[1] * std::pow(K4 + rY[1], -1);
    }

    // reaction8
    double reaction8 = 0.0;
    {
        double a1 = 0.05;
        reaction8 = a1 * rY[0] * rY[3];
    }

    // reaction9
    double reaction9 = 0.0;
    {
        double a2 = 0.05;
        reaction9 = a2 * rY[4];
    }

    // desinhibition of cyclin
    double reaction10 = 0.0;
    {
        double alpha = 0.1;
        double d1 = 0.05;
        reaction10 = alpha * d1 * rY[4];
    }

    // degradation of inhibited cyclin
    double reaction11 = 0.0;
    {
        double kd = 0.02;
        double alpha = 0.1;
        reaction11 = alpha * kd * rY[4];
    }

    // creation of cyclin inhibitor
    double reaction12 = 0.0;
    {
        double vs = 0.2;
        reaction12 = vs;
    }

    // degradation of cyclin inhibitor
    double reaction13 = 0.0;
    {
        double d1 = 0.05;
        reaction13 = d1 * rY[3];
    }


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
    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("M");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Y");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Z");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SbmlCellCycleWrapperModel using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Gardner1998OdeSystem, 5> Gardner1998CellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Gardner1998OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Gardner1998OdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998CellCycleModel)