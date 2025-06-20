#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Goldbeter1991SbmlOdeSystemAndSrnModel.hpp"

namespace sm = sbmlmath;

Goldbeter1991SbmlOdeSystem::Goldbeter1991SbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>);

    // COMPARTMENTS:
    cell = 1.0;

    // STATE VARIABLES:
    C = 0.01; // Cyclin
    M = 0.01; // cdc_2_kinase
    X = 0.01; // Cyclin Protease

    SetDefaultInitialCondition(0, C);
    SetDefaultInitialCondition(1, M);
    SetDefaultInitialCondition(2, X);

    if (stateVariables.size() == 3)
    {
        C = stateVariables[0];
        M = stateVariables[1];
        X = stateVariables[2];
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("Goldbeter1991SbmlOdeSystem: Expected 3 state variables, got " + std::to_string(stateVariables.size()));
    }

    mStateVariables.push_back(C);
    mStateVariables.push_back(M);
    mStateVariables.push_back(X);

    // STATE PARAMETERS:

    V1 = 0.0;
    V3 = 0.0;


    mParameters.push_back(V1);
    mParameters.push_back(V3);

}

Goldbeter1991SbmlOdeSystem::~Goldbeter1991SbmlOdeSystem()
{
}

void Goldbeter1991SbmlOdeSystem::RefreshState(const std::vector<double> &rY)
{
    // STATE VARIABLES:
    C = rY[0];
    M = rY[1];
    X = rY[2];

    // STATE PARAMETERS:

    V1 = GetParameter("V1");
    V3 = GetParameter("V3");
}

void Goldbeter1991SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RefreshState(rY);

    // RULES:
    V1 = C * VM1 * std::pow(C + Kc, -1.0);
    V3 = M * VM3;

    // UPDATE STATE PARAMETERS:

    SetParameter("V1", V1);
    SetParameter("V3", V3);

    // REACTIONS:

    // creation of cyclin
    double reaction1 = 0.0;
    {
        double vi = 0.025;
        reaction1 = cell * vi;
    }

    // default degradation of cyclin
    double reaction2 = 0.0;
    {
        double kd = 0.01;
        reaction2 = C * cell * kd;
    }

    // cdc2 kinase triggered degration of cyclin
    double reaction3 = 0.0;
    {
        double vd = 0.25;
        double Kd = 0.02;
        reaction3 = C * cell * vd * X * std::pow(C + Kd, -1.0);
    }

    // activation of cdc2 kinase
    double reaction4 = 0.0;
    {
        double K1 = 0.005;
        reaction4 = cell * (1.0 + -1.0 * M) * V1 * std::pow(K1 + -1.0 * M + 1.0, -1.0);
    }

    // deactivation of cdc2 kinase
    double reaction5 = 0.0;
    {
        double V2 = 1.5;
        double K2 = 0.005;
        reaction5 = cell * M * V2 * std::pow(K2 + M, -1.0);
    }

    // activation of cyclin protease
    double reaction6 = 0.0;
    {
        double K3 = 0.005;
        reaction6 = cell * V3 * (1.0 + -1.0 * X) * std::pow(K3 + -1.0 * X + 1.0, -1.0);
    }

    // deactivation of cyclin protease
    double reaction7 = 0.0;
    {
        double K4 = 0.005;
        double V4 = 0.5;
        reaction7 = cell * V4 * X * std::pow(K4 + X, -1.0);
    }

    // ODES:
    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // d[Cyclin]/dt
    rDY[1] = (reaction4 - reaction5) / cell; // d[cdc_2_kinase]/dt
    rDY[2] = (reaction6 - reaction7) / cell; // d[Cyclin Protease]/dt

    // Scale time appropriately
}


// FUNCTION DEFINITIONS:

template <>
void CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES:
    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("M");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    // STATE PARAMETERS:
    this->mParameterNames.push_back("V1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("V3");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Goldbeter1991SbmlOdeSystem, 3> Goldbeter1991SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Goldbeter1991SbmlOdeSystem, 3)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SbmlSrnModel)