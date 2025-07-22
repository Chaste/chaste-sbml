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

    // STATE VARIABLES
    C = 0.01;
    M = 0.01;
    X = 0.01;

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

    // DERIVED QUANTITIES

    // VARIABLE PARAMETERS
    cell = 1.0;

    mParameters.push_back(cell);

    // RULE-BASED PARAMETERS
    V1 = 0.0;
    V3 = 0.0;

    // REACTIONS
    reaction1 = 0.0;
    reaction2 = 0.0;
    reaction3 = 0.0;
    reaction4 = 0.0;
    reaction5 = 0.0;
    reaction6 = 0.0;
    reaction7 = 0.0;


    ProcessRules(0.0, mStateVariables);
}

Goldbeter1991SbmlOdeSystem::~Goldbeter1991SbmlOdeSystem()
{
}

void Goldbeter1991SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    ProcessRules(time, rY);

    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // d[C]/dt
    rDY[1] = (reaction4 - reaction5) / cell; // d[M]/dt
    rDY[2] = (reaction6 - reaction7) / cell; // d[X]/dt

    // Scale time appropriately
}


void Goldbeter1991SbmlOdeSystem::ProcessRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
    C = rY[0];
    M = rY[1];
    X = rY[2];

    // VARIABLE PARAMETERS
    cell = GetParameter(0);

    // ASSIGNMENT RULES
    V1 = C * VM1 * std::pow(C + Kc, -1.0);
    V3 = M * VM3;

    // REACTIONS
    // creation of cyclin
    reaction1 = 0.0;
    {
        double vi = 0.025;
        reaction1 = cell * vi;
    }

    // default degradation of cyclin
    reaction2 = 0.0;
    {
        double kd = 0.01;
        reaction2 = C * cell * kd;
    }

    // cdc2 kinase triggered degration of cyclin
    reaction3 = 0.0;
    {
        double vd = 0.25;
        double Kd = 0.02;
        reaction3 = C * cell * vd * X * std::pow(C + Kd, -1.0);
    }

    // activation of cdc2 kinase
    reaction4 = 0.0;
    {
        double K1 = 0.005;
        reaction4 = cell * (1.0 + -1.0 * M) * V1 * std::pow(K1 + -1.0 * M + 1.0, -1.0);
    }

    // deactivation of cdc2 kinase
    reaction5 = 0.0;
    {
        double V2 = 1.5;
        double K2 = 0.005;
        reaction5 = cell * M * V2 * std::pow(K2 + M, -1.0);
    }

    // activation of cyclin protease
    reaction6 = 0.0;
    {
        double K3 = 0.005;
        reaction6 = cell * V3 * (1.0 + -1.0 * X) * std::pow(K3 + -1.0 * X + 1.0, -1.0);
    }

    // deactivation of cyclin protease
    reaction7 = 0.0;
    {
        double K4 = 0.005;
        double V4 = 0.5;
        reaction7 = cell * V4 * X * std::pow(K4 + X, -1.0);
    }

}




// FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("M");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    // DERIVED QUANTITIES

    // PARAMETERS
    this->mParameterNames.push_back("cell");
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