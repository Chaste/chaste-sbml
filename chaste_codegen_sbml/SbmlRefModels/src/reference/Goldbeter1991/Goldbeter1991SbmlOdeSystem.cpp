#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Goldbeter1991SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Goldbeter1991SbmlOdeSystem::Goldbeter1991SbmlOdeSystem()
        : AbstractSbmlOdeSystem(3, 1, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>);

    // STATE VARIABLES
    C = 0.01;
    M = 0.01;
    X = 0.01;

    // DERIVED QUANTITIES

    // VARIABLE PARAMETERS
    cell = 1.0;

    // RULE-BASED PARAMETERS
    V1 = 0.0;
    V3 = 0.0;

    // INITIAL ASSIGNMENTS

    // ODE SYSTEM INFORMATION
    SetDefaultInitialCondition(0, C);
    SetDefaultInitialCondition(1, M);
    SetDefaultInitialCondition(2, X);

    mStateVariables.push_back(C);
    mStateVariables.push_back(M);
    mStateVariables.push_back(X);

    mParameters.push_back(cell);

    // REACTIONS
    reaction1 = 0.0;
    reaction2 = 0.0;
    reaction3 = 0.0;
    reaction4 = 0.0;
    reaction5 = 0.0;
    reaction6 = 0.0;
    reaction7 = 0.0;

    // EVENTS

    // Run model rules to complete state initialisation
    RunModelRules(0.0, mStateVariables);
}

Goldbeter1991SbmlOdeSystem::~Goldbeter1991SbmlOdeSystem()
{
}

std::vector<double> Goldbeter1991SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    RunModelRules(time, rY);

    std::vector<double> dqs;
    return dqs;
}

void Goldbeter1991SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // d[C]/dt
    rDY[1] = (reaction4 - reaction5) / cell;             // d[M]/dt
    rDY[2] = (reaction6 - reaction7) / cell;             // d[X]/dt

    // TODO: Scale time appropriately
}

double Goldbeter1991SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

void Goldbeter1991SbmlOdeSystem::RunModelRules(double time, const std::vector<double>& rY)
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

// MODEL FUNCTIONS

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

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlOdeSystem)