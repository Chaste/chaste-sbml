#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Gardner1998SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Gardner1998SbmlOdeSystem::Gardner1998SbmlOdeSystem()
        : AbstractSbmlOdeSystem(5, 1, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Gardner1998SbmlOdeSystem>);

    // VARIABLE PARAMETERS
    Cell = 1.0;

    // STATE VARIABLES
    C = 0.0;
    X = 0.0;
    M = 0.0;
    Y = 1.0;
    Z = 1.0;

    // DERIVED QUANTITIES
    V1 = 0.0;
    V3 = 0.0;

    // INITIAL ASSIGNMENTS

    // ODE SYSTEM INFORMATION
    SetDefaultInitialCondition(0, C);
    SetDefaultInitialCondition(1, X);
    SetDefaultInitialCondition(2, M);
    SetDefaultInitialCondition(3, Y);
    SetDefaultInitialCondition(4, Z);

    mStateVariables.push_back(C);
    mStateVariables.push_back(X);
    mStateVariables.push_back(M);
    mStateVariables.push_back(Y);
    mStateVariables.push_back(Z);

    mParameters.push_back(Cell);

    // REACTIONS
    reaction1 = 0.0;
    reaction2 = 0.0;
    reaction3 = 0.0;
    reaction4 = 0.0;
    reaction5 = 0.0;
    reaction6 = 0.0;
    reaction7 = 0.0;
    reaction8 = 0.0;
    reaction9 = 0.0;
    reaction10 = 0.0;
    reaction11 = 0.0;
    reaction12 = 0.0;
    reaction13 = 0.0;

    // EVENTS

    // Run model rules to complete state initialisation
    RunModelRules(0.0, mStateVariables);
}

Gardner1998SbmlOdeSystem::~Gardner1998SbmlOdeSystem()
{
}

std::vector<double> Gardner1998SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    RunModelRules(time, rY);

    std::vector<double> dqs;
    dqs.push_back(V1);
    dqs.push_back(V3);
    return dqs;
}

void Gardner1998SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (reaction1 - reaction2 - reaction3 - reaction8 + reaction9 + reaction10) / Cell; // d[C]/dt
    rDY[1] = (reaction6 - reaction7) / Cell;                                                  // d[X]/dt
    rDY[2] = (reaction4 - reaction5) / Cell;                                                  // d[M]/dt
    rDY[3] = (-reaction8 + reaction9 + reaction11 + reaction12 - reaction13) / Cell;          // d[Y]/dt
    rDY[4] = (reaction8 - reaction9 - reaction10 - reaction11) / Cell;                        // d[Z]/dt

    // TODO: Scale time appropriately
}

double Gardner1998SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

void Gardner1998SbmlOdeSystem::RunModelRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
    C = rY[0];
    X = rY[1];
    M = rY[2];
    Y = rY[3];
    Z = rY[4];

    // VARIABLE PARAMETERS
    Cell = GetParameter(0);

    // ASSIGNMENT RULES
    V1 = C * V1p * std::pow(C + K6, -1.0);
    V3 = M * V3p;

    // REACTIONS
    // creation of cyclin
    reaction1 = 0.0;
    {
        double vi = 0.1;
        reaction1 = vi;
    }

    // cdc2 kinase triggered degration of cyclin
    reaction2 = 0.0;
    {
        double k1 = 0.5;
        double K5 = 0.02;
        reaction2 = C * k1 * X * std::pow(C + K5, -1.0);
    }

    // default degradation of cyclin
    reaction3 = 0.0;
    {
        double kd = 0.02;
        reaction3 = C * kd;
    }

    // activation of cdc2 kinase
    reaction4 = 0.0;
    {
        double K1 = 0.1;
        reaction4 = (1.0 + -1.0 * M) * V1 * std::pow(K1 + -1.0 * M + 1.0, -1.0);
    }

    // deactivation of cdc2 kinase
    reaction5 = 0.0;
    {
        double V2 = 0.25;
        double K2 = 0.1;
        reaction5 = M * V2 * std::pow(K2 + M, -1.0);
    }

    // activation of cyclin protease
    reaction6 = 0.0;
    {
        double K3 = 0.2;
        reaction6 = V3 * (1.0 + -1.0 * X) * std::pow(K3 + -1.0 * X + 1.0, -1.0);
    }

    // deactivation of cyclin protease
    reaction7 = 0.0;
    {
        double K4 = 0.1;
        double V4 = 0.1;
        reaction7 = V4 * X * std::pow(K4 + X, -1.0);
    }

    // reaction8
    reaction8 = 0.0;
    {
        double a1 = 0.05;
        reaction8 = a1 * C * Y;
    }

    // reaction9
    reaction9 = 0.0;
    {
        double a2 = 0.05;
        reaction9 = a2 * Z;
    }

    // desinhibition of cyclin
    reaction10 = 0.0;
    {
        double alpha = 0.1;
        double d1 = 0.05;
        reaction10 = alpha * d1 * Z;
    }

    // degradation of inhibited cyclin
    reaction11 = 0.0;
    {
        double kd = 0.02;
        double alpha = 0.1;
        reaction11 = alpha * kd * Z;
    }

    // creation of cyclin inhibitor
    reaction12 = 0.0;
    {
        double vs = 0.2;
        reaction12 = vs;
    }

    // degradation of cyclin inhibitor
    reaction13 = 0.0;
    {
        double d1 = 0.05;
        reaction13 = d1 * Y;
    }
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Gardner1998SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
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

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("V1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("V3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("Cell");
    this->mParameterUnits.push_back("volume");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Gardner1998SbmlOdeSystem)