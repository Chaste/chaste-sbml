#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Goldbeter1991SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Goldbeter1991SbmlOdeSystem::Goldbeter1991SbmlOdeSystem()
        : AbstractSbmlOdeSystem(3, 3, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>);

    Initialise();

    // EVENTS
}

Goldbeter1991SbmlOdeSystem::~Goldbeter1991SbmlOdeSystem()
{
}

std::vector<double> Goldbeter1991SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNTS
    double amt__C = C * cell; //
    double amt__M = M * cell; //
    double amt__X = X * cell; //

    dqs.push_back(cell);
    dqs.push_back(amt__C);
    dqs.push_back(amt__M);
    dqs.push_back(amt__X);
    dqs.push_back(V1);
    dqs.push_back(V3);

    return dqs;
}

void Goldbeter1991SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void Goldbeter1991SbmlOdeSystem::Initialise(double time)
{
    cell = 1.0;                              //
    C = 0.01;                                //
    M = 0.01;                                //
    X = 0.01;                                //
    VM1 = 3.0;                               //
    VM3 = 1.0;                               //
    Kc = 0.5;                                //
    V1 = C * VM1 * std::pow((C + Kc), -1.0); //
    V3 = M * VM3;                            //
    // reaction1:
    {
        double vi = 0.025;
        reaction1 = cell * vi;
    }
    // reaction2:
    {
        double kd = 0.01;
        reaction2 = C * cell * kd;
    }
    // reaction3:
    {
        double vd = 0.25;
        double Kd = 0.02;
        reaction3 = C * cell * vd * X * std::pow((C + Kd), -1.0);
    }
    // reaction4:
    {
        double K1 = 0.005;
        reaction4 = cell * (1.0 + -1.0 * M) * V1 * std::pow((K1 + -1.0 * M + 1.0), -1.0);
    }
    // reaction5:
    {
        double V2 = 1.5;
        double K2 = 0.005;
        reaction5 = cell * M * V2 * std::pow((K2 + M), -1.0);
    }
    // reaction6:
    {
        double K3 = 0.005;
        reaction6 = cell * V3 * (1.0 + -1.0 * X) * std::pow((K3 + -1.0 * X + 1.0), -1.0);
    }
    // reaction7:
    {
        double K4 = 0.005;
        double V4 = 0.5;
        reaction7 = cell * V4 * X * std::pow((K4 + X), -1.0);
    }
    d_C_dt = (reaction1 - reaction2 - reaction3) / cell; //
    d_M_dt = (reaction4 - reaction5) / cell;             //
    d_X_dt = (reaction6 - reaction7) / cell;             //

    mStateVariables.push_back(C);
    mStateVariables.push_back(M);
    mStateVariables.push_back(X);

    SetDefaultInitialCondition(0, C);
    SetDefaultInitialCondition(1, M);
    SetDefaultInitialCondition(2, X);

    mParameters.push_back(VM1);
    mParameters.push_back(VM3);
    mParameters.push_back(Kc);
}

double Goldbeter1991SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void Goldbeter1991SbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void Goldbeter1991SbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> Goldbeter1991SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    C = rStateVariables[0];
    M = rStateVariables[1];
    X = rStateVariables[2];

    VM1 = GetParameter(0);
    VM3 = GetParameter(1);
    Kc = GetParameter(2);

    V1 = C * VM1 * std::pow((C + Kc), -1.0); //
    V3 = M * VM3;                            //
    // reaction1:
    {
        double vi = 0.025;
        reaction1 = cell * vi;
    }
    // reaction2:
    {
        double kd = 0.01;
        reaction2 = C * cell * kd;
    }
    // reaction3:
    {
        double vd = 0.25;
        double Kd = 0.02;
        reaction3 = C * cell * vd * X * std::pow((C + Kd), -1.0);
    }
    // reaction4:
    {
        double K1 = 0.005;
        reaction4 = cell * (1.0 + -1.0 * M) * V1 * std::pow((K1 + -1.0 * M + 1.0), -1.0);
    }
    // reaction5:
    {
        double V2 = 1.5;
        double K2 = 0.005;
        reaction5 = cell * M * V2 * std::pow((K2 + M), -1.0);
    }
    // reaction6:
    {
        double K3 = 0.005;
        reaction6 = cell * V3 * (1.0 + -1.0 * X) * std::pow((K3 + -1.0 * X + 1.0), -1.0);
    }
    // reaction7:
    {
        double K4 = 0.005;
        double V4 = 0.5;
        reaction7 = cell * V4 * X * std::pow((K4 + X), -1.0);
    }
    d_C_dt = (reaction1 - reaction2 - reaction3) / cell; //
    d_M_dt = (reaction4 - reaction5) / cell;             //
    d_X_dt = (reaction6 - reaction7) / cell;             //

    std::vector<double> derivatives(3);
    derivatives[0] = d_C_dt;
    derivatives[1] = d_M_dt;
    derivatives[2] = d_X_dt;
    return derivatives;
}

// REACTIONS
void Goldbeter1991SbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void Goldbeter1991SbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void Goldbeter1991SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    C = rStateVariables[0];
    M = rStateVariables[1];
    X = rStateVariables[2];
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
    this->mDerivedQuantityNames.push_back("cell");
    this->mDerivedQuantityUnits.push_back("volume");

    this->mDerivedQuantityNames.push_back("amt__C");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__M");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__X");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("V1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("V3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("VM1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("VM3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Kc");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlOdeSystem)