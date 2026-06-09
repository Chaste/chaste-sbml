#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Gardner1998SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Gardner1998SbmlOdeSystem::Gardner1998SbmlOdeSystem()
        : AbstractSbmlOdeSystem(5, 3, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Gardner1998SbmlOdeSystem>);

    Initialise();

    // EVENTS
}

Gardner1998SbmlOdeSystem::~Gardner1998SbmlOdeSystem()
{
}

std::vector<double> Gardner1998SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNT / CONCENTRATION CONVERSIONS
    double amt__C = C * Cell; //
    double amt__X = X * Cell; //
    double amt__M = M * Cell; //
    double amt__Y = Y * Cell; //
    double amt__Z = Z * Cell; //

    dqs.push_back(Cell);
    dqs.push_back(V1);
    dqs.push_back(V3);
    dqs.push_back(amt__C);
    dqs.push_back(amt__X);
    dqs.push_back(amt__M);
    dqs.push_back(amt__Y);
    dqs.push_back(amt__Z);

    return dqs;
}

void Gardner1998SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void Gardner1998SbmlOdeSystem::Initialise(double time)
{
    Cell = 1.0;                              //
    C = 0.0;                                 //
    X = 0.0;                                 //
    M = 0.0;                                 //
    Y = 1.0;                                 //
    Z = 1.0;                                 //
    K6 = 0.3;                                //
    V1p = 0.75;                              //
    V3p = 0.3;                               //
    C = C / Cell;                            //
    X = X / Cell;                            //
    M = M / Cell;                            //
    Y = Y / Cell;                            //
    Z = Z / Cell;                            //
    V1 = C * V1p * std::pow((C + K6), -1.0); //
    V3 = M * V3p;                            //
    // reaction1:
    {
        double vi = 0.1;
        reaction1 = vi;
    }
    // reaction2:
    {
        double k1 = 0.5;
        double K5 = 0.02;
        reaction2 = C * k1 * X * std::pow((C + K5), -1.0);
    }
    // reaction3:
    {
        double kd = 0.02;
        reaction3 = C * kd;
    }
    // reaction4:
    {
        double K1 = 0.1;
        reaction4 = (1.0 + -1.0 * M) * V1 * std::pow((K1 + -1.0 * M + 1.0), -1.0);
    }
    // reaction5:
    {
        double V2 = 0.25;
        double K2 = 0.1;
        reaction5 = M * V2 * std::pow((K2 + M), -1.0);
    }
    // reaction6:
    {
        double K3 = 0.2;
        reaction6 = V3 * (1.0 + -1.0 * X) * std::pow((K3 + -1.0 * X + 1.0), -1.0);
    }
    // reaction7:
    {
        double K4 = 0.1;
        double V4 = 0.1;
        reaction7 = V4 * X * std::pow((K4 + X), -1.0);
    }
    // reaction8:
    {
        double a1 = 0.05;
        reaction8 = a1 * C * Y;
    }
    // reaction9:
    {
        double a2 = 0.05;
        reaction9 = a2 * Z;
    }
    // reaction10:
    {
        double alpha = 0.1;
        double d1 = 0.05;
        reaction10 = alpha * d1 * Z;
    }
    // reaction11:
    {
        double kd = 0.02;
        double alpha = 0.1;
        reaction11 = alpha * kd * Z;
    }
    // reaction12:
    {
        double vs = 0.2;
        reaction12 = vs;
    }
    // reaction13:
    {
        double d1 = 0.05;
        reaction13 = d1 * Y;
    }
    d_C_dt = ((reaction1 - reaction2 - reaction3 - reaction8) + reaction9 + reaction10) / Cell; //
    d_X_dt = (reaction6 - reaction7) / Cell;                                                    //
    d_M_dt = (reaction4 - reaction5) / Cell;                                                    //
    d_Y_dt = (-reaction8 + reaction9 + reaction11 + reaction12 - reaction13) / Cell;            //
    d_Z_dt = (reaction8 - reaction9 - reaction10 - reaction11) / Cell;                          //

    mStateVariables.push_back(C);
    mStateVariables.push_back(X);
    mStateVariables.push_back(M);
    mStateVariables.push_back(Y);
    mStateVariables.push_back(Z);

    SetDefaultInitialCondition(0, C);
    SetDefaultInitialCondition(1, X);
    SetDefaultInitialCondition(2, M);
    SetDefaultInitialCondition(3, Y);
    SetDefaultInitialCondition(4, Z);

    mParameters.push_back(K6);
    mParameters.push_back(V1p);
    mParameters.push_back(V3p);
}

double Gardner1998SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    // Ensure all member variables (state vars, parameters, derived quantities) reflect
    // the rY passed in. Without this, event triggers and assignments would use stale
    // values from the last EvaluateYDerivatives call, which may differ from rY when
    // called from CalculateRootFunction or CalculateStoppingEvent with a different state.
    RunModelEquations(time, rY);

    // Do NOT clear mEventAdjustedStateVars/Parameters here. Once set by an event fire,
    // they must persist across all CVODE bisection calls until AdjustParameters applies
    // them. Clearing here would erase the stored assignment when a later bisection call
    // lands in the clamped state (mEventSatisfied=true), causing the halving to be lost.
    // CalculateStoppingEvent (BackwardEuler path) clears these itself before calling.

    // Root function for CVODE: the maximum signed event distance, where each distance is
    // positive exactly when its event's trigger condition holds. Taking the MAXIMUM (not the
    // minimum absolute value) means the combined function crosses zero the moment ANY event
    // becomes triggered, and cannot be masked by another event that happens to sit just below
    // its own boundary (a small negative distance). A min-abs combination misses an event
    // whose rising edge coincides with another event re-arming near its threshold.
    double max_dist = -std::numeric_limits<double>::max();

    return max_dist; // Signed distance of the event closest to triggering
}

// ASSIGNMENT RULES
void Gardner1998SbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void Gardner1998SbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> Gardner1998SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    C = rStateVariables[0];
    X = rStateVariables[1];
    M = rStateVariables[2];
    Y = rStateVariables[3];
    Z = rStateVariables[4];

    K6 = GetParameter(0);
    V1p = GetParameter(1);
    V3p = GetParameter(2);

    V1 = C * V1p * std::pow((C + K6), -1.0); //
    V3 = M * V3p;                            //
    // reaction1:
    {
        double vi = 0.1;
        reaction1 = vi;
    }
    // reaction2:
    {
        double k1 = 0.5;
        double K5 = 0.02;
        reaction2 = C * k1 * X * std::pow((C + K5), -1.0);
    }
    // reaction3:
    {
        double kd = 0.02;
        reaction3 = C * kd;
    }
    // reaction4:
    {
        double K1 = 0.1;
        reaction4 = (1.0 + -1.0 * M) * V1 * std::pow((K1 + -1.0 * M + 1.0), -1.0);
    }
    // reaction5:
    {
        double V2 = 0.25;
        double K2 = 0.1;
        reaction5 = M * V2 * std::pow((K2 + M), -1.0);
    }
    // reaction6:
    {
        double K3 = 0.2;
        reaction6 = V3 * (1.0 + -1.0 * X) * std::pow((K3 + -1.0 * X + 1.0), -1.0);
    }
    // reaction7:
    {
        double K4 = 0.1;
        double V4 = 0.1;
        reaction7 = V4 * X * std::pow((K4 + X), -1.0);
    }
    // reaction8:
    {
        double a1 = 0.05;
        reaction8 = a1 * C * Y;
    }
    // reaction9:
    {
        double a2 = 0.05;
        reaction9 = a2 * Z;
    }
    // reaction10:
    {
        double alpha = 0.1;
        double d1 = 0.05;
        reaction10 = alpha * d1 * Z;
    }
    // reaction11:
    {
        double kd = 0.02;
        double alpha = 0.1;
        reaction11 = alpha * kd * Z;
    }
    // reaction12:
    {
        double vs = 0.2;
        reaction12 = vs;
    }
    // reaction13:
    {
        double d1 = 0.05;
        reaction13 = d1 * Y;
    }
    d_C_dt = ((reaction1 - reaction2 - reaction3 - reaction8) + reaction9 + reaction10) / Cell; //
    d_X_dt = (reaction6 - reaction7) / Cell;                                                    //
    d_M_dt = (reaction4 - reaction5) / Cell;                                                    //
    d_Y_dt = (-reaction8 + reaction9 + reaction11 + reaction12 - reaction13) / Cell;            //
    d_Z_dt = (reaction8 - reaction9 - reaction10 - reaction11) / Cell;                          //

    std::vector<double> derivatives(5);
    derivatives[0] = d_C_dt;
    derivatives[1] = d_X_dt;
    derivatives[2] = d_M_dt;
    derivatives[3] = d_Y_dt;
    derivatives[4] = d_Z_dt;
    return derivatives;
}

// REACTIONS
void Gardner1998SbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void Gardner1998SbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void Gardner1998SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    C = rStateVariables[0];
    X = rStateVariables[1];
    M = rStateVariables[2];
    Y = rStateVariables[3];
    Z = rStateVariables[4];
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
    this->mDerivedQuantityNames.push_back("Cell");
    this->mDerivedQuantityUnits.push_back("volume");

    this->mDerivedQuantityNames.push_back("V1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("V3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__X");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__M");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__Y");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__Z");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("K6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("V1p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("V3p");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Gardner1998SbmlOdeSystem)