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

    // AMOUNT / CONCENTRATION CONVERSIONS
    double amt__C = C * cell; //
    double amt__M = M * cell; //
    double amt__X = X * cell; //

    dqs.push_back(cell);
    dqs.push_back(V1);
    dqs.push_back(V3);
    dqs.push_back(reaction1);
    dqs.push_back(reaction2);
    dqs.push_back(reaction3);
    dqs.push_back(reaction4);
    dqs.push_back(reaction5);
    dqs.push_back(reaction6);
    dqs.push_back(reaction7);
    dqs.push_back(amt__C);
    dqs.push_back(amt__M);
    dqs.push_back(amt__X);

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
        [[maybe_unused]] double vi = 0.025;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction1 = cell * vi;
    }
    // reaction2:
    {
        [[maybe_unused]] double kd = 0.01;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction2 = C * cell * kd;
    }
    // reaction3:
    {
        [[maybe_unused]] double vd = 0.25;
        [[maybe_unused]] double Kd = 0.02;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction3 = C * cell * vd * X * std::pow((C + Kd), -1.0);
    }
    // reaction4:
    {
        [[maybe_unused]] double K1 = 0.005;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction4 = cell * (1.0 + -1.0 * M) * V1 * std::pow((K1 + -1.0 * M + 1.0), -1.0);
    }
    // reaction5:
    {
        [[maybe_unused]] double V2 = 1.5;
        [[maybe_unused]] double K2 = 0.005;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction5 = cell * M * V2 * std::pow((K2 + M), -1.0);
    }
    // reaction6:
    {
        [[maybe_unused]] double K3 = 0.005;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction6 = cell * V3 * (1.0 + -1.0 * X) * std::pow((K3 + -1.0 * X + 1.0), -1.0);
    }
    // reaction7:
    {
        [[maybe_unused]] double K4 = 0.005;
        [[maybe_unused]] double V4 = 0.5;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction7 = cell * V4 * X * std::pow((K4 + X), -1.0);
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
        [[maybe_unused]] double vi = 0.025;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction1 = cell * vi;
    }
    // reaction2:
    {
        [[maybe_unused]] double kd = 0.01;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction2 = C * cell * kd;
    }
    // reaction3:
    {
        [[maybe_unused]] double vd = 0.25;
        [[maybe_unused]] double Kd = 0.02;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction3 = C * cell * vd * X * std::pow((C + Kd), -1.0);
    }
    // reaction4:
    {
        [[maybe_unused]] double K1 = 0.005;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction4 = cell * (1.0 + -1.0 * M) * V1 * std::pow((K1 + -1.0 * M + 1.0), -1.0);
    }
    // reaction5:
    {
        [[maybe_unused]] double V2 = 1.5;
        [[maybe_unused]] double K2 = 0.005;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction5 = cell * M * V2 * std::pow((K2 + M), -1.0);
    }
    // reaction6:
    {
        [[maybe_unused]] double K3 = 0.005;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction6 = cell * V3 * (1.0 + -1.0 * X) * std::pow((K3 + -1.0 * X + 1.0), -1.0);
    }
    // reaction7:
    {
        [[maybe_unused]] double K4 = 0.005;
        [[maybe_unused]] double V4 = 0.5;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->reaction7 = cell * V4 * X * std::pow((K4 + X), -1.0);
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

    this->mDerivedQuantityNames.push_back("V1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("V3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction4");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("reaction7");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__M");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__X");
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