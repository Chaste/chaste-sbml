#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Semantic00928L3V2SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Semantic00928L3V2SbmlOdeSystem::Semantic00928L3V2SbmlOdeSystem()
        : AbstractSbmlOdeSystem(2, 1, 1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Semantic00928L3V2SbmlOdeSystem>);

    Initialise();

    // EVENTS
    mEventType.resize(1, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    // mEventType[0] = SbmlEventType::CELL_DIVISION; //

    mEventSatisfied = { false }; // From SBML trigger initialValue
    mEventTriggered.resize(1, false);

    mEventAdjustedParameters.resize(1, false);
    mEventAdjustedParameterValues.resize(1, 0.0);

    mEventAdjustedStateVars.resize(2, false);
    mEventAdjustedStateValues.resize(2, 0.0);

    // SBML trigger initialValue="false": fire this event at t=0 if its trigger is true.
    {
        double time = 0.0;
        if (time >= 0.0)
        {
            S1 = 0.00015 / C;
            SetStateVariable(0, S1);
            SetDefaultInitialCondition(0, S1);
            mEventSatisfied[0] = true;
        }
    }
}

Semantic00928L3V2SbmlOdeSystem::~Semantic00928L3V2SbmlOdeSystem()
{
}

std::vector<double> Semantic00928L3V2SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNTS
    double amt__S1 = S1 * C; //
    double amt__S2 = S2 * C; //

    dqs.push_back(C);
    dqs.push_back(amt__S1);
    dqs.push_back(amt__S2);

    return dqs;
}

void Semantic00928L3V2SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void Semantic00928L3V2SbmlOdeSystem::Initialise(double time)
{
    C = 1.0;                  //
    S1 = 0.0;                 //
    S2 = 0.0;                 //
    k1 = 1.0;                 //
    S1 = S1 / C;              //
    S2 = S2 / C;              //
    reaction1 = C * k1 * S1;  //
    d_S1_dt = -reaction1 / C; //
    d_S2_dt = reaction1 / C;  //

    mStateVariables.push_back(S1);
    mStateVariables.push_back(S2);

    SetDefaultInitialCondition(0, S1);
    SetDefaultInitialCondition(1, S2);

    mParameters.push_back(k1);
}

double Semantic00928L3V2SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    //========================================
    // EVENT:
    //========================================
    {
        double event_dist = (time) - (0.0) - std::numeric_limits<double>::epsilon();

        // Once an event has fired and its trigger remains active, force a large positive
        // distance so CVODE does not keep detecting the same root.
        if (mEventSatisfied[0] && (time >= 0.0))
        {
            event_dist = std::abs(event_dist) + 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if (time >= 0.0)
        {
            if (!mEventSatisfied[0])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[0] = true;

                // Adjust relevant state variables and parameters
                // S1 = 0.00015 / C
                mEventAdjustedStateVars[0] = true;
                mEventAdjustedStateValues[0] = 0.00015 / C;
            }
            mEventSatisfied[0] = true;
        }
        else
        {
            mEventSatisfied[0] = false;
            mEventTriggered[0] = false;
        }
    }

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void Semantic00928L3V2SbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void Semantic00928L3V2SbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> Semantic00928L3V2SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    S1 = rStateVariables[0];
    S2 = rStateVariables[1];

    k1 = GetParameter(0);

    reaction1 = C * k1 * S1;  //
    d_S1_dt = -reaction1 / C; //
    d_S2_dt = reaction1 / C;  //

    std::vector<double> derivatives(2);
    derivatives[0] = d_S1_dt;
    derivatives[1] = d_S2_dt;
    return derivatives;
}

// REACTIONS
void Semantic00928L3V2SbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void Semantic00928L3V2SbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void Semantic00928L3V2SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    S1 = rStateVariables[0];
    S2 = rStateVariables[1];
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Semantic00928L3V2SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("S1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("S2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("C");
    this->mDerivedQuantityUnits.push_back("litre");

    this->mDerivedQuantityNames.push_back("amt__S1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__S2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("k1");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Semantic00928L3V2SbmlOdeSystem)