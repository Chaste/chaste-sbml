#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Semantic00991L2V5SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Semantic00991L2V5SbmlOdeSystem::Semantic00991L2V5SbmlOdeSystem()
        : AbstractSbmlOdeSystem(1, 3, 1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Semantic00991L2V5SbmlOdeSystem>);

    Initialise();

    // EVENTS
    mEventType.resize(1, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    // mEventType[0] = SbmlEventType::CELL_DIVISION; //

    mEventSatisfied = { true }; // From SBML trigger initialValue
    mEventTriggered.resize(1, false);

    mEventAdjustedParameters.resize(3, false);
    mEventAdjustedParameterValues.resize(3, 0.0);

    mEventAdjustedStateVars.resize(1, false);
    mEventAdjustedStateValues.resize(1, 0.0);
}

Semantic00991L2V5SbmlOdeSystem::~Semantic00991L2V5SbmlOdeSystem()
{
}

std::vector<double> Semantic00991L2V5SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNTS
    double amt__X = X * default_compartment; //

    dqs.push_back(default_compartment);
    dqs.push_back(amt__X);

    return dqs;
}

void Semantic00991L2V5SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void Semantic00991L2V5SbmlOdeSystem::Initialise(double time)
{
    default_compartment = 1.0;                //
    X = 1.0;                                  //
    p1 = 1.0;                                 //
    k1 = 1.0;                                 //
    Xref = p1;                                //
    J0 = k1;                                  //
    d_X_dt = Xref * J0 / default_compartment; //

    mStateVariables.push_back(X);

    SetDefaultInitialCondition(0, X);

    mParameters.push_back(Xref);
    mParameters.push_back(p1);
    mParameters.push_back(k1);
}

double Semantic00991L2V5SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    //========================================
    // EVENT:
    //========================================
    {
        double event_dist = (X) - (2.0) - std::numeric_limits<double>::epsilon();

        // Once an event has fired and its trigger remains active, force a large positive
        // distance so CVODE does not keep detecting the same root.
        if (mEventSatisfied[0] && (X >= 2.0))
        {
            event_dist = std::abs(event_dist) + 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if (X >= 2.0)
        {
            if (!mEventSatisfied[0])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[0] = true;

                // Adjust relevant state variables and parameters
                // p1 = 2.0
                mEventAdjustedParameters[1] = true;
                mEventAdjustedParameterValues[1] = 2.0;
                SetParameter(1, 2.0);
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
void Semantic00991L2V5SbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void Semantic00991L2V5SbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> Semantic00991L2V5SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    X = rStateVariables[0];

    Xref = GetParameter(0);
    p1 = GetParameter(1);
    k1 = GetParameter(2);

    Xref = p1;                                //
    J0 = k1;                                  //
    d_X_dt = Xref * J0 / default_compartment; //

    std::vector<double> derivatives(1);
    derivatives[0] = d_X_dt;
    return derivatives;
}

// REACTIONS
void Semantic00991L2V5SbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void Semantic00991L2V5SbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void Semantic00991L2V5SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    X = rStateVariables[0];
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Semantic00991L2V5SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("default_compartment");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__X");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("Xref");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("p1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k1");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Semantic00991L2V5SbmlOdeSystem)