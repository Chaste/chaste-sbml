#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "{{ ode_hpp_file }}"

namespace sm = sbmlmath;

{{ ode_class_name }}::{{ ode_class_name }}()
    : AbstractSbmlOdeSystem({{ state_variables|length }}, {{parameters|length }}, {{ events|length }})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ ode_class_name }}>);

    Initialise();

    // EVENTS
{% if events %}
    mEventType.resize({{ events|length }}, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

{% for event in events %}
{% if event["type"] == EventType.CELL_DIVISION %}
    mEventType[{{ event["index"] }}] = SbmlEventType::CELL_DIVISION; // {{ event["label"] }}
{% else %}
    // mEventType[{{ event["index"] }}] = SbmlEventType::CELL_DIVISION; // {{ event["label"] }}
{% endif %}
{% endfor %}

    mEventSatisfied = { {% for event in events %}{{ "true" if event["initial_satisfied"] else "false" }}{% if not loop.last %}, {% endif %}{% endfor %} }; // From SBML trigger initialValue
    mEventTriggered.resize({{ events|length }}, false);

    mEventAdjustedParameters.resize({{ parameters|length }}, false);
    mEventAdjustedParameterValues.resize({{ parameters|length }}, 0.0);

    mEventAdjustedStateVars.resize({{ state_variables|length }}, false);
    mEventAdjustedStateValues.resize({{ state_variables|length }}, 0.0);

{% for event in events %}
{% if not event["initial_satisfied"] %}
    // SBML trigger initialValue="false": fire this event at t=0 if its trigger is true.
    {
        [[maybe_unused]] double time = 0.0; // May be unused if no expression references time
        if ({{ event["trigger"] }})
        {
{% for assignment in event["assignments"] %}
{% if ( assignment["type"] == VarType.STATE_VARIABLE ) %}
            {{ assignment["lhs"] }} = {{ assignment["rhs"] }};
            SetStateVariable({{ assignment["index"] }}, {{ assignment["lhs"] }});
            SetDefaultInitialCondition({{ assignment["index"] }}, {{ assignment["lhs"] }});
{% elif ( assignment["type"] == VarType.PARAMETER ) %}
            {{ assignment["lhs"] }} = {{ assignment["rhs"] }};
            SetParameter({{ assignment["index"] }}, {{ assignment["lhs"] }});
{% else %}
            {{ assignment["lhs"] }} = {{ assignment["rhs"] }};
{% endif %}
{% endfor %}
            mEventSatisfied[{{ event["index"] }}] = true;
        }
    }
{% endif %}
{% endfor %}
{% endif %} {# 'if events' #}
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

std::vector<double> {{ ode_class_name }}::ComputeDerivedQuantities(double time, const std::vector<double> &rY)
{
    std::vector<double> dqs;
{% if derived_quantities %}
    dqs.reserve({{ derived_quantities | length }});
{% if scale_time %}
    time *= {{ time_multiplier }}; // Chaste integrates in hours; use the model's native time units
{% endif %}
    RunModelEquations(time, rY);

    // AMOUNT / CONCENTRATION CONVERSIONS
{% for eq in equations %}
{% if ( eq["type"] == EquationType.CONVERSION ) %}
    double {{ eq["var"] }} = {{ eq["rhs"] }}; // {{ eq["label"] }}
{% endif %}
{% endfor %}

{% for dq in derived_quantities %}
    dqs.push_back({{ dq["id"] }});
{% endfor %}
{% endif %} {# 'if derived_quantities' #}

    return dqs;
} // LCOV_EXCL_LINE (gcov marks this closing brace of a std::vector-returning function as uncovered)

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
{% if scale_time %}
    // Chaste integrates in hours; convert to the model's native time units and scale the
    // resulting derivatives (dY/d(hours) = {{ time_multiplier }} * dY/d(native)).
    time *= {{ time_multiplier }};
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = {{ time_multiplier }} * derivatives[i];
    }
{% else %}
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
{% endif %}
}

void {{ ode_class_name }}::Initialise(double time)
{
{% for eq in equations %}
{% if ( eq["type"] != EquationType.CONVERSION ) %}
{% if eq["local_parameters"] %}
    // {{ eq["var"] }}: {{ eq["label"] }}
    {
{% for local_parameter in eq["local_parameters"] %}
        [[maybe_unused]] double {{ local_parameter["id"] }} = {{ local_parameter["value"] }};
{% endfor %}
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->{{ eq["var"] }} = {{ eq["rhs"] }};
    }
{% else %}
    {{ eq["var"] }} = {{ eq["rhs"] }};  // {{ eq["label"] }}
{% endif %}
{% endif %}
{% endfor %}

{% for var in state_variables %}
    mStateVariables.push_back({{ var["id"] }});
{% endfor %}

{% for var in state_variables %}
    SetDefaultInitialCondition({{ var["index"] }}, {{ var["id"] }});
{% endfor %}

{% for param in parameters %}
    mParameters.push_back({{ param["id"] }});
{% endfor %}
}

double {{ ode_class_name }}::ProcessModelEvents(double time, const std::vector<double> &rY)
{
{% if scale_time %}
    time *= {{ time_multiplier }}; // Chaste integrates in hours; use the model's native time units
{% endif %}
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

{% for event in events %}
    //========================================
    // EVENT: {{ event["label"] }}
    //========================================
    {
        double event_dist = {{ event["distance"] }};

        // active: the raw SBML trigger condition. detected: the signed distance has reached zero,
        // the point CVODE roots on. For a >=/<= trigger detected lags active by an epsilon at the
        // boundary; for a non-relational trigger the distance is a constant >= 0 so detected is
        // always true and the logic below reduces to the raw condition.
        // The fire is gated on detected so an event whose crossing lands on a sample grid point is
        // not latched by an uncommitted evaluation at the grid point before it can be applied. The
        // satisfied/re-arm state tracks active, so a trigger that is already true at the initial
        // condition (initialValue=true) stays satisfied and does not spuriously fire.
        bool active = {{ event["trigger"] }};
        bool detected = event_dist >= 0.0;

        // Suppress an event whose trigger was already active when this Solve segment started
        // (a carried-over trigger) by forcing a large negative distance, so CVODE reports no
        // spurious root at the initial condition. mEventClampActive is frozen at segment start
        // (CalculateStoppingEvent) and cleared below the instant the trigger first goes false.
        // Using this monotonic per-segment flag rather than the live, in-step-mutated
        // mEventSatisfied keeps the root function stable across CVODE's root bracketing, so an
        // event localizes at its true crossing instead of the integration step endpoint.
        if (mEventClampActive[{{ event["index"] }}] && active)
        {
            event_dist = -(std::abs(event_dist) + 1.0);
        }

        // Update max_dist (closest event to triggering)
        if (event_dist > max_dist)
        {
            max_dist = event_dist;
        }

        // Process the event
        if (active)
        {
            if (!mEventSatisfied[{{ event["index"] }}] && detected)
            {
                // The condition is transitioning from false -> true: trigger the event. The
                // assignment values are recorded below (not here) so they can be refreshed at the
                // committed event point.
                mEventTriggered[{{ event["index"] }}] = true;
            }
            // Latch only once the distance has crossed zero. Until then (active but not yet
            // detected, i.e. within the epsilon boundary) leave the satisfied state untouched, so
            // a crossing exactly on a grid point is neither prematurely latched nor, when the
            // trigger is already true at t=0, re-armed.
            if (detected)
            {
                mEventSatisfied[{{ event["index"] }}] = true;
            }
        }
        else if (!mEventTriggered[{{ event["index"] }}])
        {
            // Trigger is false and the event has not fired in this segment, so it (re-)arms:
            // clear the satisfied latch and the clamp (the clamp permanently, monotonically,
            // so it stays stable across CVODE's in-step root bracketing and the next rising
            // edge is detected). Once the event HAS fired this segment we leave these sticky,
            // so a later root-bracketing evaluation that lands on the trigger-false side cannot
            // undo the fire and leave the event spuriously unsatisfied (which would re-fire it
            // at the next segment's initial condition).
            mEventSatisfied[{{ event["index"] }}] = false;
            mEventClampActive[{{ event["index"] }}] = false;
        }

        // Record this event's assignments while it has fired this segment - re-evaluated on every
        // call rather than only at first detection. The assignment is still deferred (applied by
        // AdjustParameters at the committed point), but recording it here lets the harness re-run
        // ProcessModelEvents at the localized root before applying, so a state-dependent value
        // (e.g. a compartment-resize rescale S * C_old / C_new) uses the root state rather than the
        // integration step where the event was first detected. event_priority orders simultaneous
        // events: an assignment only overwrites one already recorded this firing if its event has
        // lower-or-equal priority, so the lowest-priority event - which SBML executes last - wins a
        // conflict. Events with no priority use +inf, reducing to last-writer-wins.
        if (mEventTriggered[{{ event["index"] }}])
        {
{% if event["priority"] is not none %}
            [[maybe_unused]] double event_priority = {{ event["priority"] }};
{% else %}
            [[maybe_unused]] double event_priority = std::numeric_limits<double>::max();
{% endif %}
{% for assignment in event["assignments"] %}
{% if ( assignment["type"] == VarType.STATE_VARIABLE ) %}
            // {{ assignment["lhs"] }} = {{ assignment["rhs"] }}
            if (!mEventAdjustedStateVars[{{ assignment["index"] }}]
                || event_priority <= mEventAdjustedStatePriority[{{ assignment["index"] }}])
            {
                mEventAdjustedStateVars[{{ assignment["index"] }}] = true;
                mEventAdjustedStateValues[{{ assignment["index"] }}] = {{ assignment["rhs"] }};
                mEventAdjustedStatePriority[{{ assignment["index"] }}] = event_priority;
            }

{% elif ( assignment["type"] == VarType.PARAMETER ) %}
            // {{ assignment["lhs"] }} = {{ assignment["rhs"] }}
            if (!mEventAdjustedParameters[{{ assignment["index"] }}]
                || event_priority <= mEventAdjustedParameterPriority[{{ assignment["index"] }}])
            {
                mEventAdjustedParameters[{{ assignment["index"] }}] = true;
                mEventAdjustedParameterValues[{{ assignment["index"] }}] = {{ assignment["rhs"] }};
                mEventAdjustedParameterPriority[{{ assignment["index"] }}] = event_priority;
            }

{% else %}
            {{ assignment["lhs"] }} = {{ assignment["rhs"] }}; {# TODO: does this case exist? #}

{% endif %}
{% endfor %} {# 'for assignment in event["assignments"]' #}
        }
    }

{% endfor %} {# 'for event in events' #}

    return max_dist; // Signed distance of the event closest to triggering
}

std::vector<double> {{ ode_class_name }}::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
{% for var in state_variables %}
    {{ var["id"] }} = rStateVariables[{{ var["index"] }}];
{% endfor %}

{% for param in parameters %}
    {{ param["id"] }} = GetParameter({{ param['index'] }});
{% endfor %}

{% for eq in equations %}
{% if ( eq["type"] not in [EquationType.INITIAL_VALUE, EquationType.INITIAL_ASSIGNMENT, EquationType.CONVERSION] ) %}
{% if eq["local_parameters"] %}
    // {{ eq["var"] }}: {{ eq["label"] }}
    {
{% for local_parameter in eq["local_parameters"] %}
        [[maybe_unused]] double {{ local_parameter["id"] }} = {{ local_parameter["value"] }};
{% endfor %}
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->{{ eq["var"] }} = {{ eq["rhs"] }};
    }
{% else %}
    {{ eq["var"] }} = {{ eq["rhs"] }};  // {{ eq["label"] }}
{% endif %}
{% endif %}
{% endfor %}

    std::vector<double> derivatives({{ state_variables|length }});
{% for var in state_variables %}
    derivatives[{{ var["index"] }}] = {{ var["derivative_id"] }};
{% endfor %}
    return derivatives;
}

// MODEL FUNCTIONS
{% for func in functions %}
inline double {{ ode_class_name }}::{{ func["id"] }}({{ func["args"] }})
{
    return {{ func["body"] }};
}

{% endfor %}

template <>
void CellwiseOdeSystemInformation<{{ ode_class_name }}>::Initialise()
{
    // STATE VARIABLES
{% for var in state_variables %}
    this->mVariableNames.push_back("{{ var['name'] }}");
    this->mVariableUnits.push_back("{{ var['units'] }}");
{% if var["initial_value"] is not none %}
    this->mInitialConditions.push_back({{ var['initial_value'] }});
{% else %}
    this->mInitialConditions.push_back(0.0);
{% endif %}

{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    this->mDerivedQuantityNames.push_back("{{ dq['name'] }}");
    this->mDerivedQuantityUnits.push_back("{{ dq['units'] }}");

{% endfor %}

    // PARAMETERS
{% for var in parameters %}
    this->mParameterNames.push_back("{{ var['name'] }}");
    this->mParameterUnits.push_back("{{ var['units'] }}");

{% endfor %}
    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
