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

    mEventSatisfied.resize({{ events|length }}, true); // Prevent events from triggering at the start
    mEventTriggered.resize({{ events|length }}, false);

    mEventAdjustedParameters.resize({{ parameters|length }}, false);
    mEventAdjustedParameterValues.resize({{ parameters|length }}, 0.0);

    mEventAdjustedStateVars.resize({{ state_variables|length }}, false);
    mEventAdjustedStateValues.resize({{ state_variables|length }}, 0.0);
{% endif %} {# 'if events' #}
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

std::vector<double> {{ ode_class_name }}::ComputeDerivedQuantities(double time, const std::vector<double> &rY)
{
    std::vector<double> dqs;
{% if derived_quantities %}
    RunModelEquations(time, rY);

    // AMOUNTS
{% for eq in equations %}
{% if ( eq["type"] == EquationType.AMOUNT ) %}
    double {{ eq["lhs"] }} = {{ eq["rhs"] }}; // {{ eq["label"] }}
{% endif %}
{% endfor %}

{% for dq in derived_quantities %}
    dqs.push_back({{ dq["id"] }});
{% endfor %}
{% endif %} {# 'if derived_quantities' #}

    return dqs;
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void {{ ode_class_name }}::Initialise(double time)
{
{% for eq in equations %}
{% if ( eq["type"] != EquationType.AMOUNT ) %}
    {{ eq["lhs"] }} = {{ eq["rhs"] }}; // {{ eq["label"] }}
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
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

{% for event in events %}
    //========================================
    // EVENT: {{ event["label"] }}
    //========================================
    {
        double event_dist = {{ event["distance"] }};

        // Avoid oscillation by ensuring event_dist is not close to 0 unless triggered
        if (std::abs(event_dist) < 1.0)
        {
            event_dist = 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if ({{ event["trigger"] }})
        {
            if (!mEventSatisfied[{{ event["index"] }}])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[{{ event["index"] }}] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
{% for assignment in event["assignments"] %}
{% if ( assignment["type"] == VarType.STATE_VARIABLE ) %}
                // {{ assignment["lhs"] }} = {{ assignment["rhs"] }}
                mEventAdjustedStateVars[{{ assignment["index"] }}] = true;
                mEventAdjustedStateValues[{{ assignment["index"] }}] = {{ assignment["rhs"] }};

{% elif ( assignment["type"] == VarType.VARIABLE_PARAMETER ) %}
                // {{ assignment["lhs"] }} = {{ assignment["rhs"] }}
                mEventAdjustedParameters[{{ assignment["index"] }}] = true;
                mEventAdjustedParameterValues[{{ assignment["index"] }}] = {{ assignment["rhs"] }};

{% else %}
                {{ assignment["lhs"] }} = {{ assignment["rhs"] }}; {# TODO: does this case exist? #}

{% endif %}
{% endfor %} {# 'for assignment in event["assignments"]' #}

            }
            mEventSatisfied[{{ event["index"] }}] = true;
        }
        else
        {
            mEventSatisfied[{{ event["index"] }}] = false;
            mEventTriggered[{{ event["index"] }}] = false;
        }
    }

{% endfor %} {# 'for event in events' #}

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void {{ ode_class_name }}::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void {{ ode_class_name }}::RunInitialAssignments(double time)
{
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
{% if ( eq["type"] not in [EquationType.INITIAL_VALUE, EquationType.INITIAL_ASSIGNMENT, EquationType.AMOUNT] ) %}
{% if eq["local_parameters"] %}
    // {{ eq["lhs"] }}: {{ eq["label"] }}
    {
{% for local_parameter in eq["local_parameters"] %}
        double {{ local_parameter["id"] }} = {{ local_parameter["value"] }};
{% endfor %}
        {{ eq["lhs"] }} = {{ eq["rhs"] }};
    }
{% else %}
    {{ eq["lhs"] }} = {{ eq["rhs"] }};  // {{ eq["label"] }}
{% endif %}
{% endif %}
{% endfor %}

    std::vector<double> derivatives({{ state_variables|length }});
{% for var in state_variables %}
    derivatives[{{ var["index"] }}] = {{ var["derivative_id"] }};
{% endfor %}
    return derivatives;
}

// REACTIONS
void {{ ode_class_name }}::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void {{ ode_class_name }}::UpdateParameters(double time)
{
{% for param in parameters %}
{% if param["fixed"] is false() %}
    {{ param["id"] }} = GetParameter({{ param['index'] }});
{% endif %}
{% endfor %}
}

// STATE VARIABLES
void {{ ode_class_name }}::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
{% for var in state_variables %}
    {{ var["id"] }} = rStateVariables[{{ var["index"] }}];
{% endfor %}
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
    this->mVariableNames.push_back("{{ var['id'] }}");
    this->mVariableUnits.push_back("{{ var['units'] }}");
{% if var["initial_value"] is not none %}
    this->mInitialConditions.push_back({{ var['initial_value'] }});
{% else %}
    this->mInitialConditions.push_back(0.0);
{% endif %}

{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    this->mDerivedQuantityNames.push_back("{{ dq['id'] }}");
    this->mDerivedQuantityUnits.push_back("{{ dq['units'] }}");

{% endfor %}

    // PARAMETERS
{% for var in parameters %}
    this->mParameterNames.push_back("{{ var['id'] }}");
    this->mParameterUnits.push_back("{{ var['units'] }}");

{% endfor %}
    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
