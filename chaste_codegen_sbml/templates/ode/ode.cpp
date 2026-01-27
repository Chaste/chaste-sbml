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

    // PARAMETERS
{% for param in parameters %}
{% if param["value"] is not none %}
    {{ param["id"] }} = {{ param["value"] }};
{% endif %}
{% endfor %}

    // STATE VARIABLES
{% for var in state_variables %}
{% if var["initial_value"] is not none %}
    {{ var["id"] }} = {{ var["initial_value"] }};
{% endif %}
{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
{% if dq["initial_value"] is not none %}
    {{ dq["id"] }} = {{ dq["initial_value"] }};
{% endif %}
{% endfor %}

    // STOICHIOMETRY VARIABLES
{% for svar in stoichiometry_variables %}
{% if svar["rhs"] is not none %}
   {{ svar["id"] }} = {{ svar["rhs"] }};
{% elif svar["initial_value"] is not none %}
    {{ svar["id"] }} = {{ svar["initial_value"] }};
{% endif %}
{% endfor %}

    // INITIAL ASSIGNMENTS
    RunInitialAssignments(0.0);

{% for ia in initial_assignments %}
{% if ia["custom"] is true() %}
    {{ ia["lhs"] }} = {{ ia["rhs"] }}; // {{ ia["label"] }}
{% endif %}
{% endfor %}

    // ODE SYSTEM INFORMATION
{% for var in state_variables %}
    SetDefaultInitialCondition({{ var["index"] }}, {{ var["id"] }});
{% endfor %}

{% for var in state_variables %}
    mStateVariables.push_back({{ var["id"] }});
{% endfor %}

{% for param in parameters %}
    mParameters.push_back({{ param["id"] }});
{% endfor %}

    // REACTIONS
    RunReactions(0.0);

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

    RunModelRules(time, rY);

{% for dq in derived_quantities %}
    dqs.push_back({{ dq["id"] }});
{% endfor %}

    // AMOUNTS
{% for amount in amounts %}
    double {{ amount["id"] }} = {{ amount["rhs"] }};
{% endfor %}

{% for amount in amounts %}
    dqs.push_back({{ amount["id"] }});
{% endfor %}
    return dqs;
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RunModelRules(time, rY);

{% for var in state_variables %}
    rDY[{{ var["index"] }}] = {{ var["rhs"] }}; // d[{{ var["id"] }}]/dt
{% endfor %}

    // TODO: Scale time appropriately
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
{% for rule in assignment_rules %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endfor %}
}

// INITIAL ASSIGNMENTS
void {{ ode_class_name }}::RunInitialAssignments(double time)
{
{% for ia in initial_assignments %}
{% if ia["custom"] is false() %}
    {{ ia["lhs"] }} = {{ ia["rhs"] }}; // {{ ia["label"] }}
{% endif %}
{% endfor %}
}

// REACTIONS
void {{ ode_class_name }}::RunReactions(double time)
{
{% for reaction in reactions %}
  {% if reaction["label"] %}
    // {{ reaction["label"] }}
  {% endif %}
  {% if reaction["parameters"] %}
    {{ reaction["id"] }} = 0.0;
    {
  {% for param in reaction["parameters"] %}
        double {{ param["id"] }} = {{ param["value"] }};
  {% endfor %}
        {{ reaction["id"] }} = {{ reaction["rhs"] }};
    }
  {% else %}
    {{ reaction["id"] }} = {{ reaction["rhs"] }};
  {% endif %}

{% endfor %}
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

{% for amount in amounts %}
    this->mDerivedQuantityNames.push_back("{{ amount['id'] }}");
    this->mDerivedQuantityUnits.push_back("{{ amount['units'] }}");

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
