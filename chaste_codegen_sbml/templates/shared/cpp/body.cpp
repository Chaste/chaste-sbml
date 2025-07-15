#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "{{ model_hpp_file }}"

namespace sm = sbmlmath;

{{ ode_class_name }}::{{ ode_class_name }}(std::vector<double> stateVariables)
    : AbstractOdeSystem({{ state_variables|length }})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ ode_class_name }}>);

    // STATE VARIABLES
{% for var in state_variables %}
    {{ var["id"] }} = {{ var["initial_value"] }};
{% endfor %}

{% for var in state_variables %}
    SetDefaultInitialCondition({{ var["index"] }}, {{ var["id"] }});
{% endfor %}

    if (stateVariables.size() == {{ state_variables|length }})
    {
{% for var in state_variables %}
        {{ var["id"] }} = stateVariables[{{ var["index"] }}];
{% endfor %}
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("{{ ode_class_name }}: Expected {{ state_variables|length }} state variables, got " + std::to_string(stateVariables.size()));
    }

{% for var in state_variables %}
    mStateVariables.push_back({{ var["id"] }});
{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    {{ dq["id"] }} = 0.0;
{% endfor %}

    // VARIABLE PARAMETERS
{% for param in variable_parameters %}
    {{ param["id"] }} = {{ param["initial_value"] }};
{% endfor %}

{% for param in variable_parameters %}
    mParameters.push_back({{ param["id"] }});
{% endfor %}

    // RULE-BASED PARAMETERS
{% for var in rule_based_parameters %}
    {{ var['id'] }} = 0.0;
{% endfor %}

    // REACTIONS
{% for reaction in reactions %}
    {{ reaction["id"] }} = 0.0;
{% endfor %}

{% if events %}
    // EVENTS
    mEventsSatisfied.resize({{ events|length }}, false);
    mEventsInitialised = false;
{% endif %}

    ProcessRules(0.0, mStateVariables);
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    ProcessRules(time, rY);

{% for var in state_variables %}
    rDY[{{ var["index"] }}] = {{ var["rhs"] }}; // d[{{ var["id"] }}]/dt
{% endfor %}

    // Scale time appropriately
}

{% if derived_quantities %}
std::vector<double> {{ ode_class_name }}::ComputeDerivedQuantities(double time, const std::vector<double> &rY)
{
    ProcessRules(time, rY);

    std::vector<double> dqs;
{% for dq in derived_quantities %}
    dqs.push_back({{ dq["id"] }});
{% endfor %}
    return dqs;
}
{% endif %}

void {{ ode_class_name }}::ProcessRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
{% for var in state_variables %}
    {{ var["id"] }} = rY[{{ var["index"] }}];
{% endfor %}

    // VARIABLE PARAMETERS
{% for var in variable_parameters %}
    {{ var["id"] }} = GetParameter({{ var['index'] }});
{% endfor %}

    // ASSIGNMENT RULES
{% for rule in assignment_rules %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endfor %}

    // REACTIONS
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

{% if events %}
double {{ ode_class_name }}::ProcessEvents(double time, const std::vector<double> &rY)
{
    ProcessRules(time, rY);

    double min_dist = std::numeric_limits<double>::max();
    double event_dist = min_dist;

    std::vector<bool> state_vars_updated({{ state_variables|length }}, false);

{% for event in events %}
    // EVENT: {{ event["label"] }}
    event_dist = {{ event["distance"] }};

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
        if (!mEventsSatisfied[{{ loop.index0 }}] && mEventsInitialised)
        {
            // The condition is transitioning from false to true,
            // and this is not the first time-step => trigger the event.
            event_dist = 0.0;
            min_dist = 0.0;

{% for assignment in event["assignments"] %}
{% if ( assignment["type"] == VarType.STATE_VARIABLE ) %}
            SetStateVariable({{ assignment["index"] }}, {{ assignment["rhs"] }});
            state_vars_updated[{{ assignment["index"] }}] = true;
{% elif ( assignment["type"] == VarType.VARIABLE_PARAMETER ) %}
            SetParameter({{ assignment["index"] }}, {{ assignment["rhs"] }});
{% else %}
            {{ assignment["lhs"] }} = {{ assignment["rhs"] }}; 
{% endif %}
{% endfor %}
        }
        mEventsSatisfied[{{ loop.index0 }}] = true; // Flag the condition true
    }
    else
    {
        mEventsSatisfied[{{ loop.index0 }}] = false; // Flag the condition false
    }
{% endfor %}

    // Event triggered, update state variables if necessary
    if (min_dist == 0.0)
    {
        for (unsigned i = 0; i < {{ state_variables|length }}; ++i)
        {
            if (!state_vars_updated[i])
            {
                SetStateVariable(i, rY[i]);
            }
        }
    }

    mEventsInitialised = true; // Flag that events have been processed at least once

    // Distance to closest event
    return min_dist;
}
{% endif %}

{% if events %}
double {{ ode_class_name }}::CalculateRootFunction(double time, const std::vector<double> &rY)
{
    return ProcessEvents(time, rY);
}
{% endif %}

{% if events %}
bool {{ ode_class_name }}::CalculateStoppingEvent(double time, const std::vector<double> &rY)
{
    return ProcessEvents(time, rY) == 0.0;
}
{% endif %}

// FUNCTIONS
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
    this->mInitialConditions.push_back({{ var['initial_value'] }});

{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    this->mDerivedQuantityNames.push_back("{{ dq['id'] }}");
    this->mDerivedQuantityUnits.push_back("{{ dq['units'] }}");

{% endfor %}

    // PARAMETERS
{% for var in variable_parameters %}
    this->mParameterNames.push_back("{{ var['id'] }}");
    this->mParameterUnits.push_back("{{ var['units'] }}");

{% endfor %}
    this->mInitialised = true;
}

// Define {{ wrapper_class_name }} using wrappers
#include "{{ wrapper_class_name }}.hpp"
#include "{{ wrapper_class_name }}.cpp"

typedef {{ wrapper_class_name }}<{{ ode_class_name }}, {{ state_variables|length }}> {{ model_class_name }};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
EXPORT_TEMPLATE_CLASS2({{ wrapper_class_name }}, {{ ode_class_name }}, {{ state_variables|length }})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ model_class_name }})
