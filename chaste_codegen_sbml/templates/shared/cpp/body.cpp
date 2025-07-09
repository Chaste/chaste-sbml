#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "{{ model_hpp_file }}"

namespace sm = sbmlmath;

{{ ode_class_name }}::{{ ode_class_name }}(std::vector<double> stateVariables)
    : AbstractOdeSystem({{ num_state_vars }})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ ode_class_name }}>);

    // STATE VARIABLES
{% for var in state_variables %}
    {{ var["id"] }} = {{ var["initial_value"] }};
{% endfor %}

{% for var in state_variables %}
    SetDefaultInitialCondition({{ var["index"] }}, {{ var["id"] }});
{% endfor %}

    if (stateVariables.size() == {{ num_state_vars }})
    {
{% for var in state_variables %}
        {{ var["id"] }} = stateVariables[{{ var["index"] }}];
{% endfor %}
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("{{ ode_class_name }}: Expected {{ num_state_vars }} state variables, got " + std::to_string(stateVariables.size()));
    }

{% for var in state_variables %}
    mStateVariables.push_back({{ var["id"] }});
{% endfor %}

    // PARAMETERS
{% for param in variable_parameters %}
    {{ param["id"] }} = {{ param["initial_value"] }};
{% endfor %}

    ProcessRules(0.0, mStateVariables);

{% for param in variable_parameters %}
    mParameters.push_back({{ param["id"] }});
{% endfor %}

{% if events %}
    // EVENTS
    mEventsSatisfied.resize({{ events|length }}, false);
    mEventsInitialised = false;
{% endif %}
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    ProcessRules(time, rY);

{% for var in state_variables %}
    rDY[{{ var["index"] }}] = {{ var["rhs"] }}; // d[{{ var["name"] }}]/dt
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

void {{ ode_class_name }}::ProcessRules(const std::vector<double> &rY)
{
    // STATE VARIABLES
{% for var in state_variables %}
    {{ var["id"] }} = rY[{{ var["index"] }}];
{% endfor %}

    // RULES
{% for rule in rules %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endfor %}

    // PARAMETERS
{% for var in variable_parameters %}
{% if var["rhs"] %}
    SetParameter("{{ var['id'] }}", {{ var['id'] }});
{% endif %}
{% endfor %}

    // REACTIONS
{% for reaction in reactions %}
  {% if reaction["name"] %}
    // {{ reaction["name"] }}
  {% endif %}
  {% if reaction["parameters"] %}
    double {{ reaction["id"] }} = 0.0;
    {
  {% for param in reaction["parameters"] %}
        double {{ param["id"] }} = {{ param["value"] }};
  {% endfor %}
        {{ reaction["id"] }} = {{ reaction["rhs"] }};
    }
  {% else %}
    double {{ reaction["id"] }} = {{ reaction["rhs"] }};
  {% endif %}

{% endfor %}
}

{% if events %}
double {{ ode_class_name }}::ProcessEvents(double time, const std::vector<double> &rY)
{
    ProcessRules(time, rY);

    double min_dist = std::numeric_limits<double>::max();
    double event_dist = min_dist;

{% for event in events %}
    // EVENT: {{ event["name"] }}
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

            UpdateDefaultInitialConditions(rY); // TODO: Make this work for multiple events
{% for assignment in event["assignments"] %}
            {{ assignment }};
{% endfor %}
        }
        mEventsSatisfied[{{ loop.index0 }}] = true; // Flag the condition true
    }
    else
    {
        mEventsSatisfied[{{ loop.index0 }}] = false; // Flag the condition false
    }
{% endfor %}

    mEventsInitialised = true; // Flag that events have been processed at least once

    // Distance to closest event
    return min_dist;
}

double {{ ode_class_name }}::CalculateRootFunction(double time, const std::vector<double> &rY)
{
    return ProcessEvents(time, rY);
}

bool {{ ode_class_name }}::CalculateStoppingEvent(double time, const std::vector<double> &rY)
{
    return ProcessEvents(time, rY) == 0.0;
}

void {{ ode_class_name }}::UpdateDefaultInitialConditions(const std::vector<double> &rY)
{
{% for var in state_variables %}
    SetDefaultInitialCondition({{ var["state_variable_index"] }}, rY[{{ var["state_variable_index"] }}]); // {{ var["id"] }}
{% endfor %}
}

// FUNCTION DEFINITIONS:
{% for fd in function_definitions %}
double {{ ode_class_name }}::{{ fd["id"] }}({{ fd["args"] }})
{
    return {{ fd["body"]}};
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

typedef {{ wrapper_class_name }}<{{ ode_class_name }}, {{ num_state_vars }}> {{ model_class_name }};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
EXPORT_TEMPLATE_CLASS2({{ wrapper_class_name }}, {{ ode_class_name }}, {{ num_state_vars }})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ model_class_name }})
