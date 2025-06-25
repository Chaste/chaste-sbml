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

{% if compartments %}
    // COMPARTMENTS:
{% for comp in compartments %}
    {{ comp["id"] }} = {{ comp["size"] }};
{% endfor %}
{% endif %}

    // STATE VARIABLES:
{% for sp in species %}
{% if sp["has_ode"] is true() %}
    {{ sp["id"] }} = {{ sp["concentration"] }}; // {{ sp["name"] }}
{% endif %}
{% endfor %}

{% for rule in rules %}
{% if rule["is_state_variable"] is true() %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endif %}
{% endfor %}

{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    SetDefaultInitialCondition({{ sp["state_variable_index"] }}, {{ sp["id"] }});
{% endif %}
{% endfor %}

    if (stateVariables.size() == {{ num_state_vars }})
    {
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
        {{ sp["id"] }} = stateVariables[{{ sp["state_variable_index"] }}];
{% endif %}
{% endfor %}
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("{{ ode_class_name }}: Expected {{ num_state_vars }} state variables, got " + std::to_string(stateVariables.size()));
    }

{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    mStateVariables.push_back({{ sp["id"] }});
{% endif %}
{% endfor %}

    // STATE PARAMETERS:
{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    {{ sp["id"] }} = {{ sp["concentration"] }};
{% endif %}
{% endfor %}

{% for parameter in parameters %}
{% if parameter["is_state_parameter"] is true() %}
    {{ parameter["id"] }} = {{ parameter["value"] }};
{% endif %}
{% endfor %}

{% for rule in rules %}
{% if rule["is_state_variable"] is false() %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endif %}
{% endfor %}

{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    mParameters.push_back({{ sp["id"] }});
{% endif %}
{% endfor %}

{% for param in parameters %}
{% if param["is_state_parameter"] is true() %}
    mParameters.push_back({{ param["id"] }});
{% endif %}
{% endfor %}

{% if events %}
    // EVENTS:
    eventsSatisfied.resize({{ events|length }}, false);
    eventsInitialised = false;
{% endif %}
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

void {{ ode_class_name }}::RefreshState(const std::vector<double> &rY)
{
    // STATE VARIABLES:
{% for sp in species %}
{% if sp["has_ode"] is true() %}
    {{ sp["id"] }} = rY[{{ sp["state_variable_index"] }}];
{% endif %}
{% endfor %}

    // STATE PARAMETERS:
{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    {{ sp["id"] }} = GetParameter("{{ sp['id'] }}");
{% endif %}
{% endfor %}

{% for param in parameters %}
{% if param["is_state_parameter"] is true() %}
    {{ param["id"] }} = GetParameter("{{ param['id'] }}");
{% endif %}
{% endfor %}
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RefreshState(rY);

    // RULES:
{% for rule in rules %}
{% if rule["is_state_variable"] is true() %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endif %}
{% endfor %}

{% for rule in rules %}
{% if rule["is_state_variable"] is false() %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endif %}
{% endfor %}

    // UPDATE STATE PARAMETERS:
{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    SetParameter("{{ sp['id'] }}", {{ sp['id'] }});
{% endif %}
{% endfor %}

{% for param in parameters %}
{% if param["is_state_parameter"] is true() %}
    SetParameter("{{ param['id'] }}", {{ param['id'] }});
{% endif %}
{% endfor %}

{% if reactions %}
    // REACTIONS:

{% for reaction in reactions %}
  {% if reaction["name"] %}
    // {{ reaction["name"] }}
  {% endif %}
  {% if reaction["parameters"] %}
    double {{ reaction["varname"] }} = 0.0;
    {
  {% for param in reaction["parameters"] %}
        double {{ param["varname"] }} = {{ param["value"] }};
  {% endfor %}
        {{ reaction["varname"] }} = {{ reaction["rhs"] }};
    }
  {% else %}
    double {{ reaction["varname"] }} = {{ reaction["rhs"] }};
  {% endif %}

{% endfor %}
{% endif %}
    // ODES:
{% for sp in species %}
{% if sp["ode"] %}
    {{ sp["ode"]["lhs"] }} = {{ sp["ode"]["rhs"] }}; // d[{{ sp["name"] }}]/dt
{% endif %}
{% endfor %}

    // Scale time appropriately
}

{% if events %}
double {{ ode_class_name }}::ProcessEvents(double time, const std::vector<double> &rY)
{
    RefreshState(rY);

    double min_dist = std::numeric_limits<double>::max();
    double event_dist = min_dist;

{% for event in events %}
    // EVENT: {{ event["trigger"] }}
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
        if (!eventsSatisfied[{{ loop.index0 }}] && eventsInitialised)
        {
            // The condition is transitioning from false to true,
            // and this is not the first time-step => trigger the event.
            event_dist = 0.0;
            min_dist = 0.0;

            UpdateDefaultInitialConditions(rY);
{% for assignment in event["assignments"] %}
            {{ assignment }};
{% endfor %}
        }
        eventsSatisfied[{{ loop.index0 }}] = true; // Flag the condition true
    }
    else
    {
        eventsSatisfied[{{ loop.index0 }}] = false; // Flag the condition false
    }
{% endfor %}

    eventsInitialised = true; // Flag that events have been processed at least once

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
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    SetDefaultInitialCondition({{ sp["state_variable_index"] }}, rY[{{ sp["state_variable_index"] }}]); // {{ sp["id"] }}
{% endif %}
{% endfor %}
}
{% endif %}

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
    // STATE VARIABLES:
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    this->mVariableNames.push_back("{{ sp['id'] }}");
    this->mVariableUnits.push_back("{{ sp['units'] }}");
    this->mInitialConditions.push_back({{ sp['concentration'] }});

{% endif %}
{% endfor %}

    // STATE PARAMETERS:
{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    this->mParameterNames.push_back("{{ sp['id'] }}");
    this->mParameterUnits.push_back("{{ sp['units'] }}");

{% endif %}
{% endfor %}
{% for param in parameters %}
{% if param["is_state_parameter"] is true() %}
    this->mParameterNames.push_back("{{ param['id'] }}");
    this->mParameterUnits.push_back("{{ param['units'] }}");

{% endif %}
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
