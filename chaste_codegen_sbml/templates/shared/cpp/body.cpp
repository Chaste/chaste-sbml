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
{% if sp["is_state_variable"] is true() %}
    {{ sp["id"] }} = {{ sp["concentration"] }}; // {{ sp["name"] }}
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
{% endif %}
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

{% for fd in function_definitions %}
double {{ ode_class_name }}::{{ fd["id"] }}({{ fd["args"] }})
{
    return {{ fd["body"]}};
}
{% endfor %}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    // STATE VARIABLES:
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
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

    // RULES:
{% for rule in rules %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
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
bool {{ ode_class_name }}::CalculateStoppingEvent(double time, const std::vector<double> & rY)
{
    if (time <= 0.0)
    {
        return false;
    }

    // STATE VARIABLES:
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    double {{ sp["id"] }} = rY[{{ sp["state_variable_index"] }}];
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

    bool stop = false;

    std::vector<double> dy(rY.size()); // Initialise derivatives vector
    EvaluateYDerivatives(time, rY, dy);

{% for event in events %}
    if ({{ event["trigger"] }})
    {
        if (!eventsSatisfied[{{ loop.index0 }}])
        {
            stop = true;
{% for assignment in event["assignments"] %}
            {{ assignment }}
{% endfor %}
        }
        eventsSatisfied[{{ loop.index0 }}] = true;
    }
    else
    {
        eventsSatisfied[{{ loop.index0 }}] = false;
    }
{% endfor %}

    return stop;
}
{% endif %}

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
