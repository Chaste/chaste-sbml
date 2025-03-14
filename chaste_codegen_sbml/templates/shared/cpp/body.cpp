#include "{{ model_hpp_file }}"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
{{ ode_class_name }}::{{ ode_class_name }}(std::vector<double> stateVariables)
    : AbstractOdeSystem({{ num_state_vars }})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ ode_class_name }}>);

    Init();

{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    SetDefaultInitialCondition({{ sp["state_variable_index"] }}, {{ sp["concentration"] }}); // {{ sp["name"] }}"
{% endif %}
{% endfor %}

{% for param in parameters %}
{% if param["is_defined"] is false() %}
    this->mParameters.push_back({{ param["default"] }}); // {{ param["varname"] }}
{% endif %}
{% endfor %}

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

{% for fd in function_definitions %}
{% if fd["args"] %}{% set args = "double " ~ fd["args"] | join(", double ") %}{% endif %}
double {{ ode_class_name }}::{{ fd["id"] }}({{ args }})
{
    return {{ fd["body"]}};
}
{% endfor %}

{{ functions_impl }}

void {{ ode_class_name }}::Init()
{
{% if compartments %}
    /* Initialise model compartments. */
{% for comp in compartments %}
    {{ comp["id"] }} = {{ comp["size"] }}; // {{ comp["varname"] }}
{% endfor %}
{% endif %}

{% if parameters %}
    /* Initialise model parameters. */
{% for param in parameters %}
    {{ param["id"] }} = {{ param["value"] }};
{% endfor %}
{% endif %}

{% if events %}
    /* Initialise vector to check if events have been triggered. */
    eventsSatisfied.resize({{ events|length }}, false);
{% endif %}
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
{% if species %}
    /* Define state variables */
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    double {{ sp["id"] }} = rY[{{ sp["state_variable_index"] }}]; // {{ sp["name"] }}
{% endif %}
{% endfor %}
{% endif %}

{% if species %}
    /* Define state parameters (Species). */
{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    double {{ sp["id"] }} = this->mParameters[{{ sp["state_parameter_index"] }}]; // {{ sp["name"] }}
{% endif %}
{% endfor %}
{% endif %}

{% if parameters %}
    /* Define state parameters (Parameters). */
{% for param in parameters %}
{% if param["is_state_parameter"] is true() %}
    double {{ param["id"] }} = this->mParameters[{{ param["state_parameter_index"] }}]; // {{ param["name"] }}
{% endif %}
{% endfor %}
{% endif %}

{% if rules %}
    /* Define algebraic rules. */
{% for rule in rules %}
    {{ rule["id"] }} = {{ rule["formula"] }};
{% endfor %}
{% endif %}

{% if reactions %}
    /* Define the reactions in this model. */
{% for reaction in reactions %}
{% if reaction["name"] %}
    // {{ reaction["name"] }}
{% endif %}
{% for param in reaction["parameters"] %}
    double {{ param["varname"] }} = {{ param["value"] }}; // {{ param["name"] }}
{% endfor %}
    double {{ reaction["varname"] }} = {{ reaction["rhs"] }};

{% endfor %}
{% endif %}

{% for sp in species %}
{% if sp["ode"] %}
    {{ sp["ode"]["lhs"] }} = {{ sp["ode"]["rhs"] }}; // d{{ sp["varname"] }}/dt
{% endif %}
{% endfor %}

    /* Account for the differences in timescales. */
    {{ ode_timescale_def }}

}

{% if events %}
bool {{ ode_class_name }}::CalculateStoppingEvent(double time, const std::vector<double> & rY)
{
    // Return true if all events have been triggered.
    return AreAllEventsSatisfied(time, rY);
}

void {{ ode_class_name }}::CheckAndUpdateEvents(double time, const std::vector<double> & rY)
{
    std::vector<double> dy(rY.size()); // Initialise derivatives vector
    EvaluateYDerivatives(time, rY, dy);
{% for event in events %}
    if ({{ event["trigger"] }})
    {
{% for assignment in event["assignments"] %}
        {{ assignment }}
{% endfor %}
        eventsSatisfied[{{ loop.index0 }}] = true;
    }
{% endfor %}
}

bool {{ ode_class_name }}::AreAllEventsSatisfied(double time, const std::vector<double>& rY)
{
    CheckAndUpdateEvents(time, rY);
    bool events_satisfied = true;
    if (std::find(eventsSatisfied.begin(), eventsSatisfied.end(), false) != eventsSatisfied.end())
    {
        events_satisfied = false;
    }
    if (events_satisfied) // Reset events vector if cell division is triggered
    {
        std::fill(eventsSatisfied.begin(), eventsSatisfied.end(), false);
    }
    return events_satisfied;
}
{% endif %}

template <>
void CellwiseOdeSystemInformation<{{ ode_class_name }}>::Initialise()
{
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    this->mVariableNames.push_back("{{ sp['name'] }}");
    this->mVariableUnits.push_back("{{ sp['units'] }}");
    this->mInitialConditions.push_back({{ sp['concentration'] }});

{% elif sp["is_state_parameter"] is true() %}
    this->mParameterNames.push_back("{{ sp['name'] }}");
    this->mParameterUnits.push_back("{{ sp['units'] }}");

{% endif %}
{% endfor %}

{% if parameters %}
    /* Define states: parameters. */
    // Parameters without set values must be externally defined
{% for param in parameters %}
{% if param["is_defined"] is false() %}
    this->mParameterNames.push_back("{{ param['varname'] }}");
    this->mParameterUnits.push_back("{{ param['units'] }}");

{% endif %}
{% endfor %}
{% endif %}
    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "{{ wrapper_class_name }}.hpp"
#include "{{ wrapper_class_name }}.cpp"

typedef {{ wrapper_class_name }}<{{ ode_class_name }}, {{ num_state_vars }}> {{ model_class_name }};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
EXPORT_TEMPLATE_CLASS2({{ wrapper_class_name }}, {{ ode_class_name }}, {{ num_state_vars }})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ model_class_name }})
