#include "{{model_hpp_file}}"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
{{ode_class_name}}::{{ode_class_name}}(std::vector<double> stateVariables)
    : AbstractOdeSystem({{num_state_vars}})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ode_class_name}}>);

    Init();

    {{species_defaults}}

{% for parameter in parameters %}
{% if parameter["is_defined"] is false() %}
    this->mParameters.push_back({{parameter["default"]}}); // {{parameter["varname"]}}
{% endif %}
{% endfor %}

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

{{ode_class_name}}::~{{ode_class_name}}()
{
}

{{functions_impl}}

void {{ode_class_name}}::Init()
{
{% if compartments %}
    /* Initialise model compartments. */
{% for compartment in compartments %}
    {{compartment["id"]}} = {{compartment["size"]}}; // {{compartment["varname"]}}
{% endfor %}
{% endif %}

{% if parameters %}
    /* Initialise model parameters. */
{% for parameter in parameters %}
    {{parameter["id"]}} = {{parameter["value"]}};
{% endfor %}
{% endif %}

    /* Initialise vector to check if events have been triggered. */
    {{event_vector_init}}

}

void {{ode_class_name}}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{

    /* Define state variables */
    {{state_var_def}}

    /* Define state parameters */
    {{species_state_param_def}}

{% if parameters %}
    /* Define states: parameters. */
{% for parameter in parameters %}
{% if parameter["is_state"] is true() %}
    double {{parameter["id"]}} = this->mParameters[{{parameter["state_index"]}}]; // {{parameter["name"]}}
{% endif %}
{% endfor %}
{% endif %}

{% if rules %}
    /* Define algebraic rules. */
{% for rule in rules %}
    {{rule["id"]}} = {{rule["formula"]}};
{% endfor %}
{% endif %}

    /* Define the reactions in this model. */
    {{reaction_def}}

    {{ode_def}}

    /* Account for the differences in timescales. */
    {{ode_timescale_def}}

}

{% if num_events > 0 %}
bool {{ode_class_name}}::CalculateStoppingEvent(double time, const std::vector<double> & rY)
{
    // Return true if all events have been triggered.
    return AreAllEventsSatisfied(time, rY);
}

void {{ode_class_name}}::CheckAndUpdateEvents(double time, const std::vector<double> & rY)
{
    std::vector<double> dy(rY.size()); // Initialise derivatives vector
    EvaluateYDerivatives(time, rY, dy);
    {{event_defs}}
}

bool {{ode_class_name}}::AreAllEventsSatisfied(double time, const std::vector<double>& rY)
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
void CellwiseOdeSystemInformation<{{ode_class_name}}>::Initialise()
{
    {{species_ode_init}}

{% if parameters %}
    /* Define states: parameters. */
    // Parameters without set values must be externally defined
{% for parameter in parameters %}
{% if parameter["is_defined"] is false() %}
    this->mParameterNames.push_back("{{parameter['varname']}}");
    this->mParameterUnits.push_back("{{parameter['units']}}");

{% endif %}
{% endfor %}
{% endif %}

    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "{{wrapper_class_name}}.hpp"
#include "{{wrapper_class_name}}.cpp"

typedef {{wrapper_class_name}}<{{ode_class_name}}, {{num_state_vars}}> {{model_class_name}};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ode_class_name}})
EXPORT_TEMPLATE_CLASS2({{wrapper_class_name}}, {{ode_class_name}}, {{num_state_vars}})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{model_class_name}})
