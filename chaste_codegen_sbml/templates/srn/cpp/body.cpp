#include "{{model_hpp_file}}"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
{{ode_class_name}}::{{ode_class_name}}(std::vector<double> stateVariables)
    : AbstractOdeSystem({{num_state_vars}})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ode_class_name}}>);

    Init();

    {{species_defaults}}

    {{parameter_defaults}}

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

    /* Initialise model parameters. */
    {{parameter_init}}

    /* Initialise vector to check if events have been triggered. */
    {{event_vector_init}}

}

void {{ode_class_name}}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{

    /* Define state variables */
    {{state_var_def}}

    /* Define state parameters */
    {{species_state_param_def}}

    {{parameter_state_param_def}}

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

template <>
void CellwiseOdeSystemInformation<{{ode_class_name}}>::Initialise()
{
    {{species_ode_init}}

    {{parameter_ode_init}}

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
