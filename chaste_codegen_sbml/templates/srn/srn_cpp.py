function_impl_template = """
    double {ode_name}::{fn_id}({fn_args})
    {{
        return {fn_body_cpp};
    }}
"""

state_param_template = "double {par_id} = this->mParameters[{par_num}]; // {par_name}"

srn_cpp_template = """
#include "{srn_model_name}.hpp"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
{ode_system_name}::{ode_system_name}(std::vector<double> stateVariables)
    : AbstractOdeSystem({num_state_vars})
{{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{ode_system_name}>);

    Init();

    {species_defaults}

    {parameter_defaults}

    if (stateVariables != std::vector<double>())
    {{
        SetStateVariables(stateVariables);
    }}
}}

{ode_system_name}::~{ode_system_name}()
{{
}}

{functions_impl}

void {ode_system_name}::Init()
 {{
    /* Initialise model compartments. */
    {compartment_init}

    /* Initialise model parameters. */
    {parameter_init}

    /* Initialise vector to check if events have been triggered. */
    {event_vector_init}

}}

void {ode_system_name}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{{

    /* Define state variables */
    {state_var_def}

    /* Define state parameters */
    {species_state_param_def}

    {parameter_state_param_def}

    /* Define algebraic rules. */
    {rule_def}

    /* Define the reactions in this model. */
    {reaction_def}

    {ode_def}

    /* Account for the differences in timescales. */
    {ode_timescale_def}

}}

template <>
void CellwiseOdeSystemInformation<{ode_system_name}>::Initialise()
{{
    {species_ode_init}

    {parameter_ode_init}

    this->mInitialised = true;
}}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{ode_system_name}, {num_state_vars}> {srn_model_name};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({ode_system_name})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {ode_system_name}, {num_state_vars})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({srn_model_name})
"""
