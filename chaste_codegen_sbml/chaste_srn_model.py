from .chaste_model import ODE_SUFFIX, TAB, ChasteModel
from .templates.srn.srn_cpp import function_impl_template
from .templates.srn.srn_hpp import srn_hpp_template

HEADER_GUARD_SUFFIX = "ODESYSTEMANDSRNMODEL_HPP_"
SRN_SUFFIX = "SrnModel"


class ChasteSRNModel(ChasteModel):

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        self.ode_name = self.model_name + ODE_SUFFIX
        self.srn_name = self.model_name + SRN_SUFFIX

        self.size = self.model.getNumReactions() + self.model.getNumRules()

    def generate_hpp(self) -> str:
        """Generate the Chaste header for an SRN model from SBML data.

        return: The generated header file as a string.
        """

        # Get inputs for the header file template
        header_guard_str = self.model_name.upper() + self.HEADER_GUARD_SUFFIX

        compartments = self.model.getListOfCompartments()
        compartment_decls_str = f"\n{TAB}".join(
            [f"double {self.get_identifier(c)};" for c in compartments]
        )

        parameters = self.model.getListOfParameters()
        parameter_decls_str = f"\n{TAB}".join(
            [f"double {self.get_identifier(p)};" for p in parameters]
        )

        functions = self.model.getListOfFunctionDefinitions()
        function_decls_list = []
        for fn in functions:
            fn_id = fn.getId()
            fn_args = self.get_function_definition_arguments(fn)
            function_decls_list.append(f"double {fn_id}({', '.join(fn_args)});")
        function_decls_str = f"\n{TAB}".join(function_decls_list)

        # Apply inputs to the header file template
        hpp = srn_hpp_template.format(
            header_guard_str=header_guard_str,
            ode_name=self.ode_name,
            srn_name=self.srn_name,
            model_name=self.model_name,
            size=self.size,
            compartment_decls_str=compartment_decls_str,
            parameter_decls_str=parameter_decls_str,
            function_decls_str=function_decls_str,
        )

        return hpp

    def generate_cpp(self) -> str:
        """Generate the Chaste source for an SRN model from SBML data.

        return: The generated source file as a string.
        """

        species = self.model.getListOfSpecies()
        species_inits = []
        for i, s in enumerate(species):
            conc = self.get_species_concentration(s)
            species_inits.append(f"SetDefaultInitialCondition({i}, {conc});")
        species_init_str = f"\n{TAB}".join(species_inits)

        compartments = self.model.getListOfCompartments()
        compartment_inits = []
        for compartment in compartments:
            c_id = self.get_identifier(compartment)
            size = compartment.getSize()
            compartment_inits.append(f"{c_id} = {size};")
        compartment_init_str = f"\n{TAB}".join(compartment_inits)

        parameters = self.model.getListOfParameters()
        parameter_vector_inits = []
        parameter_inits = []
        for parameter in parameters:
            p_id = self.get_identifier(parameter)
            value = self.get_parameter_value(parameter)
            parameter_vector_inits.append(f"this->mParameters.push_back({value});")
            parameter_inits.append(f"{p_id} = {value};")
        parameter_vector_init_str = f"\n{TAB}".join(parameter_vector_inits)
        parameter_init_str = f"\n{TAB}".join(parameter_inits)

        functions = self.model.getListOfFunctionDefinitions()
        function_impls = []
        for fn in functions:
            fn_id = fn.getId()
            args_list = self.get_function_definition_arguments(fn)
            body_cpp = self.convert_function_body(fn.getBody())
            impl = function_impl_template.format(
                ode_name=self.ode_name,
                fn=fn_id,
                fn_args=", ".join(args_list),
                fn_body_cpp=body_cpp,
            )
            function_impls.append(impl)
        functions_impl_str = "\n".join(function_impls)

        num_events = self.model.getNumEvents()
        event_vector_init_str = ""
        if num_events > 0:
            event_vector_init_str = f"eventsSatisfied.resize({num_events}, false);"

        cpp = f"""
#include "{self.srn_name}.hpp"
#include "CellwiseOdeSystemInformation.hpp"

// SBML ODE System
{self.ode_name}::{self.ode_name} (std::vector<double> stateVariables)
    : AbstractOdeSystem({self.num_species})
{{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{self.ode_name}>);

    Init();

    {species_init_str}

    {parameter_vector_init_str}

    if (stateVariables != std::vector<double>())
    {{
        SetStateVariables(stateVariables);
    }}
}}

{self.ode_name}::~{self.ode_name}()
{{
}}

{functions_impl_str}

void {self.ode_name}::Init()
 {{
    // Initialise the compartments.
    {compartment_init_str}

    // Initialise the parameters.
    {parameter_init_str}

    // Initialise vector to check if events have been triggered.
    {event_vector_init_str}
    
}}

void {self.ode_name}::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{{
    
    /* Define state variables */
    double C = rY[0]; // Cyclin
    double M = rY[1]; // cdc_2_kinase
    double X = rY[2]; // Cyclin_Protease

    /* Define state parameters */

     /* Define algebraic rules. */
    V1 = C * VM1 * pow(C + Kc, -1);
    V3 = M * VM3;

    /* Define the reactions in this model. */
    // creation of cyclin
    double vi = 0.025;
    double reaction1 = cell * vi;

    // default degradation of cyclin
    double kd = 0.01;
    double reaction2 = C * cell * kd;

    // cdc2 kinase triggered degration of cyclin
    double vd = 0.25;
    double Kd = 0.02;
    double reaction3 = C * cell * vd * X * pow(C + Kd, -1);

    // activation of cdc2 kinase
    double K1 = 0.005;
    double reaction4 = cell * (1 + -1 * M) * V1 * pow(K1 + -1 * M + 1, -1);

    // deactivation of cdc2 kinase
    double V2 = 1.5;
    double K2 = 0.005;
    double reaction5 = cell * M * V2 * pow(K2 + M, -1);

    // activation of cyclin protease
    double K3 = 0.005;
    double reaction6 = cell * V3 * (1 + -1 * X) * pow(K3 + -1 * X + 1, -1);

    // deactivation of cyclin protease
    double K4 = 0.005;
    double V4 = 0.5;
    double reaction7 = cell * V4 * X * pow(K4 + X, -1);


    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // dCyclin/dt
    rDY[1] = (reaction4 - reaction5) / cell; // dcdc_2_kinase/dt
    rDY[2] = (reaction6 - reaction7) / cell; // dCyclin Protease/dt

    /* Account for the differences in timescales. */
    // rDY[0] *= 3600.0;
    // rDY[1] *= 3600.0;
    // rDY[2] *= 3600.0;

    // std::cout << rDY[0] << ", " << rDY[1] << ", " << rDY[2] << ", " << std::endl;
}}

template<>
void CellwiseOdeSystemInformation<{self.ode_name}>::Initialise()
{{
    this->mVariableNames.push_back("Cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("cdc_2_kinase");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("Cyclin Protease");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    this->mInitialised = true;
}}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{self.ode_name}, {self.num_species}> {self.srn_name};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({self.ode_name})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {self.ode_name}, {self.num_species})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({self.srn_name})
"""

        # Get the derivative function
        derivs_fn = GetEvaluateYDerivativesVoidString(filename, model)
        source_file.write(derivs_fn)

        # Initialise function
        initialise_fn = GetInitialiseString(filename, model)
        source_file.write(initialise_fn)

        # Define SRN Model
        srn_model_defn = GetModelDefinitionString(filename, model, False)
        source_file.write(srn_model_defn)

        return cpp

    def generate_code(self) -> None:
        """Generate the Chaste code for an SRN model from SBML data."""
        self.hpp_source = self.generate_hpp()
        self.cpp_source = self.generate_cpp()
