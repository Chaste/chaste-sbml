"""Module for the ChasteSbmlSrnModel class."""

import pathlib

from ._config import ODE_SUFFIX, SRN_SUFFIX, TAB
from ._utils import (
    convert_formula,
    convert_function_body,
    get_function_definition_arguments,
    get_index_by_id,
    get_species_concentration,
)
from .chaste_sbml_model import ChasteSbmlModel


class ChasteSbmlSrnModel(ChasteSbmlModel):
    """Class for generating Chaste code for an SRN model from SBML data."""

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml_file: str, model_name: str = None, **kwargs) -> None:
        """Initialise the ChasteSbmlSrnModel."""
        super().__init__(sbml_file, model_name, **kwargs)

    @property
    def srn_cpp_filename(self) -> str:
        """Get the output {srn_model}.cpp filename."""
        return pathlib.Path(self.srn_hpp_filename).with_suffix(".cpp")

    @property
    def srn_hpp_filename(self) -> str:
        """Get the output {srn_model}.hpp filename."""
        return f"{self._model_name}{ODE_SUFFIX}And{SRN_SUFFIX}Model.hpp"

    # -- PRIVATE ---------------------------------------

    def _generate(self) -> None:
        """Generate Chaste code for the SRN model."""
        self._generate_srn_hpp()
        self._generate_srn_cpp()

    def _generate_srn_cpp(self) -> None:
        """Generate the cpp file code for the SRN model."""
        template = self._get_template("srn.cpp")
        vars = self._get_srn_cpp_vars()
        code = template.render(vars)
        self._outputs[self.srn_cpp_filename] = code

    def _generate_srn_hpp(self) -> None:
        """Generate the hpp file code for the SRN model."""
        template = self._get_template("srn.hpp")
        vars = self._get_hpp_vars(SRN_SUFFIX, self.srn_hpp_filename)
        code = template.render(vars)
        self._outputs[self.srn_hpp_filename] = code

    def _get_srn_cpp_vars(self) -> dict[str, str]:
        """Generate the template variables for the SRN model cpp file.

        return: The generated source file as a string.
        """
        # Sub-templates
        function_impl_template = self._get_template("srn/cpp/function_impl.cpp")
        state_param_template = self._get_template("srn/cpp/state_param.cpp")

        # Compartments
        compartment_inits = []
        for compartment in self._compartments:
            c_id = compartment.getId()
            c_var = self._get_varname(compartment)
            c_size = compartment.getSize()

            compartment_inits.append(f"{c_id} = {c_size}; // {c_var}")

        compartment_init_str = f"\n{TAB}".join(compartment_inits)

        # Rules
        rule_defs = []

        for rule in self._rules:
            r_id = rule.getId()
            formula = convert_formula(rule.getFormula())

            rule_defs.append(f"{r_id} = {formula};")

        rule_def_str = f"\n{TAB}".join(rule_defs)

        # Species
        # time_multiplier = self._get_timescale_multiplier()
        species_defaults = []
        state_vars = []
        state_params = []
        ode_defs = []
        ode_timescale_defs = []
        species_ode_inits = []
        for i, species in enumerate(self._species):
            s_id = species.getId()
            s_name = self._get_name(species)
            s_var = self._get_varname(species)
            s_conc = get_species_concentration(species)
            s_units = species.getSubstanceUnits()

            species_defaults.append(f"SetDefaultInitialCondition({i}, {s_conc}); // {s_var}")

            if self._is_state_parameter(species):
                state_param = state_param_template.render(
                    par_id=s_id,
                    par_num=len(state_params),
                    par_name=s_var,
                )
                state_params.append(state_param)
            else:
                state_vars.append(f"double {s_id} = rY[{i}]; // {s_var}")

            # ODE system
            c_id = species.getCompartment()
            compartment = self._compartments.get(c_id)
            c_var = self._get_varname(compartment)
            n = len(ode_defs)

            if s_id in self._odes_dict:
                if not species.getBoundaryCondition():
                    rhs = self._odes_dict[s_id]
                    ode_def = f"rDY[{n}] = ({rhs}) / {c_var}; // d{s_var}/dt"
                    ode_defs.append(ode_def)

                # TODO: Handle time scaling
                # if time_multiplier != 1.0:
                #     # This does not include species defined in algebraic rules
                #     ode_timescale_defs.append(f"rDY[{n}] *= {time_multiplier};")

            elif (s_id == "drag") or (s_var == "drag"):
                rhs = f"drag - rY[{i}]"
                ode_def = f"rDY[{n}] = ({rhs}) / {c_var}; // d{s_var}/dt;"
                ode_defs.append(ode_def)

            elif s_id in self._rules_dict:
                # Species defined by algebraic rules are not in odes_dict

                # Assuming these are assignments where variables are added together to represent a total
                rule_string = self._rules_dict[s_id]
                rule_list = rule_string.split(" ")

                rhs_list = []
                for rs_id in rule_list:
                    index = get_index_by_id(rs_id, self._species)
                    if index is not None:
                        rhs_list.append(f"rDY[{index}]")
                rhs = " + ".join(rhs_list)

                if rhs:
                    ode_def = f"rDY[{n}] = ({rhs}) / {c_var}; // d{s_var}/dt;"
                    ode_defs.append(ode_def)

            # TODO: include other rules

            # Initial conditions

            # If there's a compartment, then we would have normalised
            # the ODE, so declare it as non-dimensional
            init_units = "non-dim" if compartment else s_units

            if (s_id in self._odes_dict) or (s_id in self._rules_dict):
                species_ode_inits.append(f'this->mVariableNames.push_back("{s_name}");')
                species_ode_inits.append(f'this->mVariableUnits.push_back("{init_units}");')
                species_ode_inits.append(f"this->mInitialConditions.push_back({s_conc});")

            elif self._is_state_parameter(species):
                species_ode_inits.append(f'this->mParameterNames.push_back("{s_name}");')
                species_ode_inits.append(f'this->mParameterUnits.push_back("{init_units}");')

        species_defaults_str = f"\n{TAB}".join(species_defaults)
        species_state_param_def_str = f"\n{TAB}".join(state_params)
        state_var_def_str = f"\n{TAB}".join(state_vars)
        ode_def_str = f"\n{TAB}".join(ode_defs)
        ode_timescale_def_str = f"\n{TAB}".join(ode_timescale_defs)
        species_ode_init_str = f"\n{TAB}".join(species_ode_inits)

        # Parameters
        param_defaults = []
        param_inits = []
        state_params = []
        param_ode_inits = []
        for param in self._parameters:
            p_id = param.getId()
            p_var = self._get_varname(param)

            unset = (not param.isSetValue()) and (p_id not in self._rules_dict)
            special = p_var and any(x in p_var for x in ["wnt", "gamma", "ComplexTransit"])

            if unset or special:
                def_val = 0.0
                if special and "gamma" in p_var or "ComplexTransit" in p_var:
                    def_val = 1.0
                param_defaults.append(f"this->mParameters.push_back({def_val}); // {p_var}")

            p_value = param.getValue()
            param_inits.append(f"{p_id} = {p_value};")

            if self._is_state_parameter(param):
                state_param = state_param_template.render(
                    par_id=p_id,
                    par_num=len(state_params),
                    par_name=p_var,
                )
                state_params.append(state_param)

            # Parameters without set values must be externally defined
            if unset or special:
                init_units = species.getUnits() if param.isSetUnits() else "non-dim"
                param_ode_inits.append(f'this->mParameterNames.push_back("{p_var}");')
                param_ode_inits.append(f'this->mParameterUnits.push_back("{init_units}");')

        parameter_defaults_str = f"\n{TAB}".join(param_defaults)
        parameter_init_str = f"\n{TAB}".join(param_inits)
        parameter_state_param_def_str = f"\n{TAB}".join(state_params)
        parameter_ode_init_str = f"\n{TAB}".join(param_ode_inits)

        # Reactions
        reaction_defs = []
        for reaction in self._reactions:
            r_id = reaction.getId()
            r_name = reaction.getName()
            r_var = self._get_varname(reaction)

            kinetic_law = reaction.getKineticLaw()
            r_formula = convert_formula(kinetic_law.getFormula())

            r_params = kinetic_law.getListOfParameters()
            r_param_defs = []
            for r_param in r_params:
                rp_var = self._get_varname(r_param)
                rp_value = r_param.getValue()

                r_param_defs.append(f"double {rp_var} = {rp_value};")

            rparam_defs_str = f"\n{TAB}".join(r_param_defs)

            reaction_defs.append(f"// {r_name}")
            reaction_defs.append(rparam_defs_str)
            reaction_defs.append(f"double {r_var} = {r_formula};")

        reaction_def_str = f"\n\n{TAB}".join(reaction_defs)

        # Function Definitions
        functions = self._function_definitions
        function_impls = []
        for fn in functions:
            fn_id = fn.getId()
            args_list = get_function_definition_arguments(fn)
            body_cpp = convert_function_body(fn.getBody())
            impl = function_impl_template.render(
                ode_name=self._ode_class_name,
                fn=fn_id,
                fn_args=", ".join(args_list),
                fn_body_cpp=body_cpp,
            )
            function_impls.append(impl)
        functions_impl_str = "\n".join(function_impls)

        # Events
        num_events = self._model.getNumEvents()
        event_vector_init_str = ""
        if num_events > 0:
            event_vector_init_str = f"eventsSatisfied.resize({num_events}, false);"

        # Create inputs for the source file template
        cpp_vars = dict(
            compartment_init=compartment_init_str,
            event_vector_init=event_vector_init_str,
            functions_impl=functions_impl_str,
            model_hpp_file=self.srn_hpp_filename,
            model_name=self._model_name,
            model_suffix=SRN_SUFFIX,
            ode_def=ode_def_str,
            ode_class_name=self._ode_class_name,
            ode_timescale_def=ode_timescale_def_str,
            parameter_defaults=parameter_defaults_str,
            parameter_init=parameter_init_str,
            parameter_ode_init=parameter_ode_init_str,
            parameter_state_param_def=parameter_state_param_def_str,
            reaction_def=reaction_def_str,
            rule_def=rule_def_str,
            num_state_vars=self._num_state_vars,
            species_defaults=species_defaults_str,
            species_ode_init=species_ode_init_str,
            species_state_param_def=species_state_param_def_str,
            state_var_def=state_var_def_str,
        )

        return cpp_vars
