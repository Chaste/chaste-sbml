from ._config import ODE_SYSTEM_SUFFIX, SRN_HEADER_GUARD_SUFFIX, SRN_MODEL_SUFFIX, TAB
from ._utils import get_index_by_id
from .chaste_model import ChasteModel
from .templates.srn.srn_cpp import function_impl_template, srn_cpp_template, state_param_template
from .templates.srn.srn_hpp import srn_hpp_template


class ChasteSRNModel(ChasteModel):

    def __init__(self, sbml_file: str, model_name: str = None, **kwargs) -> None:
        super().__init__(sbml_file, model_name, **kwargs)

        self.ode_system_name = self.class_name + ODE_SYSTEM_SUFFIX
        self.srn_model_name = self.class_name + SRN_MODEL_SUFFIX

        self.hpp_filename = f"{self.class_name}{SRN_MODEL_SUFFIX}.hpp"
        self.cpp_filename = f"{self.class_name}{SRN_MODEL_SUFFIX}.cpp"

        self.size = len(self.odes_dict) + len(self.rules_dict)

    def generate_hpp(self) -> str:
        """Generate the Chaste header for an SRN model from SBML data.

        return: The generated header file as a string.
        """

        # Get inputs for the header file template
        header_guard_str = self.class_name.upper() + SRN_HEADER_GUARD_SUFFIX

        compartments = self.model.getListOfCompartments()
        compartment_decls_str = f"\n{TAB}".join(
            [f"double {c.getId()}; // {self.get_varname(c)}" for c in compartments]
        )

        parameters = self.model.getListOfParameters()
        parameter_decls_str = f"\n{TAB}".join(
            [f"double {p.getId()}; // {self.get_varname(p)}" for p in parameters]
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
            header_guard=header_guard_str,
            ode_system_name=self.ode_system_name,
            srn_name=self.srn_model_name,
            size=self.size,
            compartment_decls=compartment_decls_str,
            parameter_decls=parameter_decls_str,
            function_decls=function_decls_str,
        )

        return hpp

    def generate_cpp(self) -> str:
        """Generate the Chaste source for an SRN model from SBML data.

        return: The generated source file as a string.
        """

        # Compartments
        compartments = self.model.getListOfCompartments()
        compartment_inits = []
        for compartment in compartments:
            c_id = compartment.getId()
            c_var = self.get_varname(compartment)
            c_size = compartment.getSize()

            compartment_inits.append(f"{c_id} = {c_size}; // {c_var}")

        compartment_init_str = f"\n{TAB}".join(compartment_inits)

        # Rules
        rules = self.model.getListOfRules()
        rule_defs = []

        for rule in rules:
            r_id = rule.getId()
            formula = self.convert_formula(rule.getFormula())

            rule_defs.append(f"{r_id} = {formula};")

        rule_def_str = f"\n{TAB}".join(rule_defs)

        # Species
        specii = self.model.getListOfSpecies()
        time_multiplier = self.get_timescale_multiplier()
        species_defaults = []
        state_vars = []
        state_params = []
        ode_defs = []
        ode_timescale_defs = []
        species_ode_inits = []
        for i, species in enumerate(specii):
            s_id = species.getId()
            s_var = self.get_varname(species)
            s_conc = self.get_species_concentration(species)
            s_units = species.getSubstanceUnits()

            species_defaults.append(f"SetDefaultInitialCondition({i}, {s_conc}); // {s_var}")

            if self.is_state_parameter(species):
                state_param = state_param_template.format(
                    par_id=s_id,
                    par_num=len(state_params),
                    par_name=s_var,
                )
                state_params.append(state_param)
            else:
                state_vars.append(f"double {s_id} = rY[{i}]; // {s_var}")

            # ODE system
            c_id = species.getCompartment()
            compartment = compartments.get(c_id)
            c_var = self.get_varname(compartment)
            n = len(ode_defs)

            if s_id in self.odes_dict:
                if not species.getBoundaryCondition():
                    rhs = self.odes_dict[s_id]
                    ode_def = f"rDY[{n}] = ({rhs}) / {c_var}; // d{s_var}/dt"
                    ode_defs.append(ode_def)

                if time_multiplier != 1.0:
                    # This does not include species defined in algebraic rules
                    ode_timescale_defs.append(f"rDY[{n}] *= {time_multiplier};")

            elif (s_id == "drag") or (s_var == "drag"):
                rhs = f"drag - rY[{i}]"
                ode_def = f"rDY[{n}] = ({rhs}) / {c_var}; // d{s_var}/dt;"
                ode_defs.append(ode_def)

            elif s_id in self.rules_dict:
                # Species defined by algebraic rules are not in odes_dict

                # Assuming these are assignments where variables are added together to represent a total
                rule_string = self.rules_dict[s_id]
                rule_list = rule_string.split(" ")

                rhs_list = []
                for rs_id in rule_list:
                    index = get_index_by_id(rs_id, specii)
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

            if (s_id in self.odes_dict) or (s_id in self.rules_dict):
                species_ode_inits.append(f'this->mVariableNames.push_back("{s_var}");')
                species_ode_inits.append(f'this->mVariableUnits.push_back("{init_units}");')
                species_ode_inits.append(f"this->mInitialConditions.push_back({s_conc});")

            elif self.is_state_parameter(species):
                species_ode_inits.append(f'this->mParameterNames.push_back("{s_var}");')
                species_ode_inits.append(f'this->mParameterUnits.push_back("{init_units}");')

        species_defaults_str = f"\n{TAB}".join(species_defaults)
        species_state_param_def_str = f"\n{TAB}".join(state_params)
        state_var_def_str = f"\n{TAB}".join(state_vars)
        ode_def_str = f"\n{TAB}".join(ode_defs)
        ode_timescale_def_str = f"\n{TAB}".join(ode_timescale_defs)
        species_ode_init_str = f"\n{TAB}".join(species_ode_inits)

        # Parameters
        params = self.model.getListOfParameters()
        param_defaults = []
        param_inits = []
        state_params = []
        param_ode_inits = []
        for param in params:
            p_id = param.getId()
            p_var = self.get_varname(param)

            unset = (not param.isSetValue()) and (p_id not in self.rules_dict)
            special = p_var and any(x in p_var for x in ["wnt", "gamma", "ComplexTransit"])

            if unset or special:
                def_val = 0.0
                if special and "gamma" in p_var or "ComplexTransit" in p_var:
                    def_val = 1.0
                param_defaults.append(f"this->mParameters.push_back({def_val}); // {p_var}")

            p_value = param.getValue()
            param_inits.append(f"{p_id} = {p_value};")

            if self.is_state_parameter(param):
                state_param = state_param_template.format(
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
        reactions = self.model.getListOfReactions()
        reaction_defs = []
        for reaction in reactions:
            r_id = reaction.getId()
            r_name = reaction.getName()

            kinetic_law = reaction.getKineticLaw()
            r_formula = self.convert_formula(kinetic_law.getFormula())

            r_params = kinetic_law.getListOfParameters()
            r_param_defs = []
            for r_param in r_params:
                rp_id = r_param.getId()
                rp_value = r_param.getValue()

                r_param_defs.append(f"double {rp_id} = {rp_value};")

            rparam_defs_str = f"\n{TAB}".join(r_param_defs)

            reaction_defs.append(f"// {r_name}")
            reaction_defs.append(rparam_defs_str)
            reaction_defs.append(f"double {r_id} = {r_formula};")

        reaction_def_str = f"\n\n{TAB}".join(reaction_defs)

        # Function Definitions
        functions = self.model.getListOfFunctionDefinitions()
        function_impls = []
        for fn in functions:
            fn_id = fn.getId()
            args_list = self.get_function_definition_arguments(fn)
            body_cpp = self.convert_function_body(fn.getBody())
            impl = function_impl_template.format(
                ode_name=self.ode_system_name,
                fn=fn_id,
                fn_args=", ".join(args_list),
                fn_body_cpp=body_cpp,
            )
            function_impls.append(impl)
        functions_impl_str = "\n".join(function_impls)

        # Events
        num_events = self.model.getNumEvents()
        event_vector_init_str = ""
        if num_events > 0:
            event_vector_init_str = f"eventsSatisfied.resize({num_events}, false);"

        # Apply inputs to the source file template
        cpp = srn_cpp_template.format(
            srn_model_name=self.srn_model_name,
            ode_system_name=self.ode_system_name,
            model_name=self.class_name,
            size=self.size,
            species_defaults=species_defaults_str,
            parameter_defaults=parameter_defaults_str,
            compartment_init=compartment_init_str,
            parameter_init=parameter_init_str,
            event_vector_init=event_vector_init_str,
            state_var_def=state_var_def_str,
            species_state_param_def=species_state_param_def_str,
            parameter_state_param_def=parameter_state_param_def_str,
            functions_impl=functions_impl_str,
            ode_def=ode_def_str,
            ode_timescale_def=ode_timescale_def_str,
            species_ode_init=species_ode_init_str,
            parameter_ode_init=parameter_ode_init_str,
            rule_def=rule_def_str,
            reaction_def=reaction_def_str,
        )

        return cpp
