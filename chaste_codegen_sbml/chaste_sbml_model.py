"""Module for the ChasteSbmlModel class."""

import abc
import os
import pathlib
from typing import TYPE_CHECKING

from jinja2 import Environment, PackageLoader, select_autoescape
from libsbml import Parameter, SBMLReader, Species

from ._config import ODE_SUFFIX, SHORT_NAME_LEN
from ._utils import (
    convert_formula,
    get_function_definition_arguments,
    varname_camelcase,
    varname_sanitize,
)

if TYPE_CHECKING:
    from typing import Any

    from jinja2.environment import Template
    from libsbml import SBase


class ChasteSbmlModel:
    """Holds information about an SBML model for which Chaste code is to be generated."""

    __metaclass__ = abc.ABCMeta

    _jinja_env = Environment(
        loader=PackageLoader("chaste_codegen_sbml"),
        autoescape=select_autoescape(),
        trim_blocks=True,
        lstrip_blocks=True,
    )

    # -- PUBLIC ---------------------------------------

    def __init__(self, sbml: str, name: str = None, suffix: str = None) -> None:
        """Initialise the ChasteSbmlModel.

        :param sbml: The SBML file.
        :param name: The model name.
        :param suffix: The model type suffix e.g. "CellCycle".
        """
        self._sbml_file = sbml

        if name:
            self._model_name = name
        else:
            filename = os.path.splitext(os.path.basename(sbml))[0]
            self._model_name = varname_camelcase(filename).title()

        self._model_suffix = suffix
        self._ode_class_name = self._model_name + ODE_SUFFIX
        self._model_class_name = f"{self._model_name}{self._model_suffix}Model"
        self._wrapper_class_name = f"Sbml{self._model_suffix}WrapperModel"

        self._model = SBMLReader().readSBMLFromFile(sbml).getModel()
        self._compartments = self._model.getListOfCompartments()
        self._events = self._model.getListOfEvents()
        self._function_definitions = self._model.getListOfFunctionDefinitions()
        self._parameters = self._model.getListOfParameters()
        self._reactions = self._model.getListOfReactions()
        self._rules = self._model.getListOfRules()
        self._species = self._model.getListOfSpecies()
        self._unit_definitions = self._model.getListOfUnitDefinitions()

        self._varnames = {}

        self._odes_dict = None
        self._update_odes_dict()

        self._rules_dict = None
        self._update_rules_dict()

        self._num_state_vars = len([s for s in self._species if not self._is_state_parameter(s)])

        self._outputs = {}  # { filename: code }

    def write(self, output_directory=None):
        """Generate Chaste code and write to file."""
        # Generate the code
        self._generate()

        # Write the code to file
        if output_directory:
            root_dir = pathlib.Path(output_directory)
        else:
            root_dir = pathlib.Path().cwd()

        for filename, code in self._outputs.items():
            file_path = root_dir / filename
            with open(file_path, "w") as f:
                f.write(code)

    # -- PRIVATE ---------------------------------------

    @abc.abstractmethod
    def _generate(self) -> None:
        """Generate Chaste code for the model.
        This method should be implemented by subclasses.
        Generated code should be stored in the outputs using _add_output.
        """
        return

    def _add_output(self, filename: str, code: str) -> None:
        """Add generated code to the outputs dictionary.

        :param filename: The filename.
        :param code: The code.
        """
        self._outputs[filename] = code

    def _format_compartments(self) -> list[dict[str, str]]:
        """Get a list of compartment dictionaries for the model.

        :return: A list of compartment dictionaries.
        """
        return [
            {"id": c.getId(), "varname": self._get_varname(c), "size": c.getSize()}
            for c in self._compartments
        ]

    def _format_function_definitions(self) -> list[dict[str, str]]:
        """Get a list of function definition dictionaries for the model.

        :return: A list of function definition dictionaries.
        """
        return [
            {"id": fd.getId(), "args": get_function_definition_arguments(fd)}
            for fd in self._function_definitions
        ]

    def _format_header_guard(self, filename: str) -> str:
        """Get the header guard for a file.

        :param filename: The filename.
        :return: The header guard.
        """
        return filename.upper().replace(".", "_") + "_"

    def _format_parameters(self) -> list[dict[str, str]]:
        """Get a list of parameter dictionaries for the model.

        :return: A list of parameter dictionaries.
        """
        parameter_dicts = []
        state_index = 0
        for p in self._parameters:
            p_id = p.getId()
            p_name = self._get_name(p)
            p_var = self._get_varname(p)

            p_default = 0.0
            p_is_special = any(x in p_name for x in ["wnt", "gamma", "ComplexTransit"])
            if p_is_special:
                # Special strings in the parameter name
                # TODO: Review how to handle these special parameters
                if any(x in p_name for x in ["gamma", "ComplexTransit"]):
                    p_default = 1.0
                # 0.0 for wnt and everything else
            elif p.isSetValue():
                p_default = p.getValue()

            p_value = p.getValue() # TODO: if p.isSetValue() else p_default ?
            p_units = p.getUnits() if p.isSetUnits() else "non-dim"

            p_is_defined = (p.isSetValue() or p_id in self._rules_dict) and not p_is_special

            p_is_state = self._is_state_parameter(p)
            if p_is_state:
                p_state_index = state_index
                state_index += 1
            else:
                p_state_index = None

            parameter_dicts.append(
                {
                    "default": p_default,
                    "id": p_id,
                    "is_defined": p_is_defined,
                    "is_state": p_is_state,
                    "name": p_name,
                    "state_index": p_state_index,
                    "units": p_units,
                    "value": p_value,
                    "varname": p_var,
                }
            )

        return parameter_dicts

    def _format_rules(self) -> list[dict[str, str]]:
        """Get a list of rule dictionaries for the model.

        :return: A list of rule dictionaries.
        """
        return [{"id": r.getId(), "formula": convert_formula(r.getFormula())} for r in self._rules]

    def _get_hpp_vars(self, filename: str) -> dict[str, "Any"]:
        """Generate the template variables for the model hpp file.

        :param model_suffix: A suffix for the model type e.g. CellCycle.
        :param filename: The hpp filename for the model.
        return: The generated header file as a string.
        """
        return dict(
            compartments=self._format_compartments(),
            function_definitions=self._format_function_definitions(),
            header_guard=self._format_header_guard(filename),
            model_class_name=self._model_class_name,
            num_events=self._model.getNumEvents(),
            num_state_vars=self._num_state_vars,
            ode_class_name=self._ode_class_name,
            parameters=self._format_parameters(),
            wrapper_class_name=self._wrapper_class_name,
        )

    def _get_name(self, obj: "SBase") -> str:
        """Get the name of a libSBML object, or the ID if it doesn't have one.

        :param obj: The object.
        :return: The object name, or ID.
        """
        obj_name = obj.getName().strip()
        if obj_name:
            return obj_name
        return obj.getId()

    def _get_template(self, name: str) -> "Template":
        """Get a Jinja2 template.

        :param name: The template name.
        :return: The template object.
        """
        return self._jinja_env.get_template(name)

    def _get_timescale_multiplier(self) -> float:
        """Get the timescale multiplier.

        SBML uses seconds by default and Chaste uses hours.

        :return: The timescale multiplier.
        """
        for unit_def in self._unit_definitions:
            u_id = unit_def.getId()
            if u_id.lower() == "time":  # Do people ever call this something different?
                timescale = unit_def.getName().lower()
                if "minute" in timescale:
                    return 60.0
                elif "hour" in timescale:
                    return 1.0
        return 3600.0

    def _get_varname(self, obj: "SBase") -> str:
        """Get a suitable C++ variable name for a libSBML object.

        :param obj: The object.
        :return: The variable name.
        """
        obj_id = obj.getId()
        if obj_id in self._varnames:
            return self._varnames[obj_id]

        # Prefer the name if it is reasonably short, or shorter than the ID
        obj_name = varname_sanitize(self._get_name(obj))
        if 0 < len(obj_name) <= max(SHORT_NAME_LEN, len(obj_id)):
            var = obj_name
        else:
            var = obj_id

        # Check that all generated variable names are unique
        if var in self._varnames.values():
            i = 0
            while f"{var}_{i}" in self._varnames.values():
                i += 1
            var = f"{var}_{i}"

        self._varnames[obj_id] = var
        return var

    def _is_state_parameter(self, obj: "Species | Parameter") -> bool:
        """Check if a species or parameter is defined as a state parameter for Chaste.

        :param obj: The species or parameter to check.
        :return: True if defined as a parameter, False otherwise.
        """
        # Any species not defined by an ODE or rule is set as a state parameter
        if isinstance(obj, Species):
            species_id = obj.getId()
            return (species_id not in self._odes_dict) and (species_id not in self._rules_dict)

        if isinstance(obj, Parameter):
            # Also parameters with special strings in their name are state parameters.
            parameter_name = obj.getName()
            match_strings = ["wnt", "gamma", "ComplexTransit"]
            if parameter_name and any(s in parameter_name for s in match_strings):
                return True

            # Also parameters with unset values are state parameters.
            parameter_id = obj.getId()
            return (not obj.isSetValue()) and (parameter_id not in self._rules_dict)

        return False

    def _update_odes_dict(self) -> None:
        """Set the ODEs dictionary of equations corresponding to each species.

        Each ODE will essentially be the sum of the products minus the sum of
        the reactants divided by the compartment volume
        """
        self._odes_dict = {}
        for reaction in self._reactions:
            reaction_var = self._get_varname(reaction)

            # Decompose reaction into sum of products minus sum of reactants
            products = reaction.getListOfProducts()
            for product in products:
                # Get the species concerning the product
                species_id = product.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self._model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in self._odes_dict:
                    self._odes_dict[species_id] += " + " + reaction_var
                else:
                    self._odes_dict[species_id] = reaction_var

            reactants = reaction.getListOfReactants()
            for reactant in reactants:
                species_id = reactant.getSpecies()

                # TODO: Do we need to do something special with boundary conditions?
                # species = self._model.getSpecies(species_id)
                # if species.isSetBoundaryCondition() and not species.getBoundaryCondition():

                if species_id in self._odes_dict:
                    self._odes_dict[species_id] += " - " + reaction_var
                else:
                    self._odes_dict[species_id] = "-" + reaction_var

    def _update_rules_dict(self) -> None:
        """Set the dictionary of species defined by reaction rules."""
        self._rules_dict = {r.getId(): convert_formula(r.getFormula()) for r in self._rules}
