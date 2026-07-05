"""Loading, normalising and reading an SBML model with libsbml.

Reads an SBML file, flattens hierarchical models and runs the conversions the
generator depends on, returning a libsbml ``Model`` ready to be processed.
python-libsbml keeps the owning document alive through the returned model, so
callers need only hold the model.
Also provides small readers over individual SBML elements (compartment size,
function-definition arguments) used while building the internal representation.
"""

from typing import TYPE_CHECKING

from libsbml import LIBSBML_OPERATION_SUCCESS, ConversionProperties, SBMLReader, formulaToString

if TYPE_CHECKING:
    from libsbml import Compartment, FunctionDefinition, Model


def load_sbml_model(sbml_file: str) -> "Model":
    """Read, comp-flatten and convert an SBML file into a ready-to-process libsbml Model.

    :param sbml_file: Path to the SBML file.
    :return: The loaded libsbml ``Model``.
    :raises ValueError: if the file has read errors, or flattening/conversion fails.
    """
    reader = SBMLReader()
    doc = reader.readSBMLFromFile(sbml_file)

    if doc.getNumErrors() > 0:
        doc.printErrors()
        raise ValueError(f"Errors found while reading SBML file: {sbml_file}")

    # Flatten hierarchical (comp package) models, composing submodels into
    # a flattened model with submodelId__element names.
    if doc.getPlugin("comp") is not None:
        flatten_props = ConversionProperties()
        flatten_props.addOption("flatten comp", True, "flatten comp")
        flatten_props.addOption("leave_ports", False)
        if doc.convert(flatten_props) != LIBSBML_OPERATION_SUCCESS:
            doc.printErrors()
            raise ValueError("Errors during comp flattening")

    # Run required conversions
    config = ConversionProperties()
    # Sort assignment rules in order of dependence.
    config.addOption("sortRules")
    config.addOption("removeUnusedUnits")
    # Convert initial assignments to initial values where possible
    config.addOption("expandInitialAssignments")
    # config.addOption('replaceReactions')
    # config.addOption('expandFunctionDefinitions')

    status = doc.convert(config)
    if status != LIBSBML_OPERATION_SUCCESS:
        doc.printErrors()
        raise ValueError("Errors during conversion")

    return doc.getModel()


def get_compartment_size(compartment: "Compartment") -> float:
    """Get a compartment size.

    :return: The compartment size.
    """
    if compartment.isSetSize():
        return compartment.getSize()
    return 1.0


def get_function_definition_arguments(fn_def: "FunctionDefinition") -> list[str]:
    """Get the list of arguments in a given function definition.

    :param fn_def: The function definition
    :return: List of arguments in the function definition
    """
    n = fn_def.getNumArguments()
    return [formulaToString(fn_def.getArgument(i)) for i in range(n)]
