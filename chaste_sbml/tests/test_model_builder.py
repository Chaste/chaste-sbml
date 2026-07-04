"""Unit tests for ModelBuilder name-conflict detection."""

import pytest

from chaste_sbml._names import NameConflictError
from chaste_sbml._model_builder import ModelBuilder


def _builder_without_init() -> ModelBuilder:
    """Create a ModelBuilder without running __init__ for direct method tests."""
    return object.__new__(ModelBuilder)


def _builder_with_names(
    *,
    parameters=(),
    state_variables=(),
    derived_quantities=(),
    stoichiometry_variables=(),
    reactions=(),
    functions=(),
    initial_assignments=(),
    local_parameters=(),
) -> ModelBuilder:
    """Build a bare builder with just the id collections _check_name_conflicts inspects."""
    builder = _builder_without_init()
    builder._model_name = "Demo"
    builder._parameters = [{"id": i} for i in parameters]
    builder._state_variables = [
        {"id": i, "derivative_id": f"d_{i}_dt"} if isinstance(i, str) else i for i in state_variables
    ]
    builder._derived_quantities = [
        {"id": i, "is_conversion": False, "is_reaction": False} if isinstance(i, str) else i for i in derived_quantities
    ]
    builder._stoichiometry_variables = [{"id": i} for i in stoichiometry_variables]
    builder._reactions = [{"id": i} for i in reactions]
    builder._functions = [{"id": i} for i in functions]
    builder._initial_assignments = [{"id": i} for i in initial_assignments]
    # A single reaction equation carrying the given local-parameter ids.
    builder._equations = (
        [{"var": "J1", "local_parameters": [{"id": i} for i in local_parameters]}] if local_parameters else []
    )
    return builder


def test_check_name_conflicts_passes_for_clean_model():
    """A model with distinct, valid ids and unique derivatives raises nothing."""
    builder = _builder_with_names(parameters=["VM1"], state_variables=["C", "M"], reactions=["reaction1"])
    builder._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_synthetic_vs_real_collision():
    """A real id equal to a synthesised amount-conversion name is caught."""
    builder = _builder_with_names(
        state_variables=["X"],
        derived_quantities=[{"id": "amt__X", "is_conversion": True, "is_reaction": False}],
        parameters=["amt__X"],
    )
    with pytest.raises(NameConflictError, match="amt__X"):
        builder._check_name_conflicts()


def test_check_name_conflicts_flags_keyword_id():
    """An id that is a C++ keyword is caught."""
    builder = _builder_with_names(parameters=["double"], state_variables=["C"])
    with pytest.raises(NameConflictError, match="C\\+\\+ keyword"):
        builder._check_name_conflicts()


def test_check_name_conflicts_flags_derivative_clash():
    """A state variable whose id collides with another's derivative name is caught."""
    builder = _builder_with_names(state_variables=["C", "d_C_dt"])
    with pytest.raises(NameConflictError, match="d_C_dt"):
        builder._check_name_conflicts()


def test_check_name_conflicts_ignores_reaction_flux_duplicate():
    """A reaction also exposed as a derived-quantity output is not a false duplicate."""
    builder = _builder_with_names(
        state_variables=["C"],
        reactions=["reaction1"],
        derived_quantities=[{"id": "reaction1", "is_conversion": False, "is_reaction": True}],
    )
    builder._check_name_conflicts()  # should not raise


def test_check_name_conflicts_ignores_initial_assignment_target():
    """An initial assignment to an existing parameter is not a false duplicate.

    An initial assignment's id is the variable it assigns to (emitted as ``var = <expr>``),
    not a new declaration, so it must not be counted as a separately-emitted identifier.
    """
    builder = _builder_with_names(parameters=["q"], state_variables=["C"], initial_assignments=["q"])
    builder._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_keyword_local_parameter():
    """A kinetic-law local parameter that is a C++ keyword is caught."""
    builder = _builder_with_names(state_variables=["C"], reactions=["J1"], local_parameters=["int"])
    with pytest.raises(NameConflictError, match="local parameter 'int' is a C\\+\\+ keyword"):
        builder._check_name_conflicts()


def test_check_name_conflicts_allows_ordinary_local_parameter():
    """A validly-named local parameter (even shadowing a global) is not flagged."""
    builder = _builder_with_names(parameters=["k1"], state_variables=["C"], reactions=["J1"], local_parameters=["k1"])
    builder._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_local_parameter_shadowing_time():
    """A local parameter named 'time' would shadow the emitted time parameter and is flagged."""
    builder = _builder_with_names(state_variables=["C"], reactions=["J1"], local_parameters=["time"])
    with pytest.raises(NameConflictError, match="local parameter 'time' clashes with a reserved Chaste name"):
        builder._check_name_conflicts()
