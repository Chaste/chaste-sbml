"""Unit tests for ChasteSbmlModel internals."""

import pytest

from chaste_sbml._names import NameConflictError
from chaste_sbml.chaste_sbml_model import ChasteSbmlModel


def _model_without_init() -> ChasteSbmlModel:
    """Create an instance without running __init__ for pure string helper tests."""
    return object.__new__(ChasteSbmlModel)


def _model_with_names(
    *,
    parameters=(),
    state_variables=(),
    derived_quantities=(),
    stoichiometry_variables=(),
    reactions=(),
    functions=(),
    initial_assignments=(),
    local_parameters=(),
) -> ChasteSbmlModel:
    """Build a bare model with just the id collections _check_name_conflicts inspects."""
    model = _model_without_init()
    model._model_name = "Demo"
    model._parameters = [{"id": i} for i in parameters]
    model._state_variables = [
        {"id": i, "derivative_id": f"d_{i}_dt"} if isinstance(i, str) else i for i in state_variables
    ]
    model._derived_quantities = [
        {"id": i, "is_conversion": False, "is_reaction": False} if isinstance(i, str) else i for i in derived_quantities
    ]
    model._stoichiometry_variables = [{"id": i} for i in stoichiometry_variables]
    model._reactions = [{"id": i} for i in reactions]
    model._functions = [{"id": i} for i in functions]
    model._initial_assignments = [{"id": i} for i in initial_assignments]
    # A single reaction equation carrying the given local-parameter ids.
    model._equations = (
        [{"var": "J1", "local_parameters": [{"id": i} for i in local_parameters]}] if local_parameters else []
    )
    return model


def test_check_name_conflicts_passes_for_clean_model():
    """A model with distinct, valid ids and unique derivatives raises nothing."""
    model = _model_with_names(parameters=["VM1"], state_variables=["C", "M"], reactions=["reaction1"])
    model._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_synthetic_vs_real_collision():
    """A real id equal to a synthesised amount-conversion name is caught."""
    model = _model_with_names(
        state_variables=["X"],
        derived_quantities=[{"id": "amt__X", "is_conversion": True, "is_reaction": False}],
        parameters=["amt__X"],
    )
    with pytest.raises(NameConflictError, match="amt__X"):
        model._check_name_conflicts()


def test_check_name_conflicts_flags_keyword_id():
    """An id that is a C++ keyword is caught."""
    model = _model_with_names(parameters=["double"], state_variables=["C"])
    with pytest.raises(NameConflictError, match="C\\+\\+ keyword"):
        model._check_name_conflicts()


def test_check_name_conflicts_flags_derivative_clash():
    """A state variable whose id collides with another's derivative name is caught."""
    model = _model_with_names(state_variables=["C", "d_C_dt"])
    with pytest.raises(NameConflictError, match="d_C_dt"):
        model._check_name_conflicts()


def test_check_name_conflicts_ignores_reaction_flux_duplicate():
    """A reaction also exposed as a derived-quantity output is not a false duplicate."""
    model = _model_with_names(
        state_variables=["C"],
        reactions=["reaction1"],
        derived_quantities=[{"id": "reaction1", "is_conversion": False, "is_reaction": True}],
    )
    model._check_name_conflicts()  # should not raise


def test_check_name_conflicts_ignores_initial_assignment_target():
    """An initial assignment to an existing parameter is not a false duplicate.

    An initial assignment's id is the variable it assigns to (emitted as ``var = <expr>``),
    not a new declaration, so it must not be counted as a separately-emitted identifier.
    """
    model = _model_with_names(parameters=["q"], state_variables=["C"], initial_assignments=["q"])
    model._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_keyword_local_parameter():
    """A kinetic-law local parameter that is a C++ keyword is caught."""
    model = _model_with_names(state_variables=["C"], reactions=["J1"], local_parameters=["int"])
    with pytest.raises(NameConflictError, match="local parameter 'int' is a C\\+\\+ keyword"):
        model._check_name_conflicts()


def test_check_name_conflicts_allows_ordinary_local_parameter():
    """A validly-named local parameter (even shadowing a global) is not flagged."""
    model = _model_with_names(parameters=["k1"], state_variables=["C"], reactions=["J1"], local_parameters=["k1"])
    model._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_local_parameter_shadowing_time():
    """A local parameter named 'time' would shadow the emitted time parameter and is flagged."""
    model = _model_with_names(state_variables=["C"], reactions=["J1"], local_parameters=["time"])
    with pytest.raises(NameConflictError, match="local parameter 'time' clashes with a reserved Chaste name"):
        model._check_name_conflicts()
