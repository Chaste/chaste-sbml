"""Unit tests for ModelBuilder: SBML translation and name-conflict detection."""

import pytest

from chaste_sbml._config import ROOT_DIR, EquationType, VarType
from chaste_sbml._model_builder import ModelBuilder
from chaste_sbml._names import NameConflictError, NameManager
from chaste_sbml._records import (
    DerivedQuantity,
    Equation,
    EventAssignment,
    Function,
    InitialAssignment,
    LocalParameter,
    Parameter,
    Reaction,
    StateVariable,
)
from chaste_sbml._sbml_loader import load_sbml_model

REFERENCE = ROOT_DIR / "SbmlRefModels" / "src" / "reference"


def _build(sbml_path) -> ModelBuilder:
    """Load, resolve names for and build a real SBML model, returning the ModelBuilder."""
    sbml_model = load_sbml_model(str(sbml_path))
    names = NameManager(sbml_model)
    names.resolve_real_id_conflicts()
    builder = ModelBuilder(sbml_model, names, "Test")
    builder.build()
    return builder


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

    def _state_var(i):
        return i if isinstance(i, StateVariable) else _make_state_variable(id=i, derivative_id=f"d_{i}_dt")

    def _derived(i):
        if isinstance(i, DerivedQuantity):
            return i
        if isinstance(i, dict):
            return _make_derived_quantity(**i)
        return _make_derived_quantity(id=i)

    builder = _builder_without_init()
    builder._model_name = "Demo"
    builder._parameters = [_make_parameter(id=i) for i in parameters]
    builder._state_variables = [_state_var(i) for i in state_variables]
    builder._derived_quantities = [_derived(i) for i in derived_quantities]
    builder._stoichiometry_variables = [_make_parameter(id=i) for i in stoichiometry_variables]
    builder._reactions = [Reaction(index=0, id=i, label="") for i in reactions]
    builder._functions = [Function(id=i, label="", index=0, args="", body="") for i in functions]
    builder._initial_assignments = [InitialAssignment(id=i, label="", var=i) for i in initial_assignments]
    # A single reaction equation carrying the given local-parameter ids.
    builder._equations = (
        [
            Equation(
                var="J1",
                math=None,
                local_parameters=[LocalParameter(id=i, label="", value="0") for i in local_parameters],
            )
        ]
        if local_parameters
        else []
    )
    return builder


def _make_parameter(*, id, is_const=False):
    return Parameter(index=0, id=id, is_const=is_const, label="", initial_value=None, units="")


def _make_state_variable(*, id, derivative_id):
    return StateVariable(index=0, id=id, derivative_id=derivative_id, label="", initial_value=None, units="")


def _make_derived_quantity(*, id, is_conversion=False, is_reaction=False):
    return DerivedQuantity(
        id=id, label="", index=0, initial_value=None, is_conversion=is_conversion, is_reaction=is_reaction, units=""
    )


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


def test_build_extracts_state_variables_and_derivatives():
    """Each SBML species with an ODE becomes a state variable with a d_<id>_dt derivative."""
    builder = _build(REFERENCE / "Goldbeter1991" / "Goldbeter1991.xml")
    state_vars = {s.id: s.derivative_id for s in builder._state_variables}
    assert state_vars == {"C": "d_C_dt", "M": "d_M_dt", "X": "d_X_dt"}


def test_build_extracts_parameters_and_reactions():
    """Global parameters and reactions are collected under their SBML ids."""
    builder = _build(REFERENCE / "Goldbeter1991" / "Goldbeter1991.xml")
    assert {p.id for p in builder._parameters} == {"VM1", "VM3", "Kc"}
    assert {r.id for r in builder._reactions} == {f"reaction{i}" for i in range(1, 8)}


def test_build_adds_amount_conversions_for_amount_species():
    """Amount species get an amt__<id> conversion derived quantity."""
    builder = _build(REFERENCE / "Goldbeter1991" / "Goldbeter1991.xml")
    conversions = {d.id for d in builder._derived_quantities if d.is_conversion}
    assert conversions == {"amt__C", "amt__M", "amt__X"}


def test_build_produces_one_derivative_equation_per_state_variable():
    """A derivative equation is emitted for each state variable."""
    builder = _build(REFERENCE / "Goldbeter1991" / "Goldbeter1991.xml")
    deriv_vars = {e.var for e in builder._equations if e.type == EquationType.DERIVATIVE}
    assert deriv_vars == {"d_C_dt", "d_M_dt", "d_X_dt"}


def test_build_extracts_cell_division_event():
    """An event is translated with its trigger and state-variable assignment."""
    builder = _build(REFERENCE / "TysonNovak2001" / "TysonNovak2001.xml")
    assert len(builder._events) == 1
    event = builder._events[0]
    assert event.trigger == "CycB < 0.1"
    assert event.assignments == [EventAssignment(index=3, lhs="m", rhs="m / 2.0", type=VarType.STATE_VARIABLE)]
