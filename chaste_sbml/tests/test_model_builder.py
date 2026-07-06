"""Unit tests for ModelBuilder: SBML translation and name-conflict detection."""

import libsbml
import pytest
from libsbml import parseL3Formula

from chaste_sbml._config import (
    CHASTE_PREFIX,
    PLACEHOLDER_STATE_ID,
    PREFIX_SEP,
    ROOT_DIR,
    DerivedQuantityKind,
    EquationType,
    VarType,
)
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
    Rule,
    StateVariable,
)
from chaste_sbml._sbml_reader import load_sbml_model

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


def _make_derived_quantity(*, id, kind=DerivedQuantityKind.NORMAL):
    return DerivedQuantity(id=id, label="", index=0, initial_value=None, units="", kind=kind)


def test_check_name_conflicts_passes_for_clean_model():
    """A model with distinct, valid ids and unique derivatives raises nothing."""
    builder = _builder_with_names(parameters=["VM1"], state_variables=["C", "M"], reactions=["reaction1"])
    builder._check_name_conflicts()  # should not raise


def test_check_name_conflicts_flags_synthetic_vs_real_collision():
    """A real id equal to a synthesised amount-conversion name is caught."""
    builder = _builder_with_names(
        state_variables=["X"],
        derived_quantities=[{"id": "amt__X", "kind": DerivedQuantityKind.CONVERSION}],
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
        derived_quantities=[{"id": "reaction1", "kind": DerivedQuantityKind.REACTION}],
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


PLACEHOLDER_ID = PREFIX_SEP.join([CHASTE_PREFIX, PLACEHOLDER_STATE_ID])


def _no_ode_model(tmp_path):
    """Write a no-ODE SBML model (constant species + time-dependent assignment rule) to disk.

    The model has no reactions, ODEs or rate rules, so it has no state variables of its own.
    """
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    model.setId("NoOde")
    compartment = model.createCompartment()
    compartment.setId("cell")
    compartment.setConstant(True)
    compartment.setSize(1.0)
    compartment.setSpatialDimensions(3)
    # A parameter that is an explicit function of time via an assignment rule.
    parameter = model.createParameter()
    parameter.setId("P")
    parameter.setConstant(False)
    parameter.setValue(0.0)
    rule = model.createAssignmentRule()
    rule.setVariable("P")
    rule.setMath(parseL3Formula("2 * time"))
    # A constant species.
    species = model.createSpecies()
    species.setId("S")
    species.setCompartment("cell")
    species.setConstant(True)
    species.setBoundaryCondition(True)
    species.setHasOnlySubstanceUnits(False)
    species.setInitialConcentration(3.0)

    path = tmp_path / "NoOde.xml"
    libsbml.writeSBMLToFile(doc, str(path))
    return path


def test_build_synthesises_placeholder_state_variable_for_no_ode_model(tmp_path):
    """A model with no ODEs gets a single zero-derivative placeholder state variable."""
    builder = _build(_no_ode_model(tmp_path))

    assert [s.id for s in builder._state_variables] == [PLACEHOLDER_ID]
    placeholder = builder._state_variables[0]
    assert placeholder.initial_value == 0.0

    deriv_equations = [e for e in builder._equations if e.type == EquationType.DERIVATIVE]
    assert len(deriv_equations) == 1
    assert deriv_equations[0].var == placeholder.derivative_id
    assert deriv_equations[0].rhs == "0.0"

    # The placeholder's member must be initialised to zero: Initialise() reads it to set the
    # solver's initial condition, so without this equation the state would be uninitialised memory.
    init_equations = [e for e in builder._equations if e.type == EquationType.INITIAL_VALUE]
    assert placeholder.id in [e.var for e in init_equations]
    assert next(e for e in init_equations if e.var == placeholder.id).rhs == "0.0"


def test_build_no_ode_model_preserves_real_outputs(tmp_path):
    """The placeholder does not displace the model's real outputs (assignment rule, species)."""
    builder = _build(_no_ode_model(tmp_path))

    # The time-dependent assignment-rule parameter is still emitted as a derived quantity, and
    # the constant boundary species is still emitted (as a parameter). Neither is the placeholder.
    assert "P" in {d.id for d in builder._derived_quantities}
    assert "S" in {p.id for p in builder._parameters}


def test_build_event_modified_constant_species_is_a_parameter(tmp_path):
    """An otherwise-constant species modified by an event is a parameter, not a derived quantity.

    A derived quantity is recomputed from the current member each step, so an event's change would
    be reported at every earlier time point; a parameter is time-resolved in the recorded solution.
    """
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    compartment = model.createCompartment()
    compartment.setId("C")
    compartment.setConstant(True)
    compartment.setSize(1.0)
    compartment.setSpatialDimensions(3)
    species = model.createSpecies()
    species.setId("S")
    species.setCompartment("C")
    species.setConstant(False)
    species.setBoundaryCondition(False)
    species.setHasOnlySubstanceUnits(False)
    species.setInitialAmount(0.0)
    event = model.createEvent()
    event.setUseValuesFromTriggerTime(True)
    trigger = event.createTrigger()
    trigger.setMath(parseL3Formula("time >= 1"))
    trigger.setInitialValue(True)
    trigger.setPersistent(True)
    assignment = event.createEventAssignment()
    assignment.setVariable("S")
    assignment.setMath(parseL3Formula("5"))
    path = tmp_path / "event_species.xml"
    libsbml.writeSBMLToFile(doc, str(path))

    builder = _build(path)

    assert "S" in {p.id for p in builder._parameters}
    assert "S" not in {d.id for d in builder._derived_quantities}


def test_build_does_not_add_placeholder_when_model_has_odes():
    """A model with genuine ODEs is left untouched: no placeholder state variable is added."""
    builder = _build(REFERENCE / "Goldbeter1991" / "Goldbeter1991.xml")
    assert PLACEHOLDER_ID not in {s.id for s in builder._state_variables}


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
    conversions = {d.id for d in builder._derived_quantities if d.kind == DerivedQuantityKind.CONVERSION}
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


def _event_builder():
    """A bare builder with an empty model, wired just enough for _format_events."""
    builder = _builder_without_init()
    doc = libsbml.SBMLDocument(3, 2)
    builder._sbml_model = doc.createModel()
    builder._keep_doc = doc  # keep the owning document (and thus the model) alive
    builder._sbml_events = builder._sbml_model.getListOfEvents()
    builder._sbml_function_definitions = builder._sbml_model.getListOfFunctionDefinitions()
    builder._state_variables = []
    builder._parameters = []
    builder._derived_quantities = []
    builder._reactions = []
    builder._events = []
    builder._functions = []
    return builder


def test_format_events_skips_event_with_missing_trigger_math():
    """An event whose trigger has no MathML can never fire, so no event is emitted."""
    builder = _event_builder()
    event = builder._sbml_model.createEvent()
    event.createTrigger()  # trigger element present but with no math
    event.createEventAssignment().setVariable("p")

    builder._format_events()

    assert builder._events == []


def test_format_events_skips_event_with_no_trigger():
    """An event with no trigger element at all can never fire, so no event is emitted."""
    builder = _event_builder()
    builder._sbml_model.createEvent()  # no trigger

    builder._format_events()

    assert builder._events == []


def test_reject_unsupported_delay_allows_missing_delay_math():
    """A delay element with no MathML is treated as no delay, so it is not rejected."""
    builder = _event_builder()
    event = builder._sbml_model.createEvent()
    event.createDelay()  # delay element present but with no math

    builder._reject_unsupported_delay(event)  # should not raise


def test_reject_unsupported_packages_rejects_fbc():
    """A model using the flux-balance-constraints (fbc) package is rejected."""
    doc = libsbml.SBMLDocument(libsbml.SBMLNamespaces(3, 1, "fbc", 2))
    builder = _builder_without_init()
    builder._sbml_model = doc.createModel()
    builder._keep_doc = doc

    with pytest.raises(NotImplementedError, match="fbc"):
        builder._reject_unsupported_packages()


def test_reject_unsupported_packages_allows_plain_model():
    """A model without an unsupported package is not rejected."""
    doc = libsbml.SBMLDocument(3, 2)
    builder = _builder_without_init()
    builder._sbml_model = doc.createModel()
    builder._keep_doc = doc

    builder._reject_unsupported_packages()  # should not raise


def test_format_function_definitions_skips_missing_body():
    """A function definition with no MathML body is skipped rather than crashing."""
    builder = _event_builder()
    builder._sbml_model.createFunctionDefinition().setId("foo")  # no body

    builder._format_function_definitions()

    assert builder._functions == []


def _distance_builder() -> ModelBuilder:
    """A bare builder wired just enough for _event_trigger_distance and trigger formatting."""
    builder = _builder_without_init()
    doc = libsbml.SBMLDocument(3, 2)
    builder._sbml_model = doc.createModel()
    builder._keep_doc = doc  # keep the owning document (and thus the model) alive
    builder._state_variables = []
    builder._parameters = []
    builder._derived_quantities = []
    builder._reactions = []
    builder._functions = []
    return builder


def _trigger_event(builder: ModelBuilder, formula: str):
    """Create an event on the builder's model whose trigger is ``formula``."""
    event = builder._sbml_model.createEvent()
    event.createTrigger().setMath(parseL3Formula(formula))
    return event


@pytest.mark.parametrize(
    ("trigger", "expected"),
    [
        ("A < 5", "(5.0) - (A) - std::numeric_limits<double>::epsilon()"),
        ("A > 5", "(A) - (5.0) - std::numeric_limits<double>::epsilon()"),
        ("A >= 5", "(A) - (5.0) - std::numeric_limits<double>::epsilon()"),
        ("A <= 5", "(5.0) - (A) - std::numeric_limits<double>::epsilon()"),
        ("A == 5", "-std::abs((A) - (5.0))"),
        ("A != 5", "std::abs((A) - (5.0)) - std::numeric_limits<double>::epsilon()"),
    ],
)
def test_event_trigger_distance_relational_operators(trigger, expected):
    """Each relational operator yields the signed distance CVODE roots on (only < is in the ref models)."""
    builder = _distance_builder()
    event = _trigger_event(builder, trigger)
    assert builder._event_trigger_distance(event) == expected


def test_event_trigger_distance_non_relational_is_constant():
    """A non-relational trigger has no crossing to locate, so the distance is the constant 1.0."""
    builder = _distance_builder()
    event = _trigger_event(builder, "and(A > 1, B > 1)")
    assert builder._event_trigger_distance(event) == "1.0"


def _resize_builder(*, assignment_rules=()) -> ModelBuilder:
    """A bare builder with one compartment C and one concentration species S in it."""
    builder = _builder_without_init()
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    model.createCompartment().setId("C")
    species = model.createSpecies()
    species.setId("S")
    species.setCompartment("C")

    builder._sbml_model = model
    builder._keep_doc = doc
    builder._sbml_compartments = model.getListOfCompartments()
    builder._sbml_species = model.getListOfSpecies()
    builder._assignment_rules = list(assignment_rules)
    builder._state_variables = [_make_state_variable(id="S", derivative_id="d_S_dt")]
    builder._parameters = []
    builder._derived_quantities = []
    builder._reactions = []
    builder._functions = []
    builder._index_of = {"S": 0}
    return builder


def test_compensate_compartment_resizes_directly_resized_compartment():
    """A species in a directly-resized compartment gets a compensating amount-conserving assignment."""
    builder = _resize_builder()
    event = builder._sbml_model.createEvent()
    assignments = [EventAssignment(index=0, lhs="C", rhs="2.0", type=VarType.PARAMETER)]

    builder._compensate_compartment_resizes(event, assignments)

    added = [a for a in assignments if a.lhs == "S"]
    assert added == [EventAssignment(index=0, lhs="S", rhs="S * C / (2.0)", type=VarType.STATE_VARIABLE)]


def test_compensate_compartment_resizes_indirect_assignment_rule():
    """A compartment resized indirectly (rule C = fakeC, event assigns fakeC) still rescales its species."""
    builder = _resize_builder(assignment_rules=[Rule(id="r", label="", var="C", math=parseL3Formula("fakeC"))])
    event = builder._sbml_model.createEvent()
    event_assignment = event.createEventAssignment()
    event_assignment.setVariable("fakeC")
    event_assignment.setMath(parseL3Formula("3.0"))
    assignments = [EventAssignment(index=1, lhs="fakeC", rhs="3", type=VarType.PARAMETER)]

    builder._compensate_compartment_resizes(event, assignments)

    added = [a for a in assignments if a.lhs == "S"]
    assert added == [EventAssignment(index=0, lhs="S", rhs="S * C / (3.0)", type=VarType.STATE_VARIABLE)]


def test_compensate_compartment_resizes_scales_existing_species_assignment():
    """When the event already assigns the species, its rhs is scaled rather than duplicated."""
    builder = _resize_builder()
    event = builder._sbml_model.createEvent()
    species_assignment = EventAssignment(index=0, lhs="S", rhs="0.5", type=VarType.STATE_VARIABLE)
    assignments = [EventAssignment(index=1, lhs="C", rhs="2.0", type=VarType.PARAMETER), species_assignment]

    builder._compensate_compartment_resizes(event, assignments)

    assert species_assignment.rhs == "(0.5) * C / (2.0)"
    assert [a.lhs for a in assignments] == ["C", "S"]  # not duplicated


def test_compensate_compartment_resizes_skips_amount_species():
    """An amount-tracked species conserves its amount automatically, so it is not rescaled."""
    builder = _resize_builder()
    builder._sbml_model.getSpecies("S").setHasOnlySubstanceUnits(True)
    event = builder._sbml_model.createEvent()
    assignments = [EventAssignment(index=0, lhs="C", rhs="2.0", type=VarType.PARAMETER)]

    builder._compensate_compartment_resizes(event, assignments)

    assert [a.lhs for a in assignments] == ["C"]
