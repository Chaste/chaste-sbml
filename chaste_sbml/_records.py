"""Typed records for the internal model representation.

ModelBuilder populates lists of these dataclasses and the Jinja templates consume them via
their attributes. Jinja falls back from ``record["field"]`` to ``record.field``, so the
templates keep their subscript syntax unchanged. Field names therefore mirror the keys the
templates read; changing one means updating the templates too.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from ._config import DerivedQuantityKind, EquationType, EventType, VarType

if TYPE_CHECKING:
    from libsbml import ASTNode


@dataclass
class Rule:
    """An SBML assignment or rate rule: ``var`` defined by ``math``."""

    id: str
    label: str
    var: str
    math: "ASTNode"


@dataclass
class InitialAssignment:
    """An SBML initial assignment: ``var`` initialised from ``math``."""

    id: str
    label: str
    var: str
    math: Optional["ASTNode"] = None


@dataclass
class StateVariable:
    """A time-integrated variable, paired with the id of its derivative."""

    index: int
    id: str
    derivative_id: str
    label: str
    initial_value: Optional[float]
    units: str
    # SBML id reported to Chaste; equals ``id`` unless the id was escaped for C++ (see NameManager).
    name: str = ""


@dataclass
class Parameter:
    """A model parameter (constant, or varied by a rule or event)."""

    index: int
    id: str
    is_const: bool
    label: str
    initial_value: Optional[float]
    units: str
    # SBML id reported to Chaste; equals ``id`` unless the id was escaped for C++ (see NameManager).
    name: str = ""


@dataclass
class DerivedQuantity:
    """A quantity recomputed from state each step (incl. amount/concentration conversions)."""

    id: str
    label: str
    index: int
    initial_value: Optional[float]
    units: str
    kind: DerivedQuantityKind = DerivedQuantityKind.NORMAL
    # SBML id reported to Chaste; equals ``id`` unless the id was escaped for C++ (see NameManager).
    name: str = ""


@dataclass
class Reaction:
    """A reaction flux, declared and computed as a member."""

    index: int
    id: str
    label: str


@dataclass
class Function:
    """A model function definition."""

    id: str
    label: str
    index: int
    args: str
    body: str


@dataclass
class LocalParameter:
    """A kinetic-law local parameter scoped to a single reaction equation."""

    id: str
    label: str
    value: str


@dataclass
class Equation:
    """A ``var = rhs`` equation; ``rhs`` is filled in from ``math`` during formatting."""

    var: str
    math: Optional["ASTNode"]
    type: EquationType = EquationType.UNKNOWN
    local_parameters: Optional[list[LocalParameter]] = None
    rhs: str = ""


@dataclass
class EventAssignment:
    """A single assignment an event applies to a variable when it fires."""

    index: int
    lhs: str
    rhs: str
    type: VarType


@dataclass
class Event:
    """An SBML event: a trigger and the assignments it applies when fired."""

    label: str
    index: int
    trigger: str
    assignments: list[EventAssignment]
    distance: str
    type: EventType
    initial_satisfied: bool = True
    priority: Optional[str] = None
