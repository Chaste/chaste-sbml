"""C++ identifier naming for generated model code.

Provides the pieces that keep generated C++ identifiers valid and unique:

* ``find_name_conflicts``: detect clashes and fail loudly, so the generator
  never silently emits C++ in which two entities share a name, or a name is a
  C++ keyword, a Chaste base-class member or an invalid identifier.
* ``unique_name``: allocate a collision-free name for a generator-synthesised
  identifier.
* ``resolve_cpp_name``: make a real SBML id safe to emit verbatim by escaping
  keyword and reserved-name clashes.

It also derives C++ names from arbitrary strings, using ``to_cpp_name`` and
``to_camel_case`` for the model class name, and ``generate_header_guard`` for
include guards.

SBML already guarantees its own object ids are unique, so in practice conflicts
arise from generator-synthesised names (e.g. amount/concentration conversions
``amt__x``/``conc__x``, and state derivatives ``d_x_dt``) colliding with a real
id of the same spelling, or from an id that happens to be a C++ keyword or a
name inherited from the base classes.
"""

import re

# C++ keywords (and alternative operator spellings); none may be used as an identifier.
CPP_KEYWORDS = frozenset(
    {
        "alignas",
        "alignof",
        "and",
        "and_eq",
        "asm",
        "auto",
        "bitand",
        "bitor",
        "bool",
        "break",
        "case",
        "catch",
        "char",
        "char16_t",
        "char32_t",
        "class",
        "compl",
        "concept",
        "const",
        "consteval",
        "constexpr",
        "constinit",
        "const_cast",
        "continue",
        "co_await",
        "co_return",
        "co_yield",
        "decltype",
        "default",
        "delete",
        "do",
        "double",
        "dynamic_cast",
        "else",
        "enum",
        "explicit",
        "export",
        "extern",
        "false",
        "float",
        "for",
        "friend",
        "goto",
        "if",
        "inline",
        "int",
        "long",
        "mutable",
        "namespace",
        "new",
        "noexcept",
        "not",
        "not_eq",
        "nullptr",
        "operator",
        "or",
        "or_eq",
        "private",
        "protected",
        "public",
        "register",
        "reinterpret_cast",
        "requires",
        "return",
        "short",
        "signed",
        "sizeof",
        "static",
        "static_assert",
        "static_cast",
        "struct",
        "switch",
        "template",
        "this",
        "thread_local",
        "throw",
        "true",
        "try",
        "typedef",
        "typeid",
        "typename",
        "union",
        "unsigned",
        "using",
        "virtual",
        "void",
        "volatile",
        "wchar_t",
        "while",
        "xor",
        "xor_eq",
    }
)

# Identifiers inherited from the Chaste base classes (AbstractSbmlOdeSystem and its bases)
# that a generated subclass must not shadow with a member/local of the same name.
CHASTE_RESERVED_NAMES = frozenset(
    {
        # Data members
        "mNumberOfParameters",
        "mNumberOfEvents",
        "mNumberOfStateVariables",
        "mStateVariables",
        "mParameters",
        "mpSystemInfo",
        "mEventSatisfied",
        "mEventClampActive",
        "mEventTriggered",
        "mEventType",
        "mEventAdjustedStateVars",
        "mEventAdjustedStateValues",
        "mEventAdjustedStatePriority",
        "mEventAdjustedParameters",
        "mEventAdjustedParameterValues",
        "mEventAdjustedParameterPriority",
        # Frequently-referenced methods
        "GetParameter",
        "SetParameter",
        "GetStateVariable",
        "SetStateVariable",
        "ComputeDerivedQuantities",
        "EvaluateYDerivatives",
        "Initialise",
        "RunModelEquations",
        "ProcessModelEvents",
        # The SBML time symbol is emitted as a bare `time` identifier in equations.
        "time",
    }
)

# C++ macros (from <cmath>/<cstdlib>) that expand when used as a bare identifier, so an SBML id
# spelled like one cannot be emitted verbatim (e.g. a parameter NAN becomes `double NAN;`, which
# the preprocessor rewrites into the macro's expansion).
CPP_MACROS = frozenset({"NAN", "INFINITY", "HUGE_VAL", "EOF", "NULL"})

# Every name the generator may not use for one of its own identifiers.
RESERVED_NAMES = CPP_KEYWORDS | CHASTE_RESERVED_NAMES | CPP_MACROS

_IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class NameConflictError(ValueError):
    """Raised when generated C++ identifiers clash with each other or with reserved names."""


def is_valid_cpp_identifier(name: str) -> bool:
    """Return True if ``name`` is a syntactically valid C++ identifier.

    This does not reject keywords (which are valid *syntax* but reserved); use
    :data:`CPP_KEYWORDS` for that.

    :param name: The candidate identifier.
    :return: True if ``name`` matches ``[A-Za-z_][A-Za-z0-9_]*``.
    """
    return bool(_IDENTIFIER_RE.match(name))


def find_name_conflicts(identifiers, reserved=RESERVED_NAMES) -> list[str]:
    """Find C++ identifier conflicts among the names the generator will emit.

    :param identifiers: an iterable of ``(name, kind)`` pairs, one per C++ identifier the
        generator will emit; ``kind`` is a short human-readable category (e.g. ``"parameter"``)
        used only in the returned messages.
    :param reserved: names that must not be used (C++ keywords and Chaste base-class members).
    :return: a sorted list of human-readable conflict messages; empty if there are none.
    """
    messages = []
    first_seen: dict[str, str] = {}

    for name, kind in identifiers:
        if not is_valid_cpp_identifier(name):
            messages.append(f"{kind} '{name}' is not a valid C++ identifier")
        elif name in CPP_KEYWORDS:
            messages.append(f"{kind} '{name}' is a C++ keyword")
        elif name in reserved:
            messages.append(f"{kind} '{name}' clashes with a reserved Chaste name")

        if name in first_seen:
            messages.append(f"'{name}' is emitted twice: as {first_seen[name]} and as {kind}")
        else:
            first_seen[name] = kind

    return sorted(messages)


def unique_name(base: str, taken) -> str:
    """Return a name based on ``base`` that is not in ``taken``.

    Used to keep generator-synthesised identifiers (amount/concentration conversions, state
    derivatives, initial-assignment intermediates) collision-free. ``base`` is returned
    unchanged when it is free, so names stay clean unless a real id genuinely occupies them;
    otherwise the smallest ``_N`` (N >= 2) suffix that is free is appended.

    :param base: The desired identifier.
    :param taken: A container of names already in use.
    :return: ``base`` if free, else ``base`` with the smallest free ``_N`` suffix.
    """
    if base not in taken:
        return base
    n = 2
    while f"{base}_{n}" in taken:
        n += 1
    return f"{base}_{n}"


def resolve_cpp_name(base: str, taken) -> str:
    """Return a C++-safe, unique replacement for a real SBML id.

    An SBML ``SId`` already matches C++ identifier syntax, so the only clashes to fix are with
    C++ keywords and reserved Chaste names: those are escaped by appending ``_`` until the name
    is neither. Uniqueness against ``taken`` is then ensured. A name that needs no change is
    returned unchanged.

    :param base: The real SBML id to make safe.
    :param taken: A container of names already in use (other ids and reserved names).
    :return: A C++-safe identifier not present in ``taken``.
    """
    candidate = base
    while candidate in CPP_KEYWORDS or candidate in CHASTE_RESERVED_NAMES or candidate in CPP_MACROS:
        candidate += "_"
    return unique_name(candidate, taken)


class NameManager:
    """Manages C++ identifier naming for one SBML model.

    It resolves real ids that clash with C++ keywords or reserved names and hands out
    collision-free names for generator-synthesised identifiers such as state derivatives and
    amount/concentration conversions. The residual conflict *detection* stays
    with the model, which validates its own built collections.
    """

    def __init__(self, sbml_model) -> None:
        """:param sbml_model: The libsbml Model whose ids are being turned into C++."""
        self._sbml_model = sbml_model
        self._taken: set = set()
        # C++ identifier (the renamed id) -> original SBML id, for ids renamed to dodge a C++
        # keyword or reserved name. Lets the generated code report a variable under its real SBML
        # id even though its C++ identifier had to be escaped.
        self._sbml_names: dict[str, str] = {}

    def resolve_real_id_conflicts(self) -> None:
        """Rename real SBML ids that are C++ keywords or reserved Chaste names.

        An SBML ``SId`` is already valid C++ identifier syntax, so the only ids that cannot be
        emitted verbatim are those equal to a C++ keyword (e.g. a compartment called ``default``)
        or a Chaste base-class member. Each such id is renamed in place to a safe unique name and
        every reference to it is updated via libsbml's per-element ``renameSIdRefs``, so the rest
        of generation -- equations, events, initial assignments, the templates -- sees only clean
        ids. Ids that are already safe are left untouched, so conflict-free models are unchanged.
        Only global SId-namespace entities that become C++ identifiers are considered.
        """
        model = self._sbml_model
        # Replacements must avoid every existing id plus the reserved names.
        taken = self._collect_ids() | RESERVED_NAMES

        elements = []
        for lst in (
            model.getListOfSpecies(),
            model.getListOfParameters(),
            model.getListOfCompartments(),
            model.getListOfReactions(),
            model.getListOfFunctionDefinitions(),
        ):
            elements.extend(lst)
        for reaction in model.getListOfReactions():
            elements.extend(list(reaction.getListOfReactants()) + list(reaction.getListOfProducts()))

        for element in elements:
            if not element.isSetId():
                continue
            old_id = element.getId()
            if old_id not in CPP_KEYWORDS and old_id not in CHASTE_RESERVED_NAMES and old_id not in CPP_MACROS:
                continue
            new_id = resolve_cpp_name(old_id, taken)
            element.setId(new_id)
            for referrer in model.getListOfAllElements():
                referrer.renameSIdRefs(old_id, new_id)
            taken.add(new_id)
            self._sbml_names[new_id] = old_id

    def sbml_name(self, cpp_id: str) -> str:
        """Return the original SBML id for a (possibly renamed) C++ identifier.

        A variable is reported to Chaste under this name, so it is looked up by its real SBML id
        even when the emitted C++ identifier had to be escaped (e.g. a parameter ``time`` emitted
        as ``time_``). Ids that were never renamed are returned unchanged.

        :param cpp_id: The C++ identifier as it appears in the generated code.
        :return: The original SBML id if ``cpp_id`` was a renamed id, else ``cpp_id``.
        """
        return self._sbml_names.get(cpp_id, cpp_id)

    def reset(self) -> None:
        """Recompute the taken-name set from the (resolved) model, ready for a build.

        Call before allocating synthetic names, so they avoid every real id and reserved name.
        """
        self._taken = self._collect_ids() | RESERVED_NAMES

    def reserve(self, base: str) -> str:
        """Reserve a unique C++ identifier for a generator-synthesised variable.

        Returns ``base`` unchanged when it is free (the usual case, so names stay clean),
        otherwise the smallest ``_N`` suffix that avoids every real id, reserved name and
        previously-reserved synthetic. The chosen name is recorded so later synthetics avoid it.

        :param base: The desired synthetic identifier (e.g. ``d_C_dt`` or ``amt__X``).
        :return: A unique C++ identifier.
        """
        name = unique_name(base, self._taken)
        self._taken.add(name)
        return name

    def _collect_ids(self) -> set[str]:
        """Collect the real SBML ids that become C++ identifiers in the generated code.

        Covers the categories emitted as C++ names -- species, global and local (kinetic-law)
        parameters, compartments, reactions, function definitions and species-reference
        stoichiometry ids. Units, events and rules are excluded as they never become identifiers,
        so a synthetic name is only escaped when it truly collides.

        :return: The set of real ids that names must avoid.
        """
        model = self._sbml_model
        names = set()
        for lst in (
            model.getListOfSpecies(),
            model.getListOfParameters(),
            model.getListOfCompartments(),
            model.getListOfReactions(),
            model.getListOfFunctionDefinitions(),
        ):
            names.update(elem.getId() for elem in lst if elem.isSetId())
        for reaction in model.getListOfReactions():
            for ref in list(reaction.getListOfReactants()) + list(reaction.getListOfProducts()):
                if ref.isSetId():
                    names.add(ref.getId())
            kinetic_law = reaction.getKineticLaw()
            if kinetic_law is not None:
                names.update(lp.getId() for lp in kinetic_law.getListOfParameters() if lp.isSetId())
        return names


def to_cpp_name(name: str) -> str:
    """Sanitize an input name to a C++ compatible alphanumeric name.

    C++ keyword / reserved-name clashes are handled separately by :func:`resolve_cpp_name` and
    :class:`NameManager`; this only fixes the character set and leading character.

    :param name: The variable name.
    :return: The sanitized variable name.
    """
    name_ = name.strip()
    if not name_:
        return ""

    cpp_name = []
    # Prefix with "_" if name doesn't start with a letter or "_"
    if name_[0] != "_" and not name_[0].isalpha():
        cpp_name.append("_")

    for char in name_:
        if char.isalpha() or char.isdigit() or char == "_":
            cpp_name.append(char)
        else:
            # Replace other chars with "_"
            cpp_name.append("_")

    return "".join(cpp_name)


def to_camel_case(name: str) -> str:
    """Convert an input name to an alphanumeric name in camel case.

    :param name: The variable name.
    :return: The variable name in camel case.
    """
    camel = []

    caps = False  # True: capitalize next letter
    for char in name.strip():
        if char.isalpha():
            if caps:
                camel.append(char.upper())
                caps = False
            else:
                camel.append(char)
        elif char.isdigit():
            camel.append(char)
            caps = True
        else:
            # Skip other chars
            caps = True

    return "".join(camel)


def generate_header_guard(filename: str) -> str:
    """Generate a C++ header guard from a filename.

    :param filename: The filename.
    :return: The header guard.
    """
    name = to_cpp_name(filename)
    if not name:
        return ""

    prev = name[0]
    guard = [prev.upper()]

    for char in name[1:]:
        if char.isupper():
            # Add _ before uppercase chars if not sequence of uppercase chars
            if not prev.isupper():
                guard.append("_")
            guard.append(char)
        elif char.islower():
            guard.append(char.upper())
        elif char.isdigit():
            # Add _ between lowercase char and digit
            if prev.islower():
                guard.append("_")
            guard.append(char)
        elif guard[-1] != "_":
            # Replace non-alphanumeric chars with "_"; merging successive "_"s
            guard.append("_")
        prev = char

    return "".join(guard) + "_"
