"""C++ identifier conflict detection for generated model code (issue #35).

Phase A: detect identifier conflicts and fail loudly, so the generator never silently
emits C++ in which two entities share a name or a name clashes with a C++ keyword or a
Chaste base-class member. Later phases will *resolve* such conflicts by renaming; this
module only diagnoses them.

SBML already guarantees its own object ids are unique, so in practice conflicts arise from
generator-synthesised names (amount/concentration conversions ``amt__x``/``conc__x``, state
derivatives ``d_x_dt``, initial assignments) colliding with a real id of the same spelling,
or from an id that happens to be a C++ keyword or a name inherited from the base classes.
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

# Every name the generator may not use for one of its own identifiers.
RESERVED_NAMES = CPP_KEYWORDS | CHASTE_RESERVED_NAMES

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
    while candidate in CPP_KEYWORDS or candidate in CHASTE_RESERVED_NAMES:
        candidate += "_"
    return unique_name(candidate, taken)
