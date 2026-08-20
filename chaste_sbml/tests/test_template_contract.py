"""Static lint: a template may only reference variables its render context provides.

StrictUndefined (see chaste_sbml/_rendering.py) turns an unprovided template variable into an
error -- but only for the branches that actually execute during a given render. A variable
referenced solely in, say, the cell-cycle arm of test.hpp stays invisible until a cell-cycle model
is rendered. This lint parses every template statically and checks every variable it references,
across all branches, against the context its render path provides -- so a missing variable is
caught the moment a template is edited, without needing a model that reaches that branch.

The BASE/SRN/CELL_CYCLE key sets are anchored to what a real model actually provides by
test_context_contract_matches_model, so they cannot silently drift from _populate_template_vars;
that anchor uses ``==``, so it also flags a context key that is provided but referenced by no
template (a dead key).
"""

import pytest
from jinja2 import meta

from chaste_sbml import ChasteSbmlModel
from chaste_sbml._config import ROOT_DIR, ModelType
from chaste_sbml._rendering import CodeRenderer

REFERENCE_DIR = ROOT_DIR / "SbmlRefModels" / "src" / "reference"

# Jinja's built-ins plus the globals registered in _rendering.py: always available to every template.
GLOBALS = set(CodeRenderer._env.globals)

# Context keys ChasteSbmlModel._populate_template_vars provides for every model (the last block are
# the collections from ModelBuilder.template_data()).
BASE_KEYS = {
    "model_name",
    "ode_class_name",
    "ode_header_guard",
    "ode_hpp_file",
    "test_header_guard",
    "model_type",
    "scale_time",
    "time_multiplier",
    "time_unit_display",
    "derived_quantities",
    "equations",
    "events",
    "functions",
    "parameters",
    "reactions",
    "state_variables",
}
# Added on top of BASE_KEYS for any model with an SRN/cell-cycle wrapper (see _WRAPPER_SPECS).
WRAPPER_KEYS = {
    "wrapper_class_name",
    "wrapper_header_guard",
    "wrapper_hpp_file",
    "wrapper_base_class",
    "wrapper_abstract_class",
    "wrapper_create_method",
    "wrapper_output_method",
    "wrapper_model_noun",
}
# Added by ChasteSbmlTestSuiteModel for semantic cases (see SbmlRefModels/generate_cases.py).
CASE_KEYS = {
    "test_result_columns",
    "test_result_data",
    "test_amounts",
    "test_concentrations",
    "test_settings",
}

# The exact context each packaged template is rendered with (the union of every context it may see;
# test.hpp is rendered for all model types, so it may reference the wrapper keys).
TEMPLATE_CONTEXT = {
    "ode/ode.cpp": BASE_KEYS,
    "ode/ode.hpp": BASE_KEYS,
    "wrapper/wrapper.cpp": BASE_KEYS | WRAPPER_KEYS,
    "wrapper/wrapper.hpp": BASE_KEYS | WRAPPER_KEYS,
    "test/test.hpp": BASE_KEYS | WRAPPER_KEYS,
    "cases/semantic.hpp": BASE_KEYS | CASE_KEYS,
}


def _referenced_variables(template_name: str) -> set:
    """Return every top-level variable a template references, across all branches."""
    env = CodeRenderer._env
    source = env.loader.get_source(env, template_name)[0]
    return meta.find_undeclared_variables(env.parse(source))


def test_all_templates_are_linted():
    """Every packaged template is covered by TEMPLATE_CONTEXT, so a new one can't skip the lint."""
    assert set(CodeRenderer._env.list_templates()) == set(TEMPLATE_CONTEXT)


@pytest.mark.parametrize("template_name", sorted(TEMPLATE_CONTEXT))
def test_template_references_only_provided_variables(template_name):
    """A template only references variables its render context (or a global) provides."""
    allowed = GLOBALS | TEMPLATE_CONTEXT[template_name]
    unprovided = _referenced_variables(template_name) - allowed
    assert not unprovided, f"{template_name} references unprovided variables: {sorted(unprovided)}"


def test_cpp_metacharacters_are_not_escaped():
    """Autoescaping is off: C++ operators (< > & " ') must render verbatim, not HTML-escaped.

    These templates emit C++, so escaping would corrupt operators like ``<``, ``>`` and ``&&``.
    Rendering an expression through the real environment guards the property behaviourally,
    regardless of how autoescape happens to be configured.
    """
    expr = "if (a < b && c > d) x = &y; s = \"q'\";"
    rendered = CodeRenderer._env.from_string("{{ x }}").render(x=expr)
    assert rendered == expr


def test_none_output_is_rejected():
    """The finalize hook turns a provided-yet-None value into an error, not the text ``None``.

    StrictUndefined only catches an *unprovided* variable; a value that is present but ``None``
    (e.g. an unset ``initial_value`` rendered without its ``is not none`` guard) would otherwise
    emit the literal ``None`` into the generated C++. Rendering one through the real environment
    checks the hook is wired up and raises instead.
    """
    with pytest.raises(ValueError, match="None value"):
        CodeRenderer._env.from_string("{{ x }}").render(x=None)


@pytest.mark.parametrize(
    ("model_name", "model_type", "expected_keys"),
    [
        ("Chen2000", ModelType.GENERIC, BASE_KEYS),
        ("Goldbeter1991", ModelType.SRN, BASE_KEYS | WRAPPER_KEYS),
        ("Chen2000", ModelType.CELL_CYCLE, BASE_KEYS | WRAPPER_KEYS),
    ],
)
def test_context_contract_matches_model(model_name, model_type, expected_keys):
    """The BASE/WRAPPER contract matches what a real model provides, so it cannot drift
    from _populate_template_vars, and any provided-but-unreferenced (dead) key is flagged."""
    sbml_file = str(REFERENCE_DIR / model_name / f"{model_name}.xml")
    model = ChasteSbmlModel(sbml_file, model_type=model_type)
    assert set(model._template_vars) == expected_keys
