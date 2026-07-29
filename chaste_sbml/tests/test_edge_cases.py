"""Edge-case tests exercising codegen branches not reached by the reference/SBML-suite models.

Each test writes a minimal SBML document and drives the builder, targeting a specific branch in
``_model_builder`` / ``chaste_sbml_model`` that the broader corpus does not cover.
"""

import pytest

from chaste_sbml import ChasteSbmlModel

L3_HEADER = '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">'


def _write(tmp_path, name, body, model_attrs=""):
    """Write an SBML L3v2 document with the given <model> body and return its path."""
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f"{L3_HEADER}\n"
        f'  <model id="{name}" {model_attrs}>\n'
        f"{body}\n"
        "  </model>\n"
        "</sbml>\n"
    )
    path = tmp_path / f"{name}.xml"
    path.write_text(xml)
    return str(path)


def test_missing_file_raises(tmp_path):
    """Constructing with a nonexistent SBML path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        ChasteSbmlModel(str(tmp_path / "does_not_exist.xml"))


def test_reaction_without_kinetic_law(tmp_path):
    """A reaction with no kinetic law builds its rate from the reactant species."""
    body = """    <listOfCompartments>
      <compartment id="c" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="P" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="R" reversible="false">
        <listOfReactants>
          <speciesReference species="S" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="P" stoichiometry="1" constant="true"/>
        </listOfProducts>
      </reaction>
    </listOfReactions>"""
    out = tmp_path / "out"
    out.mkdir()
    ChasteSbmlModel(_write(tmp_path, "no_kinetic_law", body)).write(str(out))
    assert any(out.iterdir())


def test_local_parameter_without_value_raises(tmp_path):
    """A local parameter with no value is a hard error."""
    body = """    <listOfCompartments>
      <compartment id="c" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="P" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="R" reversible="false">
        <listOfReactants>
          <speciesReference species="S" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="P" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>kloc</ci><ci>S</ci></apply>
          </math>
          <listOfLocalParameters>
            <localParameter id="kloc"/>
          </listOfLocalParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>"""
    with pytest.raises(ValueError, match="has no value"):
        ChasteSbmlModel(_write(tmp_path, "no_local_value", body))


def _build(tmp_path, name, body):
    """Build and write a model, asserting it produced output."""
    out = tmp_path / f"out_{name}"
    out.mkdir()
    ChasteSbmlModel(_write(tmp_path, name, body)).write(str(out))
    assert any(out.iterdir())


def test_time_varying_compartment_conversion_factor(tmp_path):
    """Species with a conversion factor in a time-varying compartment get a scaled dilution ODE.

    ``c`` is sized by the assignment rule ``c = X`` (time-varying via X's reaction), so both a
    boundary species and a plain species in ``c`` gain a dilution derivative, each multiplied by its
    conversion factor.
    """
    body = """    <listOfCompartments>
      <compartment id="c0" size="1" constant="true"/>
      <compartment id="c" size="1" constant="false"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="X" compartment="c0" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="P" compartment="c0" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="Sb" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="true" constant="false" conversionFactor="cf"/>
      <species id="Sd" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" conversionFactor="cf"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="1" constant="true"/>
      <parameter id="cf" value="2" constant="true"/>
    </listOfParameters>
    <listOfRules>
      <assignmentRule variable="c">
        <math xmlns="http://www.w3.org/1998/Math/MathML"><ci>X</ci></math>
      </assignmentRule>
    </listOfRules>
    <listOfReactions>
      <reaction id="R" reversible="false">
        <listOfReactants>
          <speciesReference species="X" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="P" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>k</ci><ci>X</ci></apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>"""
    _build(tmp_path, "cf_dilution", body)


def test_local_parameter_shadows_variable(tmp_path):
    """A reaction whose local parameter shares a species' name does not depend on that species.

    Reaction ``R2`` has a local parameter ``A`` used in its rate, shadowing the species ``A``; the
    dependency sort must not treat the reaction as depending on species ``A``.
    """
    body = """    <listOfCompartments>
      <compartment id="c" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="A" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="B" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="1" constant="true"/>
    </listOfParameters>
    <listOfReactions>
      <reaction id="R1" reversible="false">
        <listOfReactants>
          <speciesReference species="A" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="B" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>k</ci><ci>A</ci></apply>
          </math>
        </kineticLaw>
      </reaction>
      <reaction id="R2" reversible="false">
        <listOfReactants>
          <speciesReference species="B" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="A" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>A</ci><ci>B</ci></apply>
          </math>
          <listOfLocalParameters>
            <localParameter id="A" value="1"/>
          </listOfLocalParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>"""
    _build(tmp_path, "local_shadow", body)


def test_non_differentiable_compartment_rule(tmp_path):
    """A compartment sized by a non-differentiable rule skips that term in the chain rule.

    ``c = ceil(X)`` has no analytic derivative w.r.t. X, so the total-time-derivative skips it.
    """
    body = """    <listOfCompartments>
      <compartment id="c0" size="1" constant="true"/>
      <compartment id="c" size="1" constant="false"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="X" compartment="c0" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="P" compartment="c0" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="Sd" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="1" constant="true"/>
    </listOfParameters>
    <listOfRules>
      <assignmentRule variable="c">
        <math xmlns="http://www.w3.org/1998/Math/MathML"><apply><ceiling/><ci>X</ci></apply></math>
      </assignmentRule>
    </listOfRules>
    <listOfReactions>
      <reaction id="R" reversible="false">
        <listOfReactants>
          <speciesReference species="X" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="P" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>k</ci><ci>X</ci></apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>"""
    _build(tmp_path, "nondiff", body)
