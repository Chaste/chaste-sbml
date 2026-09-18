# Changelog

Notable changes to chaste-sbml are documented here.

## [Unreleased]

## [0.0.1] - 2026-09-17

First release: a `chaste-sbml` command that turns an SBML model into Chaste C++
classes, plus the base classes those classes need.

### Added

- **Code generation from SBML.** Reads an SBML model and emits a Chaste ODE
  system (`<Model>SbmlOdeSystem.hpp` / `.cpp`).
- **Model types** (`--model-type`): `generic` for an ODE system alone, `srn` for
  a subcellular reaction network wrapper, and `cell-cycle` for a cell-cycle
  model whose division is driven by an SBML event.
- **Base classes** (`--copy-base-classes`): copies the C++ classes the generated
  code subclasses (`AbstractSbmlOdeSystem`, `AbstractSbmlSrnModel`,
  `AbstractSbmlCellCycleModel`, `SbmlMath.hpp` and friends) into a project.
- **Placeholder tests** (`--tests`, `--test-output-dir`): an optional CxxTest
  skeleton for the generated classes.
- **SBML support.** Levels 2 and 3; hierarchical models (`comp`) by flattening;
  species in amount and concentration; compartments; global and local
  parameters; reactions and kinetic laws; stoichiometry; function definitions;
  initial assignments; assignment and rate rules; events, including trigger
  initial values and priorities; and standard MathML mapped to `std::` or to the
  `sbmlmath` helpers.
- **Time handling** (`--timescale`): the model's native time unit is
  auto-detected and derivatives are scaled into Chaste's hours, with an override
  for models that declare no usable unit.

### Known limitations

- Unsupported constructs stop generation with an error rather than producing
  quietly wrong code: flux balance constraints (`fbc`), fast reactions,
  algebraic rules, non-zero event delays, the `delay()` MathML function, and a
  local parameter without a value.
- Units are not enforced: numeric unit annotations are stripped and quantities
  are treated as non-dimensional, so a model must be unit-consistent on its own.
- Of the 1821 semantic cases in the SBML Test Suite, 1345 pass (74%) and 476
  exercise unsupported features — `delay()` (262), algebraic rules (125), fast
  reactions (34), flux balance (34), random event execution (11), and event
  execution semantics (10).

[Unreleased]: https://github.com/Chaste/chaste-sbml/compare/v0.0.1...HEAD
[0.0.1]: https://github.com/Chaste/chaste-sbml/releases/tag/v0.0.1
