# SBML feature support

This page summarises which SBML constructs `chaste-sbml` can translate
into Chaste C++, and which it cannot. Unsupported constructs are rejected with a
clear error during generation rather than producing wrong code.

## SBML levels and packages

- **SBML Level 2 and Level 3** are supported.
- **Hierarchical models** (the `comp` package) are supported by *flattening* the
  model on load. If flattening fails, generation stops with an error.
- **Flux balance constraints** (the `fbc` package) are **not supported**.
- Other packages are not interpreted.

:::{note}
A model whose SBML level is outside 2–3 is not rejected outright — the generator
warns and proceeds on a best-effort basis, so the output may be incorrect.
:::

## Model components

```{list-table}
:header-rows: 1
:widths: 30 20 50

* - Construct
  - Status
  - Notes
* - Species (amount & concentration)
  - Supported
  - `hasOnlySubstanceUnits`, `boundaryCondition`, and species/model
    `conversionFactor` are all handled; both amount and concentration are
    exposed.
* - Compartments
  - Supported
  - Constant, rule-driven, and ODE-driven compartments; time-varying
    compartments add dilution terms to conserve amount.
* - Global parameters
  - Supported
  - Constant, rule-varied, event-varied, or ODE-driven.
* - Local (kinetic-law) parameters
  - Supported
  - Scoped per reaction; shadow globals of the same name. A local parameter with
    no value is an error.
* - Reactions & kinetic laws
  - Supported
  - Each reaction is exposed as a flux; a reaction with no kinetic law falls back
    to the product of its reactants.
* - Stoichiometry
  - Supported
  - Constant, variable (via rate/assignment rule), and `stoichiometryMath`.
* - Function definitions
  - Supported
  - Emitted as C++ functions taking `double` arguments.
* - Initial assignments
  - Supported
  - Expanded on load where possible; the remainder applied to species,
    parameters, compartments, and stoichiometry.
```

## Rules

```{list-table}
:header-rows: 1
:widths: 30 20 50

* - Rule type
  - Status
  - Notes
* - Assignment rules
  - Supported
  - Applied to species, compartments, parameters, and stoichiometry.
* - Rate rules
  - Supported
  - Take precedence over a species' reaction contributions.
* - Algebraic rules
  - **Not supported**
  - Generation stops with an error.
```

## Events

Events are supported, with some limits:

- **Trigger**, **event assignments**, **priority**, and the trigger's
  **initial value** (persistent-at-`t=0`) semantics are handled.
- For CVODE root-finding, a *trigger distance* is derived from a single
  top-level relational operator with two operands. Compound boolean triggers or
  n-ary relationals fall back to a constant distance (they still fire, but not
  via smooth root-finding).
- **Event delays** are **not supported**: a non-zero (or non-numeric) delay is
  rejected. A zero delay is fine.

## Time and units

- **Time-unit scaling** to Chaste's hours is supported for time units that
  reduce to seconds (`second`, `millisecond`, `minute`, `hour`). See
  [Time units](command-line.md#time-units) for the resolution rules and the
  `--timescale` override.
- **General physical units are not enforced.** Numeric-literal unit annotations
  are stripped, and species/parameter/compartment quantities are treated as
  non-dimensional. Get the units right in your model.

## MathML

Standard MathML is translated to C++:

- Arithmetic, relational, and boolean operators.
- The full set of elementary functions (`sin`, `cos`, `exp`, `log`, `sqrt`,
  `pow`, `abs`, `floor`, `ceil`, the inverse and hyperbolic trig functions, and
  so on), mapped to `std::` or to `sm::` helper functions where C++ has no direct
  equivalent (for example `cot`, `factorial`, `piecewise`, `min`/`max`).
- Constants `pi`, `exponentiale`, `avogadro`, and the `time` and `rateOf`
  csymbols; `infinity` and `notanumber`.
- Calls to a model's own function definitions.

The **`delay()`** function is **not supported** and is rejected during
generation.

:::{note}
A model function may legitimately be *named* `sin`, `min`, `log`, and so on.
Those calls are left as calls to your function definition and are **not**
rerouted to the C++ built-in of the same name.
:::

## Explicitly unsupported constructs

Each of these stops generation with an error, so an unsupported model never
yields silently-wrong code:

```{list-table}
:header-rows: 1
:widths: 40 60

* - Construct
  - Behaviour
* - Flux balance constraints (`fbc` package)
  - Rejected — flux-balance modelling is out of scope.
* - Fast reactions
  - Rejected.
* - Algebraic rules
  - Rejected.
* - Event delays (non-zero)
  - Rejected.
* - `delay()` MathML function
  - Rejected.
* - Local parameter without a value
  - Rejected.
```

## SBML Test Suite results

The generator is run against the
[SBML Test Suite](https://github.com/sbmlteam/sbml-test-suite) semantic cases in
CI (see the `sbml_test_suite` workflow). Of **1821** semantic cases:

```{list-table}
:header-rows: 1
:widths: 60 20 20

* - Outcome
  - Cases
  - Share
* - **Pass** (generated and numerically verified against Chaste)
  - 1345
  - 74%
* - Unsupported (skipped)
  - 476
  - 26%
```

The unsupported cases break down by the feature they exercise:

```{list-table}
:header-rows: 1
:widths: 60 40

* - Feature
  - Cases
* - `delay()` / event delays
  - 262
* - Algebraic rules
  - 125
* - Flux balance / steady-state (no time course)
  - 34
* - Fast reactions
  - 34
* - Random event execution
  - 11
* - Event execution semantics
  - 10
```

:::{note}
These counts come from `chaste_sbml/SbmlRefModels/test/data/sbml_test_suite_status.csv`,
which CI checks against actual generation so the figures cannot drift out of
date.
:::
