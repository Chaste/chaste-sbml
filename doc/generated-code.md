# Anatomy of generated code

This page describes the C++ files `chaste-sbml` emits and how they fit together.
The examples come from the reference model `Goldbeter1991`, generated as an
`srn` model.

## What gets generated

The files produced depend on `--model-type`. All names are derived from the
model name (the SBML filename with a `Sbml` suffix, unless you pass one):

```{list-table}
:header-rows: 1
:widths: 20 40 40

* - Model type
  - Files
  - Contents
* - `generic`
  - `<Name>OdeSystem.hpp` / `.cpp`
  - The ODE system only.
* - `srn`
  - `<Name>OdeSystem.hpp` / `.cpp`,
    `<Name>SrnModel.hpp` / `.cpp`
  - The ODE system plus an SRN model that wraps it.
* - `cell-cycle`
  - `<Name>OdeSystem.hpp` / `.cpp`,
    `<Name>CellCycleModel.hpp` / `.cpp`
  - The ODE system plus a cell-cycle model that wraps it.
```

Unless you pass `--no-tests`, a placeholder test `Test<Name>.hpp` is emitted as
well. All output is formatted with `clang-format`.

## The class hierarchy

Generated classes extend the [shipped base classes](#the-base-classes), which in
turn extend Chaste's core classes:

```text
AbstractOdeSystem              (Chaste)
└── AbstractSbmlOdeSystem      (shipped)
    └── <Name>SbmlOdeSystem    (generated)

AbstractOdeSrnModel            (Chaste)
└── AbstractSbmlSrnModel       (shipped)
    └── <Name>SbmlSrnModel     (generated, srn)

AbstractOdeBasedCellCycleModel (Chaste)
└── AbstractSbmlCellCycleModel (shipped)
    └── <Name>SbmlCellCycleModel (generated, cell-cycle)
```

## The ODE system

`<Name>SbmlOdeSystem` is the heart of the output: it turns the SBML model into a
Chaste `AbstractOdeSystem`. The header groups the model's quantities into
labelled `double` members:

```cpp
class Goldbeter1991SbmlOdeSystem : public AbstractSbmlOdeSystem
{
    // PARAMETERS
    double VM1; // VM1
    // ...
    // STATE VARIABLES
    double C; // Cyclin
    double M; // cdc_2_kinase
    double X; // Cyclin Protease

    double d_C_dt;   // one derivative member per state variable
    // ...
    // DERIVED QUANTITIES
    double cell; // cell
    // REACTIONS
    double reaction1; // creation of cyclin
    // ...
};
```

The constructor tells the base class how many state variables, parameters, and
events the model has, then initialises the system:

```cpp
Goldbeter1991SbmlOdeSystem::Goldbeter1991SbmlOdeSystem()
        : AbstractSbmlOdeSystem(3, 3, 0)   // (numStateVars, numParameters, numEvents)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>);
    Initialise();
    // EVENTS
}
```

Key methods:

`Initialise(double time = 0.0)`
: Sets every parameter, derived quantity, and reaction to its initial value, and
  seeds the state variables and their default initial conditions.

`RunModelEquations(time, rStateVariables)`
: The computational core. It reads the state, recomputes every intermediate
  equation and reaction, and returns the vector of derivatives.

`EvaluateYDerivatives(time, rY, rDY)` *(override)*
: Called by the ODE solver. It runs the model equations and writes the
  derivatives into `rDY`. When time-unit scaling applies, it multiplies time and
  each derivative by a file-local `TIMESCALE_MULTIPLIER`:

  ```cpp
  constexpr double TIMESCALE_MULTIPLIER = 3600.0; // seconds -> hours
  // ...
  rDY[i] = TIMESCALE_MULTIPLIER * derivatives[i];
  ```

`ComputeDerivedQuantities(time, rY)` *(override)*
: Returns the derived quantities, reaction fluxes, and amount/concentration
  conversions, for post-processing and output.

`ProcessModelEvents(time, rY)` *(override, private)*
: Present for every model, but only does work when the model has events. For each
  event it computes how close the trigger is to firing (used by CVODE
  root-finding) and records the deferred event assignments. See
  [Events](#events) below.

A `CellwiseOdeSystemInformation` specialization at the bottom of the `.cpp`
records variable names, units, and initial conditions.

## The wrapper (SRN / cell-cycle)

For the `srn` and `cell-cycle` types, the wrapper class owns an instance of the
ODE system and plugs it into Chaste's cell-based machinery. The two variants
share one template and differ only where the base classes genuinely differ.

The default constructor forwards to the base class; note the **SRN** wrapper
passes the state-variable count and the cell-cycle wrapper does not:

```cpp
Goldbeter1991SbmlSrnModel::Goldbeter1991SbmlSrnModel(
        boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlSrnModel(3, pOdeSolver)   // cell-cycle: AbstractSbmlCellCycleModel(pOdeSolver)
{
}
```

`Initialise()` *(override)*
: Creates the ODE system and hands it to the base class. SRN:
  `AbstractSbmlSrnModel::Initialise(new <Name>SbmlOdeSystem);` — cell-cycle:
  `mpOdeSystem = new <Name>SbmlOdeSystem; AbstractSbmlCellCycleModel::Initialise();`

`CreateSrnModel()` / `CreateCellCycleModel()` *(override)*
: Called by Chaste at cell division; returns a copy via the protected
  copy-constructor, which deep-copies the ODE system for the daughter cell:

  ```cpp
  auto* p_ode_system = static_cast<Goldbeter1991SbmlOdeSystem*>(rModel.GetOdeSystem());
  this->SetOdeSystem(new Goldbeter1991SbmlOdeSystem(*p_ode_system));
  ```

`OutputSrnModelParameters()` / `OutputCellCycleModelParameters()` *(override)*
: Chains to the base class.

The solver is chosen implicitly: the wrapper base-class constructor calls
`sbmlodesolversetup::SetUpDefaultOdeSolver`, which prefers CVODE (with sensible
tolerances) and falls back to `RungeKutta4IvpOdeSolver`.

## Events

When the model has events, the ODE system's constructor sets up per-event state
(the event type and whether its trigger is satisfied at `t = 0`), and
`ProcessModelEvents` computes, for each event:

- a signed *trigger distance* (used by CVODE to root-find the exact firing time),
- whether the trigger is currently active, and
- the deferred state/parameter assignments to apply when it fires.

For a cell-cycle model, an event tagged `SbmlEventType::CELL_DIVISION` is what
drives division: the base class's `ReadyToDivide()` returns `true` only once such
an event has fired.

## The placeholder test

`Test<Name>.hpp` is a CxxTest skeleton, not a real test — it exists so the
generated code compiles and links under Chaste's test runner. It always contains
a smoke test for the ODE system:

```cpp
class TestGoldbeter1991SbmlOdeSystem : public CxxTest::TestSuite
{
public:
    void TestOdeSystem()
    {
        TS_ASSERT_THROWS_NOTHING(Goldbeter1991SbmlOdeSystem ode_system);
    }
};
```

For `srn` and `cell-cycle` models it adds a second suite, deriving from
`AbstractCellBasedTestSuite`, that constructs the wrapper. Both suites carry a
`// TODO: Add tests` marker — fill them in with the assertions your model needs.

## The base classes

The generated code `#include`s and subclasses a small, stable set of C++ base
classes that ship with the package. Copy them into your project with
[`--copy-base-classes`](command-line.md#copying-the-base-classes):

```{list-table}
:header-rows: 1
:widths: 35 65

* - File
  - Role
* - `AbstractSbmlOdeSystem.hpp` / `.cpp`
  - Base for every generated ODE system; extends Chaste's `AbstractOdeSystem`
    and implements the event root-finding and clamping machinery.
* - `AbstractSbmlSrnModel.hpp` / `.cpp`
  - Base for SRN wrappers; extends `AbstractOdeSrnModel`.
* - `AbstractSbmlCellCycleModel.hpp` / `.cpp`
  - Base for cell-cycle wrappers; extends `AbstractOdeBasedCellCycleModel` and
    ties SBML cell-division events to `ReadyToDivide()`.
* - `SbmlEventType.hpp`
  - `enum class SbmlEventType { CELL_DIVISION, UNKNOWN }`, tagging which events
    cause division.
* - `SbmlMath.hpp`
  - The `sbmlmath` namespace (aliased `sm` in generated code) of helper functions
    for MathML operators C++ lacks (`root`, `piecewise`, `min`/`max`, and so on).
* - `SbmlOdeSolverSetup.hpp`
  - `SetUpDefaultOdeSolver`: selects CVODE where available, else Runge-Kutta 4.
```

:::{important}
The base classes are versioned with the tool. Re-run `--copy-base-classes`
whenever you upgrade `chaste-sbml`, so the copied classes match the code being
generated.
:::

:::{seealso}
[Using generated code in Chaste](using-in-chaste.md) walks through wiring these
files into a Chaste user project.
:::
