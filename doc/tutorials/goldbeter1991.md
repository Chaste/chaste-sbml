# Tutorial: Goldbeter 1991 (a simple SRN model)

This tutorial walks through generating and understanding a complete model,
end to end. We use
[Goldbeter 1991](https://www.ebi.ac.uk/biomodels/BIOMD0000000003), a minimal
model of mitotic oscillations with three state variables — cyclin (`C`), active
cdc2 kinase (`M`), and active cyclin protease (`X`) — driven by seven reactions.
It is a good first example: reactions and assignment rules, but no events or
delays.

We generate it as an `srn` (subcellular reaction network) model, so it can run
inside a cell in a Chaste cell-based simulation.

## 1. Get the model

The model ships with the repository under
`chaste_sbml/SbmlRefModels/src/reference/Goldbeter1991/Goldbeter1991.xml`, or you
can download `BIOMD0000000003` from [BioModels](https://www.ebi.ac.uk/biomodels/).

## 2. Generate the code

```bash
chaste-sbml Goldbeter1991.xml --model-type srn --output-dir out/
```

This produces four files in `out/`:

```text
Goldbeter1991SbmlOdeSystem.hpp / .cpp   # the ODE system
Goldbeter1991SbmlSrnModel.hpp  / .cpp   # the SRN wrapper
TestGoldbeter1991Sbml.hpp               # placeholder test (unless --no-tests)
```

The model name `Goldbeter1991Sbml` is derived from the filename; the classes are
named after it.

:::{note}
The model declares no explicit time unit, but it is an SBML Level 2 model, so
seconds are assumed and derivatives are scaled to Chaste's hours by 3600. Pass
`--timescale` to override this — see [Time units](../command-line.md#time-units).
:::

## 3. Read the ODE system

Open `Goldbeter1991SbmlOdeSystem.hpp`. The model's quantities appear as labelled
members, grouped by role:

```cpp
// STATE VARIABLES
double C; // Cyclin
double M; // cdc_2_kinase
double X; // Cyclin Protease
// ...
// REACTIONS
double reaction1; // creation of cyclin
```

The constructor declares the shape of the system to the base class — three state
variables, three parameters, and no events:

```cpp
Goldbeter1991SbmlOdeSystem::Goldbeter1991SbmlOdeSystem()
        : AbstractSbmlOdeSystem(3, 3, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991SbmlOdeSystem>);
    Initialise();
    // EVENTS   <- empty: this model has none
}
```

`RunModelEquations` recomputes the reactions and returns the derivatives, and
`EvaluateYDerivatives` — the method the solver calls — applies the time scaling:

```cpp
constexpr double TIMESCALE_MULTIPLIER = 3600.0; // seconds -> hours
// ...
rDY[i] = TIMESCALE_MULTIPLIER * derivatives[i];
```

See [Anatomy of generated code](../generated-code.md#the-ode-system) for the full
tour of these methods.

## 4. Read the SRN wrapper

`Goldbeter1991SbmlSrnModel` wraps the ODE system so a cell can run it. Its
`Initialise()` creates the ODE system and hands it to the base class:

```cpp
void Goldbeter1991SbmlSrnModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    AbstractSbmlSrnModel::Initialise(new Goldbeter1991SbmlOdeSystem);
}
```

`CreateSrnModel()` produces a daughter-cell copy at division, deep-copying the
ODE system. You don't call these yourself — Chaste does.

## 5. Build it in Chaste

Put the four files into a Chaste user project, copy the base classes alongside
them, and register the test. This is the standard flow described in
[Using generated code in Chaste](../using-in-chaste.md):

```bash
# in your project's src/
chaste-sbml --copy-base-classes --output-dir Chaste/projects/MyProject/src
chaste-sbml Goldbeter1991.xml --model-type srn \
  --output-dir Chaste/projects/MyProject/src \
  --test-output-dir Chaste/projects/MyProject/test
```

Add `TestGoldbeter1991Sbml.hpp` to a test pack, then build and run:

```bash
cd Chaste/build
cmake .. && cmake --build . --target project_MyProject
ctest -R Goldbeter1991
```

The placeholder test only checks that the classes construct. Replace its
`// TODO: Add tests` with the assertions your model needs — for instance,
integrating the ODE system and checking the oscillation period.

## Next steps

The [Tyson–Novak 2001 tutorial](tysonnovak2001.md) covers a more complex model
with events, a cell-division trigger, and a function definition.
