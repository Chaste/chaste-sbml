# Tutorial: Goldbeter 1991

This tutorial walks through generating and understanding a complete model,
end to end. We use
[Goldbeter 1991](https://www.ebi.ac.uk/biomodels/BIOMD0000000003), a minimal
model of mitotic oscillations with three state variables: cyclin (`C`), active
cdc2 kinase (`M`), and active cyclin protease (`X`). The model is driven by seven reactions.
It is a relatively simple example containing SBML **reactions** and **assignment rules**,
but no **events**.

We generate it as an `srn` (subcellular reaction network) model, so it can run
inside a cell in a Chaste cell-based simulation.

## 1. Get the model

You can download `BIOMD0000000003` from [BioModels](https://www.ebi.ac.uk/biomodels/).
The model also ships with the repository under
`chaste_sbml/SbmlRefModels/src/reference/Goldbeter1991/Goldbeter1991.xml`.

## 2. Generate the code

Generate into a [Chaste user project](../using-in-chaste.md):

```bash
chaste-sbml chaste_sbml/SbmlRefModels/src/reference/Goldbeter1991/Goldbeter1991.xml --model-type srn \
  --output-dir Chaste/projects/MyProject/src \
  --test-output-dir Chaste/projects/MyProject/test
```

This produces the ODE system, the SRN wrapper, and a placeholder test:

```text
src/
├── Goldbeter1991SbmlOdeSystem.hpp/.cpp   # the ODE system
└── Goldbeter1991SbmlSrnModel.hpp/.cpp    # the SRN wrapper
test/
└── TestGoldbeter1991Sbml.hpp             # placeholder test (from --test-output-dir)
```

The model name `Goldbeter1991Sbml` is derived from the filename; the classes are
named after it.

:::{note}
The model declares no explicit time unit, but it is an SBML Level 2 model, so
seconds are assumed and derivatives are scaled to Chaste's hours by 3600. Pass
`--timescale` to override this — see [Time units](../command-line.md#time-units).
:::

## 3. Review the ODE system

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

The constructor declares the shape of the system to the base class i.e. three state
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

`RunModelEquations` recomputes the reactions and returns the derivatives.
`EvaluateYDerivatives`, the method the solver calls, applies the time scaling:

```cpp
constexpr double TIMESCALE_MULTIPLIER = 3600.0; // seconds -> hours
// ...
rDY[i] = TIMESCALE_MULTIPLIER * derivatives[i];
```

:::{seealso}
[Anatomy of generated code](../generated-code.md#the-ode-system): for the full
tour of these methods.
:::

## 4. Review the SRN wrapper

`Goldbeter1991SbmlSrnModel` wraps the ODE system so a cell can run it. Its
`Initialise()` creates the ODE system and hands it to the base class:

```cpp
void Goldbeter1991SbmlSrnModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    AbstractSbmlSrnModel::Initialise(new Goldbeter1991SbmlOdeSystem);
}
```

:::{note}
Chaste calls `CreateSrnModel()` at division to produce a daughter-cell copy,
deep-copying the ODE system.
:::

## 5. Build and run

The generated classes need the base classes alongside them, so copy those into
the project's `src/` too:

```bash
chaste-sbml --copy-base-classes --output-dir Chaste/projects/MyProject/src
```

Add `TestGoldbeter1991Sbml.hpp` to a test pack, then build and run:

```bash
cd Chaste/build
cmake .. && cmake --build . --target project_MyProject
ctest -R Goldbeter1991
```

The placeholder test only checks that the classes construct. Replace its
`// TODO: Add tests` with the assertions your model needs e.g.
integrating the ODE system and checking the oscillation period.

## Next steps

The [Tyson-Novak 2001 tutorial](tysonnovak2001.md) covers a more complex model
with events, a cell-division trigger, and a function definition.
