# Tutorial: Tyson-Novak 2001

This tutorial walks through the [Tyson-Novak 2001](https://www.ebi.ac.uk/biomodels/BIOMD0000000195),
a model of the cell-cycle oscillator. It has eight state variables, twenty **reactions**,
**assignment rules**, a **function definition**, and **events**, one of
which represents cell division. This builds on the [Goldbeter 1991 tutorial](goldbeter1991.md)
with a richer model that exercises more complex SBML constructs.

We generate it as a `cell-cycle` model, so the division event drives when a
Chaste cell divides.

## 1. Generate the code

The model ships with the repository under
`chaste_sbml/SbmlRefModels/src/reference/TysonNovak2001/TysonNovak2001.xml`.

As before, generate into a [Chaste user project](../using-in-chaste.md):

```bash
chaste-sbml chaste_sbml/SbmlRefModels/src/reference/TysonNovak2001/TysonNovak2001.xml --model-type cell-cycle \
  --output-dir Chaste/projects/MyProject/src \
  --test-output-dir Chaste/projects/MyProject/test
```

This produces the ODE system, the cell-cycle wrapper, and a placeholder test:

```text
src/
├── TysonNovak2001SbmlOdeSystem.hpp/.cpp
└── TysonNovak2001SbmlCellCycleModel.hpp/.cpp
test/
└── TestTysonNovak2001Sbml.hpp
```

## 2. Function definitions

The model defines a Goldbeter-Koshland function `GK(...)`. SBML function
definitions become C++ member functions taking `double` arguments:

```cpp
// MODEL FUNCTIONS
inline double GK(double A1, double A2, double A3, double A4);
```

Calls to `GK` in the kinetic laws are emitted as ordinary calls to this method.

## 3. Events

The ODE constructor declares one event and sets up the per-event bookkeeping:

```cpp
TysonNovak2001SbmlOdeSystem::TysonNovak2001SbmlOdeSystem()
        : AbstractSbmlOdeSystem(8, 36, 1)   // 8 state vars, 36 parameters, 1 event
{
    // ...
    // EVENTS
    mEventType.resize(1, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division
    mEventType[0] = SbmlEventType::CELL_DIVISION; // Cell division

    mEventSatisfied = { true }; // From SBML trigger initialValue
    // ...
}
```

The generator tries to **guess** which events cause division from their labels
and marks them `SbmlEventType::CELL_DIVISION`.

:::{note}
Always review this block: uncomment or
comment the `mEventType[...]` lines so that exactly the events which should end
the cell cycle are tagged as division.
:::

The actual event logic lives in `ProcessModelEvents`. For each event it computes
a signed *trigger distance* so that CVODE can root-find the exact firing time,
evaluates the trigger, and records the deferred event assignments. In this case,
an event assignment halves the cell mass at division:

```cpp
double event_dist = (0.1) - (CycB) - std::numeric_limits<double>::epsilon();
bool active = CycB < 0.1;
// ...
mEventAdjustedStateValues[3] = m / 2.0;   // deferred: halve mass at division
```

:::{note}
Only a single top-level relational trigger (e.g. `CycB < 0.1`) yields a smooth
root-found firing time. Compound boolean or n-ary triggers still fire, but fall
back to a constant distance rather than smooth root-finding (e.g. `(CycB < 0.1) && (m > 2.0)`).
:::

## 4. Cell cycle

The cell-cycle wrapper `TysonNovak2001SbmlCellCycleModel` owns the ODE system.
Its base class, `AbstractSbmlCellCycleModel`, connects the SBML event to Chaste:
`ReadyToDivide()` returns `true` only once an event tagged
`SbmlEventType::CELL_DIVISION` has fired, and `ResetForDivision()` clears the
event state for the next cycle.

## 5. Build and run

Copy the base classes into the project's `src/`, register the test, and build. To
use it, pass the cell-cycle model when creating a cell:

```cpp
#include "TysonNovak2001SbmlCellCycleModel.hpp"

CellPtr p_cell(new Cell(p_state, new TysonNovak2001SbmlCellCycleModel()));
p_cell->SetCellProliferativeType(p_type);
```

The cell now advances its ODE system each timestep and divides when the SBML
division event fires.

## Things to check

- **Time units**: Confirm the assumed unit is right; pass `--timescale` if not.
- **Event tags**: Review the `mEventType[...]` lines as the cell division guess is
  heuristic.

:::{seealso}
[Supported SBML features](../sbml-features.md) for the full list of SBML
features the generator handles.
:::
