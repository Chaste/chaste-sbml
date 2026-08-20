# Using generated code in Chaste

The generated classes are ordinary Chaste C++, designed to live in a
[Chaste user project](https://chaste.github.io/docs/user-guides/user-projects/).
This page shows how to lay one out.

## Project layout

A Chaste user project sits under `projects/` in your Chaste source tree (or is
symlinked there) and has a `src/` directory for code and a `test/` directory for
tests:

```text
Chaste/projects/MyProject/
├── CMakeLists.txt
├── src/
│   ├── AbstractSbmlOdeSystem.hpp/.cpp     # copied base classes
│   ├── AbstractSbmlSrnModel.hpp/.cpp
│   ├── SbmlMath.hpp
│   ├── ...                                # (the rest of the base classes)
│   ├── MyModelSbmlOdeSystem.hpp/.cpp      # generated
│   └── MyModelSbmlSrnModel.hpp/.cpp       # generated
└── test/
    ├── CMakeLists.txt
    ├── ContinuousTestPack.txt
    └── TestMyModelSbml.hpp                # generated placeholder test
```

The two `CMakeLists.txt` files are the standard Chaste user-project boilerplate:

```cmake
# CMakeLists.txt
find_package(Chaste COMPONENTS cell_based)
chaste_do_project(MyProject)
```

```cmake
# test/CMakeLists.txt
chaste_do_test_project(MyProject)
```

## Step 1 — copy the base classes

Copy the shipped base classes into your project's `src/` so they match the
version of the tool you generated with:

```bash
chaste-sbml --copy-base-classes --output-dir Chaste/projects/MyProject/src
```

Do this once per project, and again whenever you upgrade `chaste-sbml`. See
[the base classes](generated-code.md#the-base-classes) for what gets copied.

## Step 2 — generate the model

Generate the model into `src/`, and send the placeholder test to `test/`:

```bash
chaste-sbml my_model.xml \
  --model-type srn \
  --output-dir Chaste/projects/MyProject/src \
  --test-output-dir Chaste/projects/MyProject/test
```

## Step 3 — register the test

Chaste discovers tests through *test packs* — text files in `test/` listing the
test headers to build. Add the generated placeholder to one:

```text
# test/ContinuousTestPack.txt
TestMyModelSbml.hpp
```

## Step 4 — build and run

Configure and build from your Chaste build directory as usual:

```bash
cd Chaste/build
cmake ..
cmake --build . --target project_MyProject
ctest -R MyProject
```

## Using the classes in a simulation

Once compiled, the generated classes are used like any other Chaste ODE system,
SRN model, or cell-cycle model.

For a **generic** model, construct the ODE system and hand it to a Chaste ODE
solver:

```cpp
#include "MyModelSbmlOdeSystem.hpp"

MyModelSbmlOdeSystem ode_system;
// solve with any AbstractIvpOdeSolver, then read state / derived quantities
```

For an **SRN** model, attach it to a cell so it runs as that cell's subcellular
reaction network:

```cpp
#include "MyModelSbmlSrnModel.hpp"

MAKE_PTR(WildTypeCellMutationState, p_state);
MAKE_PTR(StemCellProliferativeType, p_type);
CellPtr p_cell(new Cell(p_state, new MyModelSbmlSrnModel()));
p_cell->SetCellProliferativeType(p_type);
```

For a **cell-cycle** model, pass it as the cell's cell-cycle model; division is
driven by the SBML cell-division event:

```cpp
#include "MyModelSbmlCellCycleModel.hpp"

CellPtr p_cell(new Cell(p_state, new MyModelSbmlCellCycleModel()));
```

:::{seealso}
- [Anatomy of generated code](generated-code.md) — what each class and method
  does.
- The [Goldbeter 1991 tutorial](tutorials/goldbeter1991.md) — an end-to-end SRN
  example.
- The [Chaste user-project documentation](https://chaste.github.io) — for
  project setup details not specific to SBML.
:::
