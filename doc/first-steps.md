# First steps

This page takes you from an SBML file to Chaste C++ in a few commands.
The generated code is designed to live in a
[Chaste user project](https://chaste.github.io/docs/user-guides/user-projects/).

## Generate code from an SBML file

To generate Chaste C++ code from an SBML file:

```bash
chaste-sbml my_model.xml --output-dir src/
```

This reads `my_model.xml` and writes an
[ODE system](https://chaste.github.io/docs/user-tutorials/solvingodes/) into
`src/`, along with a [placeholder test](generated-code.md#the-placeholder-test).

## Selecting a model type

The `--model-type` option selects what kind of Chaste classes to generate:

`generic`
: Generates a Chaste ODE system only (the default).

`srn`
: Also generates a Chaste subcellular reaction network model that
  wraps the ODE system, for use inside a cell in a cell-based simulation.

`cell-cycle`
: Also generates a Chaste cell-cycle model that wraps the ODE system,
  driving cell division.

```bash
chaste-sbml my_model.xml --model-type srn --output-dir src/
```

:::{seealso}
[Anatomy of generated code](generated-code.md) explains the contents of the
generated files.
:::

## Getting the base classes

The generated code uses (via `#include` statements) a small set of extra C++
base classes that are shipped with the package, among them:

- `AbstractSbmlOdeSystem.hpp`
- `AbstractSbmlCellCycleModel.hpp`
- `AbstractSbmlSrnModel.hpp`
- `SbmlMath.hpp`

See [the base classes](generated-code.md#the-base-classes) for what each one does.

Copy these into your project with:

```bash
chaste-sbml --copy-base-classes --output-dir src/
```

:::{seealso}
- [Command-line options](command-line.md): every flag in detail.
- [Using generated code](using-in-chaste.md): wiring the output into
  a Chaste user project.
:::
