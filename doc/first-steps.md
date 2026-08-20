# First steps

This page takes you from an SBML file to compiled-ready C++ in a few commands.
For a line-by-line walkthrough of the generated files, see the
[Goldbeter 1991 tutorial](tutorials/goldbeter1991.md).

## Generate code from an SBML file

Generating Chaste C++ code from an SBML file is the default action:

```bash
chaste-sbml my_model.xml --output-dir src/
```

This reads `my_model.xml` and writes an
[ODE system](generated-code.md#the-ode-system) into `src/`, along with a
[placeholder test](generated-code.md#the-placeholder-test).

## Choose a model type

The `--model-type` option selects what kind of Chaste class to generate:

`generic`
: An ODE system only (the default). Use this for a plain system of ODEs.

`srn`
: A [subcellular reaction network](https://chaste.github.io) (SRN) model that
  wraps the ODE system, for use inside a cell in a cell-based simulation.

`cell-cycle`
: A cell-cycle model that wraps the ODE system, driving cell division.

```bash
chaste-sbml my_model.xml --model-type srn --output-dir src/
```

For the `srn` and `cell-cycle` types, the tool emits the wrapper class
*and* the ODE system it wraps. See
[Anatomy of generated code](generated-code.md) for what each file contains.

## Copy the base classes

The generated code `#include`s and subclasses a small set of C++ base classes
(for example `AbstractSbmlOdeSystem`). These ship with the package; copy them
into your project so they match the installed version of the tool:

```bash
chaste-sbml --copy-base-classes --output-dir src/
```

## Try it on a reference model

The repository ships several reference models under
`chaste_sbml/SbmlRefModels/src/reference/`. If you have a development checkout,
you can generate one straight away:

```bash
chaste-sbml \
  chaste_sbml/SbmlRefModels/src/reference/Goldbeter1991/Goldbeter1991.xml \
  --model-type srn \
  --output-dir /tmp/goldbeter
```

:::{seealso}
- [Command-line options](command-line.md) — every flag in detail.
- [Supported SBML features](sbml-features.md) — what the generator
  handles.
- [Using generated code in Chaste](using-in-chaste.md) — wiring the output into
  a Chaste user project.
:::
