# Command-line options

The `chaste-sbml` command has two modes:

- **Generate** (the default): convert an SBML file to Chaste C++ code.
- **Copy base classes** (`--copy-base-classes`): copy out the C++ base classes
  the generated code depends on, instead of generating.

```text
usage: chaste-sbml [-h] [--version] [--copy-base-classes]
                   [--output-dir OUTPUT_DIR]
                   [--model-type [{generic,srn,cell-cycle}]]
                   [--tests | --no-tests] [--timescale {ms,s,m,h}]
                   [--test-output-dir TEST_OUTPUT_DIR]
                   [sbml_file]
```

## Positional argument

`sbml_file`
: The SBML file to convert. Required when generating; must **not** be given in
  copy mode.

## General options

`-h`, `--help`
: Show the help message and exit.

`--version`
: Show the program's version number and exit.

`--output-dir OUTPUT_DIR`
: The directory to place output files in. Applies to both modes. Defaults to the
  current directory.

## Generation options

These apply only when generating code. Passing any of them together with
`--copy-base-classes` is a usage error.

`--model-type [{generic,srn,cell-cycle}]`
: The type of Chaste class to generate (default: `generic`).

  `generic`
  : An ODE system only.

  `srn`
  : A subcellular reaction network model wrapping the ODE system.

  `cell-cycle`
  : A cell-cycle model wrapping the ODE system.

  See [Anatomy of generated code](generated-code.md) for the files each type
  produces.

`--tests`, `--no-tests`
: Generate placeholder test files (default: on). Use `--no-tests` to skip them.
  The placeholder is a CxxTest skeleton — see
  [the placeholder test](generated-code.md#the-placeholder-test).

`--timescale {ms,s,m,h}`
: The model's native time unit — milliseconds, seconds, minutes, or hours — used
  to scale derivatives into Chaste's hours. Overrides auto-detection from the
  SBML; omit it to auto-detect. See [Time units](#time-units) below.

`--test-output-dir TEST_OUTPUT_DIR`
: The directory to place generated test files in. Defaults to `--output-dir`, so
  set this to keep tests separate from the model source.

## Copying the base classes

The generated code subclasses base classes such as `AbstractSbmlOdeSystem`,
`AbstractSbmlSrnModel`, and `AbstractSbmlCellCycleModel`. They ship with the
package; copy them into your project so they match the installed version of the
tool:

```bash
chaste-sbml --copy-base-classes --output-dir src/
```

:::{important}
In copy mode only `--output-dir` applies. Passing an SBML file or any generation
option (such as `--model-type` or `--no-tests`) is a usage error, so mistakes
are reported rather than silently ignored.
:::

## Time units

Chaste works in hours, so the generator scales the model's derivatives from the
model's native time unit into hours. The unit is resolved with this precedence:

1. An explicit `--timescale` override wins. If it contradicts a unit declared by
   the model, a warning is logged and the override is used.
2. Otherwise a unit declared by the model is used.
3. Otherwise the SBML default applies:
   - **Level 2** predefines the time unit as **seconds**.
   - **Level 3** leaves an unset time unit undefined, so **no conversion** is
     applied.

When no unit is declared and none is given, a warning hints at `--timescale`.
Auto-detection only recognises a time unit that reduces to seconds (`second`,
`millisecond`, `minute`, `hour`); a composite or otherwise non-second time unit
is treated as undeterminable, and you should pass `--timescale` explicitly.

:::{seealso}
[First steps](first-steps.md) for a quick tour, and the
[tutorials](tutorials/goldbeter1991.md) for full worked examples.
:::
