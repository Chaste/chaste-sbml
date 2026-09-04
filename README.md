![build_and_test](https://github.com/Chaste/chaste-sbml/actions/workflows/build_and_test.yml/badge.svg?branch=develop)
![sbml_test_suite](https://github.com/Chaste/chaste-sbml/actions/workflows/sbml_test_suite.yml/badge.svg?branch=develop)

# Import SBML models into Chaste

This project is a Python code generator that reads an [SBML](https://sbml.org) model and emits C++ classes compatible
with [Chaste](https://chaste.github.io), so you can run SBML-defined biology models inside Chaste simulations.

## Installation

### Install with `pipx` (recommended)

```sh
pipx install git+https://github.com/Chaste/chaste-sbml@develop
```

### Install with `pip`

Create and activate a virtual environment (optional)

```sh
python3 -m venv .venv
source .venv/bin/activate
```

Install the package

```sh
python3 -m pip install git+https://github.com/Chaste/chaste-sbml@develop
```

## Development

### Getting the code

Clone the repository

```sh
git clone https://github.com/Chaste/chaste-sbml
cd chaste-sbml
```

Create and activate a virtual environment

```sh
python3 -m venv .venv
source .venv/bin/activate
```

Install in editable mode with development dependencies

```sh
python3 -m pip install -e ."[dev]"
```

Run tests

```sh
python3 -m pytest
```

## Usage

Generating Chaste C++ code from an SBML file is the default action. Passing
`--copy-base-classes` switches to copying the C++ base classes instead.

```
usage: chaste-sbml [-h] [--version] [--copy-base-classes]
                   [--output-dir OUTPUT_DIR]
                   [--model-type [{generic,srn,cell-cycle}]]
                   [--tests | --no-tests] [--timescale {ms,s,m,h}]
                   [--test-output-dir TEST_OUTPUT_DIR]
                   [sbml_file]

positional arguments:
  sbml_file             The SBML file to convert

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --copy-base-classes   Copy the C++ base classes the generated code depends
                        on, instead of generating code
  --output-dir OUTPUT_DIR
                        The directory to place output files in
  --model-type [{generic,srn,cell-cycle}]
                        The type of model to generate
  --tests, --no-tests   Generate placeholder test files (default: on)
  --timescale {ms,s,m,h}
                        The model's native time unit, used to convert
                        derivatives to Chaste's hours (auto-detected if omitted)
  --test-output-dir TEST_OUTPUT_DIR
                        The directory to place generated test files in
                        (defaults to --output-dir)
```

### Generate code

Generate Chaste C++ code from an SBML file:

```sh
chaste-sbml my_model.xml --model-type srn --output-dir src/
```

By default this also emits a placeholder test (`Test<Model>Sbml.hpp`), a
CxxTest skeleton with a suite for the ODE system, and the SRN/cell-cycle
model where applicable. Pass `--no-tests` to skip it, or `--test-output-dir` to
place the placeholder test somewhere other than `--output-dir`.

### Copy the base classes

The generated code `#include`s and subclasses a set of C++ base classes (for example
`AbstractSbmlOdeSystem`). These are shipped with the package; copy them into your project so
they match the installed version of `chaste-sbml`:

```sh
chaste-sbml --copy-base-classes --output-dir src/
```

In this mode only `--output-dir` applies; passing an SBML file or a generation option (such as
`--model-type` or `--no-tests`) is a usage error.
