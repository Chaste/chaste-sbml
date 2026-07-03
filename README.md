![build_and_test](https://github.com/Chaste/chaste-codegen-sbml/actions/workflows/build_and_test.yml/badge.svg?branch=develop)
![sbml_test_suite](https://github.com/Chaste/chaste-codegen-sbml/actions/workflows/sbml_test_suite.yml/badge.svg?branch=develop)

# Import SBML models into Chaste

This project is a Python code generator that reads an [SBML](https://sbml.org) model and emits C++ classes compatible
with [Chaste](https://chaste.github.io), so you can run SBML-defined biology models inside Chaste simulations.

## Installation

### Install with `pipx` (recommended)

```sh
pipx install git+https://github.com/Chaste/chaste-codegen-sbml@develop
```

### Install with `pip`

Create and activate a virtual environment (optional)

```sh
python3 -m venv .venv
source .venv/bin/activate
```

Install the package

```sh
python3 -m pip install git+https://github.com/Chaste/chaste-codegen-sbml@develop
```

## Development

### Getting the code

Clone the repository

```sh
git clone https://github.com/Chaste/chaste-codegen-sbml
cd chaste-codegen-sbml
```

Create and activate a virtual environment

```sh
python3 -m .venv .venv
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

```
usage: chaste_codegen_sbml [-h] [--version] [--output-dir OUTPUT_DIR] [--model-type [{generic,srn,cell-cycle}]] sbml_file

Convert SBML models to Chaste C++ code

positional arguments:
  sbml_file             The SBML file to convert

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --output-dir OUTPUT_DIR
                        The directory to place output files in
  --model-type [{generic,srn,cell-cycle}]
                        The type of model to generate
```
