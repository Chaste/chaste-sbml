![build_and_test]
(https://github.com/Chaste/chaste-codegen-sbml/actions/workflows/build_and_test.yml/badge.svg?branch=develop)

# Import SBML models into Chaste

The `chaste_codegen_sbml` tool converts [SBML](https://sbml.org) models into C++
Chaste code using [libSBML](https://github.com/sbmlteam/libsbml).

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
