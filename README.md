# Chaste code generation from SBML

The `chaste_codegen_sbml` module takes [SBML](https://synonym.caltech.edu/) models as input and uses [libSBML](https://synonym.caltech.edu/software/libsbml/) to parse the models and output C++ code for Chaste.

## Installation
### User install
Installation via [`pipx`](https://pipx.pypa.io/) is recommended
```
pipx install git+https://github.com/Chaste/chaste-codegen-sbml@develop
```

Alternatively, install via `pip` (optionally in a virtual environment)
```
# Create and activate a virtual environment
python -m venv sbml-venv

# Install
pip install git+https://github.com/Chaste/chaste-codegen-sbml@develop
```

### Developer install
Clone the repository and install the package in editable mode along with the `dev` dependencies:
```
git clone https://github.com/Chaste/chaste-codegen-sbml
cd chaste-codegen-sbml
pip install -e .[dev]
```

## Usage
```
chaste_codegen_sbml [-h]
```

## Testing
