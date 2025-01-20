# SBML → Chaste

`chaste_codegen_sbml` converts [SBML](https://synonym.caltech.edu/) models to C++ Chaste code. It relies on [libSBML](https://synonym.caltech.edu/software/libsbml/) for parsing SBML.

## Installation
Installation via [`pipx`](https://pipx.pypa.io/) is recommended
```
pipx install git+https://github.com/Chaste/chaste-codegen-sbml@develop
```

Alternatively, install via `pip`
```
# Create and activate a virtual environment (optional)
python -m venv sbml-venv
source sbml-venv/bin/activate

# Install
pip install git+https://github.com/Chaste/chaste-codegen-sbml@develop
```

## Usage
```
chaste_codegen_sbml [-h]
```

## Development
### Getting the code
```
# Clone the repository
git clone https://github.com/Chaste/chaste-codegen-sbml
cd chaste-codegen-sbml

# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate

# Install in editable mode along with dev dependencies
pip install -e ."[dev]"
```

### Running tests
```
python -m pytest
```
