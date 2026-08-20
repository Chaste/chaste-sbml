# Installation

## User install

### Using pipx (recommended)

Use [pipx](https://pipx.pypa.io) to install `chaste-sbml` into its own isolated
environment and put the command on your `PATH`, without touching other Python
packages:

```bash
pipx install git+https://github.com/Chaste/chaste-sbml@0.0.1
```

:::{note}
Replace `@0.0.1` with `@x.y.z` to install a specific released version, or omit
it to install the latest development version.
:::

Verify the installation:

```bash
chaste-sbml --version
```

### Using pip

Create and activate a virtual environment (optional):

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the package:

```bash
python3 -m pip install git+https://github.com/Chaste/chaste-sbml@0.0.1
```

:::{note}
Replace `@0.0.1` with `@x.y.z` to install a specific released version, or omit
it to install the latest development version.
:::

Verify the installation:

```bash
chaste-sbml --version
```

## Developer install

Clone the repository:

```bash
git clone https://github.com/Chaste/chaste-sbml
cd chaste-sbml
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the local copy in editable mode with the `dev` dependencies:

```bash
python3 -m pip install -e ."[dev]"
```

Run the test suite:

```bash
python3 -m pytest
```

To build the documentation locally, install the `docs` extra and run `make`:

```bash
python3 -m pip install -e ."[docs]"
make -C doc livehtml      # live-reloading preview
```
