# texespy-dist

Python package for processing spectroscopic data measured with the TEXES instrument at the NASA Infrared Telescope Facility

## Installation

### Installing from GitHub (developer mode)

The latest version of code has to be downloaded from [Github](https://github.com/juanaldayparejo/texespy-dist.git). To do so, type in the command window:

```bash
git clone https://github.com/juanaldayparejo/texespy-dist.git
```

Before installing the library, we recommend users to create and load a new Python [virtual environment](https://docs.python.org/3/library/venv.html) for a clean install:

```bash
python -m venv name_of_virtual_environment/
source name_of_virtual_environment/bin/activate
```

Then move into the package directory:

```bash
cd texespy-dist
```

Finally, we need to install the library. We recommend installing the package but keeping it editable by typing:

```bash
pip install --editable .
```

After installing the library, we need to change the path where the spice kernels are located. To do so:

- Go to the file texespy/data/spice/kernels_esa/mk/esa_generic_v01.tm
- Change the variable PATH_VALUES with the path where your texespy-dist distribution is located.
