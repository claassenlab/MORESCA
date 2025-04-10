# Installation

We strongly recommend to install MORESCA into a virtual environment. Here, we use Conda:

    conda create -n <env_name> python=3.12
    conda activate <env_name>

Then, simply install MORESCA with pip:

    pip install moresca

```{important}
If you want to use Python 3.13 on MacOS, make sure to use GCC>=16. This is required for compiling scikit-misc. See [this discussion](https://stackoverflow.com/questions/48174684/fortran-codes-wont-compile-on-mac-with-gfortran) for advice.
```

## Contributing

For contribution purposes, you should clone MORESCA from GitHub and install it in dev mode:

    git clone git@github.com:claassenlab/MORESCA.git
    cd MORESCA
    pip install -e ".[dev]"
    pre-commit install

This additionally installs `ruff` and `pytest`, which we use for formatting and code style control. Please run these before you commit new code.
Additionally, it will set up a pre-commit hook to run `ruff`.