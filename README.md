# Pazuzu
Flow simulation based on the Discontinuous Galerkin method. The DG method is a compact scheme, therefore it tends to be more accurate for its stencil size. Due to the smaller stencil it is possible to use GPUs to furter increase the convergence time.

![build](https://github.com/dzilles/pazuzu/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/gh/dzilles/pazuzu/branch/project_setup/graph/badge.svg?token=)](https://codecov.io/gh/dzilles/pazuzu)

# Installation

0. For the first time create a virtual environment using `python -m venv ./.env`
1. Activate the virtual environment `source .env/bin/activate`
2. Compile code with `maturin develop`
3. Execute testcases with `cargo test` and `pythontest`


