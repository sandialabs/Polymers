# Polymers Modeling Library

[![website](https://img.shields.io/badge/GitHub-website-6e5494?logo=github)](https://sandialabs.github.io/Polymers)
[![examples](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/binder.svg)](https://mybinder.org/v2/gh/sandialabs/Polymers/main)
[![chat](https://img.shields.io/badge/Discord-chat-%237289da.svg?logo=discord&color=5865F2&logoColor=FFFFFF)](https://discord.gg/9gy8tTktD5)

The library is implemented entirely in Rust, including the Python API. The Julia API calls the Rust library.

## Python

[![docs (stable)](https://img.shields.io/badge/Docs-stable-8CA1AF?logo=readthedocs)](https://polymers.readthedocs.io/en/stable)
[![docs (latest)](https://img.shields.io/badge/Docs-latest-8CA1AF?logo=readthedocs)](https://polymers.readthedocs.io/en/latest)
[![pypi](https://img.shields.io/pypi/v/polymers?logo=pypi&logoColor=FBE072&label=PyPI&color=4B8BBE)](https://pypi.org/project/polymers)

The library can be installed as a Python package using

```shell
pip install polymers
```

If Rust is installed, the latest edition of the library can be installed from the GitHub repository:

```shell
git clone git@github.com:sandialabs/Polymers.git
cd Polymers/
pip install maturin
maturin build --features python
pip install target/wheels/*.whl
```

## Julia

[![docs (stable)](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/julia-docs-stable.svg)](https://sandialabs.github.io/Polymers/julia/docs/stable/polymers)
[![docs (latest)](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/julia-docs-latest.svg)](https://sandialabs.github.io/Polymers/julia/docs/latest/polymers)

The latest edition of the library can be installed as a Julia package using

```julia
using Pkg
Pkg.add(url="https://github.com/sandialabs/Polymers")
Pkg.build("Polymers")
```

## Rust

[![docs (stable)](https://img.shields.io/badge/Docs-stable-f46623?logo=rust&logoColor=000000)](https://docs.rs/crate/polymers)
[![docs (latest)](https://img.shields.io/badge/Docs-latest-f46623?logo=rust&logoColor=000000)](https://sandialabs.github.io/Polymers/rust/docs/latest/polymers)
[![crates](https://img.shields.io/crates/v/polymers?logo=rust&logoColor=000000&label=Crates&color=32592f)](https://crates.io/crates/polymers)

The library can be used in an existing Rust project by adding the `polymers` crate as a dependency in Cargo.toml:

```toml
[dependencies]
polymers = "*"
```
The `*` represents the newest released version of the crate, and should be changed to a specific version.

To use the latest edition of the library, add the GitHub repository to Cargo.toml:

```toml
[dependencies]
regex = { git = "https://github.com/sandialabs/polymers" }
```

## Citation

[![doi](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.7041983-blue)](https://doi.org/10.5281/zenodo.7041983)

Michael R. Buche. Polymers Modeling Library. [Zenodo (2023)](https://doi.org/10.5281/zenodo.7041983).

## Copyright

[![license](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/bsd3c.svg)](https://github.com/sandialabs/polymers/blob/main/LICENSE)

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

