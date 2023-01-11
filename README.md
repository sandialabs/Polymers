# Polymers Modeling Library

[![GitHub Pages](https://img.shields.io/badge/GitHub-website-6e5494?logo=github)](https://sandialabs.github.io/Polymers)
[![Binder](https://img.shields.io/badge/Binder-examples-579ACA.svg?logo=data:image/png;base64,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)](https://mybinder.org/v2/gh/sandialabs/Polymers/main)
[![Discord](https://img.shields.io/badge/Discord-chat-%237289da.svg?logo=discord&color=5865F2&logoColor=FFFFFF)](https://discord.gg/9gy8tTktD5)
[![Codefactor](https://img.shields.io/codefactor/grade/github/sandialabs/polymers?label=Codefactor&logo=codefactor&color=00b16a)](https://www.codefactor.io/repository/github/sandialabs/polymers)
[![Codecov](https://img.shields.io/codecov/c/github/sandialabs/polymers?label=Codecov&logo=codecov&flag=rust)](https://codecov.io/gh/sandialabs/polymers)

The library is implemented entirely in Rust, including the Python API.

## Python

[![Read the Docs](https://img.shields.io/badge/Docs-stable-8CA1AF?logo=readthedocs)](https://polymers.readthedocs.io/en/stable)
[![Read the Docs latest](https://img.shields.io/badge/Docs-latest-8CA1AF?logo=readthedocs)](https://polymers.readthedocs.io/en/latest)
[![PyPI](https://img.shields.io/pypi/v/polymers?logo=pypi&logoColor=FBE072&label=PyPI&color=4B8BBE)](https://pypi.org/project/polymers)
[![Anaconda](https://img.shields.io/conda/v/mrbuche/polymers.svg?logo=anaconda&color=3EB049&label=Anaconda)](https://anaconda.org/mrbuche/polymers/)
[![Docker Hub](https://img.shields.io/docker/v/mrbuche/polymers?color=0db7ed&label=Docker%20Hub&logo=docker&logoColor=0db7ed)](https://hub.docker.com/r/mrbuche/polymers)
[![GitHub Container Registry](https://img.shields.io/badge/GitHub-latest-6e5494?logo=github)](https://github.com/sandialabs/Polymers/pkgs/container/polymers)

The library can be installed as a Python package using

```shell
pip install polymers
```

or as Python package within a Conda environment using

```shell
conda install --channel mrbuche polymers
```

The latest edition of the library can be installed from the main branch of the GitHub repository:

```shell
git clone git@github.com:sandialabs/Polymers.git
cd Polymers/
pip install maturin
maturin build --features python
pip install target/wheels/*.whl
```

Docker images are available for all stable versions that have been released, as well as the latest edition. These images have everything necessary installed, and can be pulled from Docker Hub:

```shell
docker pull mrbuche/polymers
```

as well as from the GitHub container registry:

```shell
docker pull ghcr.io/sandialabs/polymers
```

## Rust

[![Docs.rs](https://img.shields.io/badge/Docs-stable-32592f?logo=rust&logoColor=000000)](https://docs.rs/crate/polymers)
[![Docs.rs latest](https://img.shields.io/badge/Docs-latest-32592f?logo=rust&logoColor=000000)](https://sandialabs.github.io/Polymers/rust/docs/latest/polymers)
[![Crates.io](https://img.shields.io/crates/v/polymers?logo=rust&logoColor=000000&label=Crate&color=32592f)](https://crates.io/crates/polymers)

The library can be used in an existing Rust project by adding the `polymers` crate as a dependency in Cargo.toml:

```toml
[dependencies]
polymers = "*"
```
The asterisk `*` represents the newest released version of the crate, and should be changed to a specific version.

To use the latest edition of the library, add the main branch of the GitHub repository to Cargo.toml:

```toml
[dependencies]
regex = { git = "https://github.com/sandialabs/polymers" }
```

## Citation

[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.7041983-blue)](https://doi.org/10.5281/zenodo.7041983)

Michael R. Buche. Polymers Modeling Library. [Zenodo (2023)](https://doi.org/10.5281/zenodo.7041983).

## Copyright

[![License](https://img.shields.io/github/license/sandialabs/polymers?label=License&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAPZJREFUOE+tk2ERwjAUgxMHOAAUYIEpAAngAAccDuYAJIACkIADhgMchMuu5Uq7Aj/or92a9700fSU+LEktgBHJVU3GL4CL90nOfwZIsngLIC8ybEeyh8bVO5B0BTBL/t8BHABE6/F7nGg6ktMIsHAfNtckXWDw2xEkFboIsN1zADTR5gCg0EWAO00C4EayqTgodJTkwpsDCgAHOCXZpQ6qOkk+7zJx0AE4+u4zwLBOkgA8Quo24aA8PHb3CrGqC/e+AbAIRzgBaB1kBnCApe7vk5gC82scala8hYGpTOv66StGOevqEONM5E0N6Kf07S18yuHb3hPwkpAEoqucdwAAAABJRU5ErkJggg==)](https://github.com/sandialabs/polymers/blob/main/LICENSE)

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
