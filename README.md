# Polymers Modeling Library

[![GitHub Pages](https://img.shields.io/badge/GitHub-pages-6e5494?logo=github)](https://sandialabs.github.io/Polymers)
[![Discord](https://img.shields.io/badge/Discord-chat-%237289da.svg?logo=discord&color=5865F2&logoColor=FFFFFF)](https://discord.gg/yC6dbPuc)
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

Michael R. Buche. Polymers Modeling Library. [Zenodo (2022)](https://doi.org/10.5281/zenodo.7041983).

## Copyright

[![License](https://img.shields.io/github/license/sandialabs/polymers?label=License&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAPZJREFUOE+tk2ERwjAUgxMHOAAUYIEpAAngAAccDuYAJIACkIADhgMchMuu5Uq7Aj/or92a9700fSU+LEktgBHJVU3GL4CL90nOfwZIsngLIC8ybEeyh8bVO5B0BTBL/t8BHABE6/F7nGg6ktMIsHAfNtckXWDw2xEkFboIsN1zADTR5gCg0EWAO00C4EayqTgodJTkwpsDCgAHOCXZpQ6qOkk+7zJx0AE4+u4zwLBOkgA8Quo24aA8PHb3CrGqC/e+AbAIRzgBaB1kBnCApe7vk5gC82scala8hYGpTOv66StGOevqEONM5E0N6Kf07S18yuHb3hPwkpAEoqucdwAAAABJRU5ErkJggg==)](https://github.com/sandialabs/polymers/blob/main/LICENSE)

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
