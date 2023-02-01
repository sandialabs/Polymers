# Polymers Modeling Library

[![docs (stable)](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/julia-docs-stable.svg)](https://sandialabs.github.io/Polymers/julia/docs/stable/polymers)
[![docs (latest)](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/julia-docs-latest.svg)](https://sandialabs.github.io/Polymers/julia/docs/latest/polymers)

This is the documentation for Julia API, which calls the Rust library.

## Installation

The latest edition of the library can be installed as a Julia package:

```julia
using Pkg
Pkg.add(url="https://github.com/sandialabs/Polymers")
Pkg.build("Polymers")
```

## Citation

[![doi](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.7041983-blue)](https://doi.org/10.5281/zenodo.7041983)

Michael R. Buche. Polymers Modeling Library. [Zenodo (2023)](https://doi.org/10.5281/zenodo.7041983).

## Copyright

[![license](https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/bsd3c.svg)](https://github.com/sandialabs/polymers/blob/main/LICENSE)

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

## Documentation


```@autodocs
Modules = [Polymers.Physics.SingleChain.Fjc.Thermodynamics.Isometric]
```
