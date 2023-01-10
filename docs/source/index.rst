Polymers Modeling Library
=========================

|stable| |latest|  

This is the documentation for Python API, which is implemented in Rust.

Installation
------------

|pypi| |conda| |docker| |ghcr|  

The library can be installed as a Python package using

.. code-block:: sh

   pip install polymers

or as Python package within a Conda environment using

.. code-block:: sh

   conda install --channel mrbuche polymers

The latest edition of the library can be installed from the main branch of the GitHub repository:

.. code-block:: sh

   git clone git@github.com:sandialabs/Polymers.git
   cd Polymers/
   pip install maturin
   maturin build --features python
   pip install target/wheels/*.whl

Docker images are available for all stable versions that have been released, as well as the latest edition. These images have everything necessary installed, and can be pulled from Docker Hub:

.. code-block:: sh

   docker pull mrbuche/polymers

as well as from the GitHub container registry:

.. code-block:: sh

   docker pull ghcr.io/sandialabs/polymers

Citation
--------

|zenodo|  

Michael R. Buche. Polymers Modeling Library. `Zenodo (2023) <https://doi.org/10.5281/zenodo.7041983>`_.

Copyright
---------

|license|  

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

.. |stable| image:: https://img.shields.io/badge/Docs-stable-8CA1AF?logo=readthedocs
   :target: https://polymers.readthedocs.io/en/stable

.. |latest| image:: https://img.shields.io/badge/Docs-latest-8CA1AF?logo=readthedocs
   :target: https://polymers.readthedocs.io/en/latest

.. |pypi| image:: https://img.shields.io/pypi/v/polymers?logo=pypi&logoColor=FBE072&label=PyPI&color=4B8BBE
   :target: https://pypi.org/project/polymers

.. |conda| image:: https://img.shields.io/conda/v/mrbuche/polymers.svg?logo=anaconda&color=3EB049&label=Anaconda
   :target: https://anaconda.org/mrbuche/polymers

.. |docker| image:: https://img.shields.io/docker/v/mrbuche/polymers?color=0db7ed&label=Docker%20Hub&logo=docker&logoColor=0db7ed
   :target: https://hub.docker.com/r/mrbuche/polymers

.. |ghcr| image:: https://img.shields.io/badge/GitHub-latest-6e5494?logo=github
   :target: https://github.com/sandialabs/Polymers/pkgs/container/polymers

.. |zenodo| image:: https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.7041983-blue
   :target: https://doi.org/10.5281/zenodo.7041983

.. |license| image:: https://img.shields.io/github/license/sandialabs/polymers?label=License&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAPZJREFUOE+tk2ERwjAUgxMHOAAUYIEpAAngAAccDuYAJIACkIADhgMchMuu5Uq7Aj/or92a9700fSU+LEktgBHJVU3GL4CL90nOfwZIsngLIC8ybEeyh8bVO5B0BTBL/t8BHABE6/F7nGg6ktMIsHAfNtckXWDw2xEkFboIsN1zADTR5gCg0EWAO00C4EayqTgodJTkwpsDCgAHOCXZpQ6qOkk+7zJx0AE4+u4zwLBOkgA8Quo24aA8PHb3CrGqC/e+AbAIRzgBaB1kBnCApe7vk5gC82scala8hYGpTOv66StGOevqEONM5E0N6Kf07S18yuHb3hPwkpAEoqucdwAAAABJRU5ErkJggg==
   :target: https://github.com/sandialabs/polymers/blob/main/LICENSE

.. toctree::
   :hidden:
   :caption: Modules

   Physics <polymers.physics>

.. toctree::
   :hidden:
   :caption: Examples

   example_notebook

.. toctree::
   :hidden:
   :caption: Indices and Tables

   General Index <./genindex.html#http://>
   Module Index <./py-modindex.html#http://>