Polymers Modeling Library
=========================

|stable| |latest|  

This is the documentation for Python API, which is implemented in Rust.

Installation
------------

|pypi|

The library can be installed as a Python package:

.. code-block:: sh

   pip install polymers

If Rust is installed, the latest edition of the library can be installed from the GitHub repository:

.. code-block:: sh

   git clone git@github.com:sandialabs/Polymers.git
   cd Polymers/
   pip install maturin
   maturin build --features python
   pip install target/wheels/*.whl

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

.. |zenodo| image:: https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.7041983-blue
   :target: https://doi.org/10.5281/zenodo.7041983

.. |license| image:: https://raw.githubusercontent.com/sandialabs/Polymers/main/pages/assets/images/bsd3c.svg
   :target: https://github.com/sandialabs/polymers/blob/main/LICENSE

.. toctree::
   :hidden:
   :caption: Modules

   Physics <physics>

.. toctree::
   :hidden:
   :caption: Examples

   Physics <physics/examples>

.. toctree::
   :hidden:
   :caption: Indices and Tables

   General Index <./genindex.html#http://>
   Module Index <./py-modindex.html#http://>
