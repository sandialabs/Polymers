FJC model thermodynamics (isometric)
====================================

.. toctree::
   :hidden:
   :maxdepth: 2

   Legendre <polymers.physics.single_chain.fjc.thermodynamics.isometric.legendre>

.. autoclass:: polymers.physics.single_chain.fjc.thermodynamics.isometric::FJC

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: legendre
   .. automethod:: force(end_to_end_length, temperature)
   .. automethod:: nondimensional_force(nondimensional_end_to_end_length_per_link)
   .. automethod:: equilibrium_radial_distribution(end_to_end_length)

.. raw::
 html

   <hr>

**References**

.. bibliography::
   :filter: docname in docnames