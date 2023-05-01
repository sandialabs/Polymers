WLC model thermodynamics (isometric)
====================================

.. toctree::
   :maxdepth: 1

   Legendre <isometric/legendre>

.. autoclass:: polymers.physics.single_chain.wlc.thermodynamics.isometric::WLC(number_of_links, link_length, hinge_mass, persistance_length)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: persistance_length
   .. autoattribute:: legendre
   .. automethod:: force(end_to_end_length, temperature)
   .. automethod:: nondimensional_force(nondimensional_end_to_end_length_per_link)
   .. automethod:: helmholtz_free_energy(end_to_end_length, temperature)
   .. automethod:: helmholtz_free_energy_per_link(end_to_end_length, temperature)
   .. automethod:: relative_helmholtz_free_energy(end_to_end_length, temperature)
   .. automethod:: relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy(nondimensional_end_to_end_length_per_link, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link, temperature)
   .. automethod:: nondimensional_relative_helmholtz_free_energy(nondimensional_end_to_end_length_per_link)
   .. automethod:: nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link)
   .. automethod:: equilibrium_distribution(end_to_end_length)
   .. automethod:: nondimensional_equilibrium_distribution(nondimensional_end_to_end_length_per_link)
   .. automethod:: equilibrium_radial_distribution(end_to_end_length)
   .. automethod:: nondimensional_equilibrium_radial_distribution(nondimensional_end_to_end_length_per_link)

.. raw::
 html

   <hr>

.. footbibliography::
