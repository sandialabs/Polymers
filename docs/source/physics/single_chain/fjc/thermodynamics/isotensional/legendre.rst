FJC model thermodynamics (isotensional/legendre)
================================================

.. autoclass:: polymers.physics.single_chain.fjc.thermodynamics.isotensional.legendre::FJC(number_of_links, link_length, hinge_mass)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. automethod:: helmholtz_free_energy(force, temperature)
   .. automethod:: helmholtz_free_energy_per_link(force, temperature)
   .. automethod:: relative_helmholtz_free_energy(force, temperature)
   .. automethod:: relative_helmholtz_free_energy_per_link(force, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy(nondimensional_force, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)
   .. automethod:: nondimensional_relative_helmholtz_free_energy(nondimensional_force)
   .. automethod:: nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)
