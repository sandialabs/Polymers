FJC model thermodynamics (isometric/legendre)
=============================================

.. autoclass:: polymers.physics.single_chain.fjc.thermodynamics.isometric.legendre::FJC(number_of_links, link_length, hinge_mass)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
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
   .. automethod:: gibbs_free_energy(end_to_end_length, temperature)
   .. automethod:: gibbs_free_energy_per_link(end_to_end_length, temperature)
   .. automethod:: relative_gibbs_free_energy(end_to_end_length, temperature)
   .. automethod:: relative_gibbs_free_energy_per_link(end_to_end_length, temperature)
   .. automethod:: nondimensional_gibbs_free_energy(nondimensional_end_to_end_length_per_link, temperature)
   .. automethod:: nondimensional_gibbs_free_energy_per_link(nondimensional_end_to_end_length_per_link, temperature)
   .. automethod:: nondimensional_relative_gibbs_free_energy(nondimensional_end_to_end_length_per_link)
   .. automethod:: nondimensional_relative_gibbs_free_energy_per_link(nondimensional_end_to_end_length_per_link)
