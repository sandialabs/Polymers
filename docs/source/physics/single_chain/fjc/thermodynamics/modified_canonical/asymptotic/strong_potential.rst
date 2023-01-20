FJC model thermodynamics (modified canonical/asymptotic/strong potential)
=========================================================================

.. autoclass:: polymers.physics.single_chain.fjc.thermodynamics.modified_canonical.asymptotic.strong_potential::FJC(number_of_links, link_length, hinge_mass)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. automethod:: force(potential_distance, potential_stiffness, temperature)
   .. automethod:: nondimensional_force(nondimensional_potential_distance, nondimensional_potential_stiffness)
   .. automethod:: helmholtz_free_energy(potential_distance, potential_stiffness, temperature)
   .. automethod:: helmholtz_free_energy_per_link(potential_distance, potential_stiffness, temperature)
   .. automethod:: relative_helmholtz_free_energy(potential_distance, potential_stiffness, temperature)
   .. automethod:: relative_helmholtz_free_energy_per_link(potential_distance, potential_stiffness, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
   .. automethod:: nondimensional_relative_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness)
   .. automethod:: nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness)
