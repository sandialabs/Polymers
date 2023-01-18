FJC model thermodynamics (modified canonical/asymptotic/weak potential)
=======================================================================

.. autoclass:: polymers.physics.single_chain.fjc.thermodynamics.modified_canonical.asymptotic.weak_potential::FJC(number_of_links, link_length, hinge_mass)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. automethod:: end_to_end_length(potential_distance, potential_stiffness, temperature)
   .. automethod:: end_to_end_length_per_link(potential_distance, potential_stiffness, temperature)
   .. automethod:: nondimensional_end_to_end_length(nondimensional_potential_distance, nondimensional_potential_stiffness)
   .. automethod:: nondimensional_end_to_end_length_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness)
   .. automethod:: force(potential_distance, potential_stiffness)
   .. automethod:: nondimensional_force(nondimensional_potential_distance, nondimensional_potential_stiffness)
   .. automethod:: gibbs_free_energy(potential_distance, potential_stiffness, temperature)
   .. automethod:: gibbs_free_energy_per_link(potential_distance, potential_stiffness, temperature)
   .. automethod:: relative_gibbs_free_energy(potential_distance, potential_stiffness, temperature)
   .. automethod:: relative_gibbs_free_energy_per_link(potential_distance, potential_stiffness, temperature)
   .. automethod:: nondimensional_gibbs_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
   .. automethod:: nondimensional_gibbs_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
   .. automethod:: nondimensional_relative_gibbs_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness)
   .. automethod:: nondimensional_relative_gibbs_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness)
