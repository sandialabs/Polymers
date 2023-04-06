Lennard-Jones-FJC model thermodynamics (isotensional/asymptotic/reduced/legendre)
=================================================================================

.. autoclass:: polymers.physics.single_chain.ufjc.lennard_jones.thermodynamics.isotensional.asymptotic.reduced.legendre::LENNARDJONESFJC(number_of_links, link_length, hinge_mass, link_stiffness)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: link_stiffness
   .. automethod:: helmholtz_free_energy(force, temperature)
   .. automethod:: helmholtz_free_energy_per_link(force, temperature)
   .. automethod:: relative_helmholtz_free_energy(force, temperature)
   .. automethod:: relative_helmholtz_free_energy_per_link(force, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy(nondimensional_force, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)
   .. automethod:: nondimensional_relative_helmholtz_free_energy(nondimensional_force)
   .. automethod:: nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)
