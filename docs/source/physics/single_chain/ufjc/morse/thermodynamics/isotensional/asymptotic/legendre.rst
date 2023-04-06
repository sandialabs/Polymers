Morse-FJC model thermodynamics (isotensional/asymptotic/legendre)
=================================================================

.. autoclass:: polymers.physics.single_chain.ufjc.morse.thermodynamics.isotensional.asymptotic.legendre::MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: link_stiffness
   .. autoattribute:: link_energy
   .. automethod:: helmholtz_free_energy(force, temperature)
   .. automethod:: helmholtz_free_energy_per_link(force, temperature)
   .. automethod:: relative_helmholtz_free_energy(force, temperature)
   .. automethod:: relative_helmholtz_free_energy_per_link(force, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy(nondimensional_force, temperature)
   .. automethod:: nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)
   .. automethod:: nondimensional_relative_helmholtz_free_energy(nondimensional_force)
   .. automethod:: nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)
