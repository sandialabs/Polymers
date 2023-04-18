Morse-FJC model thermodynamics (isometric/legendre)
===================================================

.. autoclass:: polymers.physics.single_chain.ufjc.morse.thermodynamics.isometric.legendre::MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: link_stiffness
   .. autoattribute:: link_energy
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
