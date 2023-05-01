WLC model thermodynamics (isometric/legendre)
=============================================

.. autoclass:: polymers.physics.single_chain.wlc.thermodynamics.isometric.legendre::WLC(number_of_links, link_length, hinge_mass, persistance_length)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: persistance_length
   .. automethod:: gibbs_free_energy(end_to_end_length, temperature)
   .. automethod:: gibbs_free_energy_per_link(end_to_end_length, temperature)
   .. automethod:: relative_gibbs_free_energy(end_to_end_length, temperature)
   .. automethod:: relative_gibbs_free_energy_per_link(end_to_end_length, temperature)
   .. automethod:: nondimensional_gibbs_free_energy(nondimensional_end_to_end_length_per_link, temperature)
   .. automethod:: nondimensional_gibbs_free_energy_per_link(nondimensional_end_to_end_length_per_link, temperature)
   .. automethod:: nondimensional_relative_gibbs_free_energy(nondimensional_end_to_end_length_per_link)
   .. automethod:: nondimensional_relative_gibbs_free_energy_per_link(nondimensional_end_to_end_length_per_link)
