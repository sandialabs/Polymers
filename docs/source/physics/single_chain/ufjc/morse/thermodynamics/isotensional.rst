Morse-FJC model thermodynamics (isotensional)
=============================================

.. toctree::
   :maxdepth: 1

   Asymptotic <isotensional/asymptotic>
   Legendre <isotensional/legendre>

.. autoclass:: polymers.physics.single_chain.ufjc.morse.thermodynamics.isotensional::MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: link_stiffness
   .. autoattribute:: link_energy
   .. autoattribute:: asymptotic
   .. autoattribute:: legendre
   .. automethod:: end_to_end_length(force, temperature)
   .. automethod:: end_to_end_length_per_link(force, temperature)
   .. automethod:: nondimensional_end_to_end_length(nondimensional_force)
   .. automethod:: nondimensional_end_to_end_length_per_link(nondimensional_force)
   .. automethod:: gibbs_free_energy(force, temperature)
   .. automethod:: gibbs_free_energy_per_link(force, temperature)
   .. automethod:: relative_gibbs_free_energy(force, temperature)
   .. automethod:: relative_gibbs_free_energy_per_link(force, temperature)
   .. automethod:: nondimensional_gibbs_free_energy(nondimensional_force, temperature)
   .. automethod:: nondimensional_gibbs_free_energy_per_link(nondimensional_force, temperature)
   .. automethod:: nondimensional_relative_gibbs_free_energy(nondimensional_force)
   .. automethod:: nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)

.. raw::
 html

   <hr>

.. footbibliography::
