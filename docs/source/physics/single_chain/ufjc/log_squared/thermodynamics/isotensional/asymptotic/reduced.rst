Log-squared-FJC model thermodynamics (isotensional/asymptotic/reduced)
======================================================================

.. toctree::
   :maxdepth: 1

   Legendre <reduced/legendre>

.. autoclass:: polymers.physics.single_chain.ufjc.log_squared.thermodynamics.isotensional.asymptotic.reduced::LOGSQUAREDFJC(number_of_links, link_length, hinge_mass, link_stiffness)

   .. autoattribute:: number_of_links
   .. autoattribute:: link_length
   .. autoattribute:: hinge_mass
   .. autoattribute:: link_stiffness
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
