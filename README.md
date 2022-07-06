# Globular_Clusters
This directory contains code pertaining to the paper 'Advancing globular cluster constraints on the axion-photon coupling'.

'Base' directories contain `run_star_extras` files and inlists for each considered scheme for modelling convective boundaries. These were setup to run on a cluster and have two inlist files: `inlist_lowinter`, which contains standard information between all simulations for that scheme, and `inlist_cluster`, which contains input physics specific to each run.

There are six base directories included:
 
  *`base_SO` - input physics for HB evolution using standard overshoot

  *`base_SC` - input physics for HB evolution using semiconvection

  *`base_PM` - input physics for HB evolution using predictive mixing

  *`base_CP` - input physics for HB evolution using convective premixing

  *`base_AGB` - input physics from the TAHB through the AGB

  *`base_RGB` - input physics from the pre-MS to the ZAHB
  
The `run_star_extras` files include axion energy-loss with the effects of electron degeneracy via the methods of G.G. Raffelt and D.S.P. Dearborn, Bounds on Hadronic Axions From Stellar Evolution, Phys. Rev. D 36 (1987) 2211, and based on the parametrisation in A. Ayala, I. Domínguez, M. Giannotti, A. Mirizzi and O. Straniero, Revisiting the bound on axion-photon coupling from Globular Clusters, Phys. Rev. Lett. 113 (2014) 191302 [1406.6053]. These are included via additional terms to the neutrino energy-loss module. 

Several features in the `run_star_extras` file and inlists originate from the MIST input physics (J. Choi, A. Dotter, C. Conroy, M. Cantiello, B. Paxton and B.D. Johnson, Mesa Isochrones and Stellar Tracks (MIST). I. Solar-scaled Models, Astrophys. J. 823 (2016) 102 [1604.08592]).

These simulations were performed using `MESA v12778`.
