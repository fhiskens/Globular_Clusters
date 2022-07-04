# Globular_Clusters
Code pertaining to paper 'Advancing globular cluster constraints on the axion-photon coupling'

'Base' directories contain run_star_extras files and inlists for each considered scheme for modelling convective boundaries. These were setup to run on a cluster and have two inlist files: 'inlist_lowinter', which contains standard information between all simulations for that scheme, and 'inlist_cluster', which contains input physics specific to each run.

There are five base directories included.
  base_SO - the base directory for the standard overshoot scheme
  base_SC - the base directory for the semiconvection scheme
  base_PM - the base directory for the predictive mixing scheme
  base_CP - the base directory fo the convective premixing scheme
  base_AGB - the base directory for AGB simulations
  
The run_star_extras.f files include axion energy-loss with the effects of electron degeneracy via the methods of G.G. Raffelt and D.S.P. Dearborn, Bounds on Hadronic
Axions From Stellar Evolution, Phys. Rev. D 36 (1987) 2211. These are included via additional terms to the neutrino energy-loss module. They also include some features from the MIST input physics J. Choi, A. Dotter, C. Conroy, M. Cantiello, B. Paxton and B.D. Johnson, Mesa Isochrones and Stellar Tracks (MIST). I. Solar-scaled Models, Astrophys. J. 823 (2016) 102 [1604.08592].
