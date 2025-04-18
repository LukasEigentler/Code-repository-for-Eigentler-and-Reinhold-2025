# Code-repository 
This repository contains code associated with the paper Names redacted (2023) DOI: TBC.
All code is provided through GNU General Public License v3.0. Please see the license file for more information

"hopf_locus_numsim_unbounded.m" calculates the location of the transition between stable and oscillatory solutions in a two-dimensional parameter plane.

"mutation_selection_balance_unbounded.m" calculates how much individual variability occurs within the prey population in the absence of predators.

"pred_prey_prey_defence_ode_unbounded.m" is a functiont that contains the expressions of the ODEs that are fed into the ode solver.

"predation_pressure_function_plots.m" visualises the predation pressure on each genotype depending on the trait distribution of the whole prey population.

"prey_defence_num_sim_para_changes_unbounded.m" calculates and visualises the bifurcation diagrams.

"prey_defence_single_run_unbounded.m" is a script that simulates the model once.

"prey_defence_single_run_fun_unbounded.m" is a function that simulates the model once and calculates output quantities of interest.

"trait_bound_test_single_para_set.m" is a script that fixes all parameters and performs simulations for this parameter set for different sizes of the trait domain.

All files named "filename_unbounded.m" have a corresponding "filename.m". These correspond to to model simulations for the model with trait domain [0,alpha_1^{-1}]. Note, however, that such simulations could also be achieved through the "filename_unbounded.m" files by changing the domain bound parameters to the required values.


