**Guidelines to measure cophylogenetic signal and phylogenetic congruence in a bipartite interaction network:**
====

<br>

This folder contains the scripts for reproducing the analysis of *Distinguishing cophylogenetic signal from phylogenetic congruence clarifies the interplay between evolutionary history and species interactions* by Benoît Perez-Lamarque and Hélène Morlon. (citation: <a href="https://doi.org/10.1093/sysbio/syae013">Perez-Lamarque & Morlon, 2024; Systematic Biology</a>).


These scripts can be used to mesure two different patterns: 
1) **cophylogenetic signal** ("do closely related hosts interact with closely related symbionts?")
2) **phylogenetic congruence** ("do host and symbiont phylogenies mirror one another?").

<br>

# Contents:

It contains 3 files:
- "script_cophylogenetic_signal.R": the script to reproduce the simulations, run global-fit approaches for measuring cophylogenetic signal, and plot the figures
- "functions_parafit_paco.R": ParaFit and PACo functions that implement the  "null model 2" that randomly shuffles species labels in one clade (see https://doi.org/10.1111/mec.16681 for more details)
- "script_run_eMPRess.sh": the script for automatically running eMPRess with the command line on our simulations for measuring phylogenetic congruence. 


Contact: benoit.perez-lamarque@utoulouse.fr

<br> 

# Guidelines for measuring cophylogenetic signals in species interaction network:

XXX


<br> 


## Assessing the effect of biogeography on the cophylogenetic signal in species interactions

This can be using biogeography-based permutations that only shuffle species that co-occur on the same area. <a href="https://github.com/BPerezLamarque/Scripts/tree/master/Biogeography_interactions">See the complete tutorial here.</a> (citation: Ibos et al., 2026, in prep.) 



<br> <br>