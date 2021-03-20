This code reproduces the following paper:

**Drivers of ESBL-producing *Escherichia coli* dynamics in calf fattening farms: A modelling study**

J. Bastard, M. Haenni, E. Gay, P. Glaser, J.Y. Madec, L. Temime, L. Opatowski

* **data_preparation.R** imports and prepares the data of ESBL-EC carriage and antimicrobial use
* **mcmc_estim_diagnost.R** performs the MCMC estimation using .bug codes from the "BUGS" folder, the chains diagnostic, posterior visualization and DIC calculation
* **simul_mod.R** defines functions for model simulations
* **mitig_strat** performs the simulation study on mitigation strategies
