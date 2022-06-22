# UNCLE-with-Analog-SWG
Importance Sampling Algorithm for Analog SWG to simulate extreme heatwaves
Analog SWG with importance sampling to simulate unprecedented climate
events (UNCLEs).
By Pascal Yiou (LSCE), April 2022.

This suite of programs enables to simulate extreme heatwaves with a stochastic
weather generator (SWG) with analogs of circulation, and importance sampling.

Analogs of circulation have to be computed first (e.g. with the "blackswan"
Web Processing Service): here $fileanalo.

A time series of temperature (or precipitation) is required in input for the
importance sampling (here $fileX)

A sample shell script (UNCLE_AnImSa_FR-warm12.sh) launches the
computations.
It contains the parameters of the simulations, which can be adjusted:
##
varname=TG
Lsim=62 ## Simulations have Lsim days
daystart=1
nsim=100 ## There are 100 simulations each year
yy0=1951 ## Annee de début
#weighTN=1 -> alphaTN=0.1 Poids sur les rangs de T
#weighTN=2 -> alphaTN=0.1 Poids sur les rangs de T, analogues de l'année de départ exclues
weightTN=2
alphaTN=0.5  ## Poids pour la selection des analogues extremes
alphacal=3 ## Poids pour la selection du jour calendaire
mostart=12
yy1=2020 ## Annee de fin
imax=1 ## Temperatures chaudes
fileX=era5_t2m_daily_fr.txt
fileanalo=ana_z500_1_-20.30.30.70.txt
##
Those variables serve as input to: UNCLE_AnSWGEN_v1.R
This R script requires to install the parallel package (with:
install.packages("parallel")). By default, the script run on 12 cores.
This can be lowered to fewer cores if the program runs on a "regular" laptop
computer.
You also need to install the lubridate package with install.packages("lubridate").

The output is an Rdata file, which contains the ensemble simulations and time averages. A sample output figure with boxplots is created with: HWgen_AnImSa_diags.R

Content:
UNCLE_AnImSa_README.txt: this file
UNCLE_AnImSa_FR-warm12.sh: shell (sh) script wrapper
UNCLE_AnSWGEN_v1.R: R script for simulation
HWgen_AnImSa_diags.R: R script for simple diagnostics
era5_t2m_daily_fr.txt: temperature daily averages over France from the ERA5 reanalysis
ana_z500_1_-20.30.30.70.txt: daily analogs of circulation from Z500 ERA5.

The paper to cite is:
Yiou, P. and Jézéquel, A.: Simulation of extreme heat waves with empirical importance sampling, Geosci. Model Dev., 13, 763–781, https://doi.org/10.5194/gmd-13-763-2020,  2020.

The code and data are provided "as is". They are distributed under a
free CECILL license (http://www.cecill.info/) for academic usage.
Please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr) for any commercial use.
