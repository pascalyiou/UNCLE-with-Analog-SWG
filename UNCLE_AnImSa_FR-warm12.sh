#!/bin/sh -l
## Simulations with Analog Importance sampling of extreme events.
## Here the simulations start in December (=12) and end in January
## The optimization is done for continental France
## The code can be launched with:
## qsub -q shortp -l nodes=1:ppn=12 ${HOME}/programmes/RStat/A2C2/IMPSAM/UNCLE_AnImSa_FR-warm12.sh
## Pascal Yiou, LSCE, Apr. 2022

## Those instructions should be modified according to the working environment
module load R/4.0.3
##PBS_NODEFILE is the number of cores. It can be specified "manually"
NCPU=`wc -l < $PBS_NODEFILE` 
export NCPU

JOBID=`echo ${PBS_JOBID} | cut -d. -f1`

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nStarting script at: ${start_date}\n"
echo -e ${JOBID}

varname=TG ## Name of the variable to be optimized
Lsim=62 ## Simulations have Lsim days
daystart=1 ## Day of the month
nsim=100 ## There are 100 simulations each year
yy0=1951 ## Stating year

## weighTN=1 -> alphaTN=0.1 Poids sur les rangs de T
## weighTN=2 -> alphaTN=0.1 Poids sur les rangs de T, analogues de l'année de départ exclues
weightTN=1
alphaTN=0.5  ## Poids pour la selection des analogues extremes

alphacal=3 ## Poids pour la selection du jour calendaire

## Winter heat waves
mostart=12
yy1=2020 ## Annee de fin
imax=1 ## Temperatures chaudes. imax=0 for cold spells

## Donnees de temperature
fileX=/home/users/ccadiou/Data/ERA5/era5_t2m_daily_fr.txt

## Analogues
## North Atl Z500 analogs 1951-2021 in 1951-2021
fileanalo=/home/users/ccadiou/Data/Analogues/ana_z500_1_-20.30.30.70.txt

jobid=${JOBID}-1951-2020

R CMD BATCH "--args ${fileX} ${varname} ${Lsim} ${mostart} ${daystart} ${nsim} ${yy0} ${yy1} ${fileanalo} ${jobid} ${alphaTN} ${weightTN} ${alphacal} ${imax}" ${HOME}/programmes/RStat/A2C2/IMPSAM/UNCLE_AnImSa_gen_v2.R ${HOME}/programmes/RStat/A2C2/IMPSAM/UNCLE_AnImSa_${jobid}.Rout

## Trace d'un diagnostic
## Attention "max" doit etre change en "min" si on calcule des valeurs min
fsim=${varname}_m${mostart}d${daystart}L${Lsim}_UNCLE-max-animsa_cal${alphacal}_alphaX${alphaTN}_meth${weightTN}_${jobid}.Rdat

R CMD BATCH "--args ${fsim}" ${HOME}/programmes/RStat/A2C2/IMPSAM/HWgen_AnImSa_diags.R ${HOME}/programmes/RStat/A2C2/IMPSAM/HWgen_AnImSa_diags.Rout

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nEnding script at: ${start_date}\n"
