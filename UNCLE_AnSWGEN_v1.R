## An extreme heatwave/cold spell stochastic weather generator
## based on analogues and a simplified "importance sampling" algorithm
## Pascal Yiou, Camille Cadiou, Flavio Pons (LSCE), June 2016, Feb 2018,
## June 2019, November 2019
## Revisions Novembre 2021, Janvier 2022, Avril 2022 (Camille Cadiou)
## Revisions June 2022 (Flavio Pons)
## This code is distributed freely and as is under a CeCill License.
## It is available in github
## Exemple. Se lance sur obelix par:
## qsub -q mediump -l nodes=1:ppn=12 ${HOME}/programmes/RStat/A2C2/IMPSAM/UNCLEgen_AnImSa.sh
## or
## qsub -q mediump -l nodes=1:ppn=12 ${HOME}/programmes/RStat/A2C2/IMPSAM/HWgen_AnImSa_demo.sh (provided on github)

## version ne prenant pas le jour de l'événement simulé mais ne filtrant pas les analogues de l'année en question.

## Computer path configurations
SI=Sys.info()
if(SI[[1]] == "Linux"){
  Tdir = "/home/users/ccadiou/Data/ERA5/" ## Needs to be adapted
  OUTdir="/home/estimr3/ccadiou/WEGE/" ## Needs to be adapted
}
if(SI[[1]] == "Darwin"){
  Tdir = "/home/users/ccadiou/Data/ERA5/" ## Needs to be adapted
  OUTdir="/home/estimr3/ccadiou/WEGE/" ## Needs to be adapted
}

library(parallel)
library(lubridate) # Required to handle dates correctly, including leap years

##R CMD BATCH "--args ${fileX} ${varname} ${Lsim} ${mostart} ${daystart} ${nsim} ${yy0} ${yy1} ${fileanalo} ${JOBID} ${alphaTN} ${weightT} ${alphacal} ${imax}" /home/users/yiou/programmes/RStat/A2C2/IMPSAM/UNCLE_AnImSa-gen.R /home/users/yiou/programmes/RStat/A2C2/IMPSAM/UNCLE_AnImSa${jobid}.Rout

## Arguments of the script
args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    fileX = args[i];i=i+1 # Fichier de variable a optimiser
    varname =args[i];i=i+1 # Nom de la variable
    Lsim =as.integer(args[i]);i=i+1 # Duree des simulations (e.g. une saison = 90j)
    mo.start =as.integer(args[i]);i=i+1 # Premier mois
    day.start =as.integer(args[i]);i=i+1 # Premier jour
    nsim =as.integer(args[i]);i=i+1 # Nombre de simulations
    yymin=as.integer(args[i]);i=i+1 # Premiere annee de simu
    yymax=as.integer(args[i]);i=i+1 # Derniere annee de simu
    fileanalo=args[i];i=i+1 # Fichier d'analogues
    jobid=args[i];i=i+1 # No de job
    alpha.TN <<- as.numeric(args[i]);i=i+1 # Poids sur le rang des temperatures
    weight.T = args[i];i=i+1 # Type de poids (par defaut: 1, sinon 2, 3 ou 4)
    alpha.cal <<- as.numeric(args[i]);i=i+1 # Poids sur les jours calendaires
##    i.max=as.numeric(args[i]);i=i+1
    i.max = ifelse(as.numeric(args[i])==1,TRUE,FALSE);i=i+1 # Valeurs elevees (i.max=TRUE) ou faibles (i.max=FALSE) de la variable
## R n'aime pas les arguments booleens, donc args[i] est 1 (vdc) ou 0 (vdf)
}else{ ## Default options
    fileX="/home/estimr2/yiou/IMPSAMP/tmax_daily_era5_WPAC.txt"
    varname="TX"
    Lsim = 30
    mo.start =06
    day.start =15
    nsim =100
    yymin=2000
    yymax=2021
    fileanalo="http://birdy.lsce.ipsl.fr:8096/outputs/d7b669a8-5287-11ec-9524-00304844f2cc/output.txt"
    jobid="test"
    alpha.TN <<- 0.5
    weight.T = 1
    alpha.cal <<- 1
    i.max=TRUE
}

if(!(weight.T %in% c(1,2))){
    print("weight.T should be in c(1,2)")
    q("no")
}

## ## Sets days in calendar year
## ## Creates l.mo, l.da et moda (list of days in year)
## l.mo=c(1:12)
## l.da=c(31,28,31,30,31,30,31,31,30,31,30,31)
## moda=c()
## for(i in l.mo){
##   for(j in 1:l.da[i]){
##     moda=c(moda,i*100+j)
##   }
## }
## ------------------------------------------------------------------------
## Read input data
## Read analog file
##example: fileanalo="http://www-lscedods.cea.fr/estimr2/DASE_NK/ana_slp_surface_base_rms_NA_latest_-80.0_50.0_22.5_70.0_1_30_20.txt"
analo = read.table(fileanalo,header=TRUE)

## date.a = analo$date
## date.a.cal=match(as.integer(substr(analo$date,5,8)),moda)

## Read time series to be optimized in the analog importance sampling
## Adjust options, for the data format in fileX
setwd(Tdir)
X.dat = read.table(file=fileX,header=TRUE,skip=2)
XX.date = X.dat[,1]
XX = X.dat[,2]

## ------------------------------------------------------------------------
## Computer parameters, especially for running in // on the LSCE batch cluster
## Please adapt to your own cluster
ncpus = as.numeric(Sys.getenv(c("NCPU")))

print(paste(detectCores(),"cores detected"))
nproc=3
if(SI[[1]] == "Darwin") nproc=2 ## Pour le mac portable
if(SI[[4]] %in% paste("obelix",2:6,sep="")) nproc=3 # Pour les machines interactives
if(SI[[4]] %in% paste("obelix",10:51,sep="")){
  nproc=  as.numeric(Sys.getenv(c("NCPU"))) # Pour les machines BATCH du LSCE
}
ncpus = max(nproc,ncpus,na.rm=TRUE)
print(paste("Calcul sur",ncpus,"CPUs"))

## ------------------------------------------------------------------------
## Simulation stochastique avec des poids sur la distance au jour calendaire
## a simuler, et le rang de la temperature.
## lanamax est le nombre de jours analogues autour du jour cible
## Version 1: poids sur les rangs des temperatures analogues
## Prefered method
"simu.extrHW1" = function(t.start=20030601,Lsim=90,alpha.cal = 0.5,
                          alpha.TN = 0.1,lanamax=30,i.max=i.max)
{
    I0 = which(analo$date == t.start)
##    t0.cal = match(as.integer(substr(t.start,5,8)),moda)
## Modification par F. Pons (2022/21/04), correctly accounts for leap years 
    ## t0.cal= yday(as.Date(paste(substr(t.start, start = 1, stop = 4),
    ##                            substr(t.start, start = 5, stop = 6), 
    ##                            substr(t.start, start = 7, stop = 8),
    ##                            sep = '-')))
    t0.cal=yday(as.Date(as.character(t.start), format='%Y%m%d'))
    t0=t.start
    T.sim=c(XX[XX.date == t0])
    t.sim=c(t0)
    ndum=c()
    for(i in 1:Lsim){
        I = which(analo$date == t0)
      I1 = I +1
      t1=analo$date[I1]
      ana.d1 = c(analo$date[I0+i],unlist(analo[I1,2:21]))
      ana.d1 = intersect(ana.d1,XX.date)
## Poids relatif au jour calendaire pour respecter le cycle saisonnier  
##        ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
## Modification par F. Pons (2022/21/04), correctly accounts for leap years   
        ## ana.d1.cal=yday(as.Date(paste(substr(ana.d1, start = 1, stop = 4),
        ##                               substr(ana.d1, start = 5, stop = 6), 
        ##                 substr(ana.d1, start = 7, stop = 8), sep = '-'))) 
       ana.d1.cal=yday(as.Date(as.character(ana.d1), format='%Y%m%d'))
       
## Modification par C. Cadiou (2021/11/24), pour tenir compte du saut
## entre decembre et janvier
      diff.ana.cal.1 = abs(ana.d1.cal - ((t0.cal+i)%%366)) ## "-1" enlevé à la fin
      diff.ana.cal.2 = 365 - diff.ana.cal.1
      diff.ana.cal = pmin(diff.ana.cal.1,diff.ana.cal.2)
      weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
      sdum=sum(weight.cal,na.rm=TRUE)
      weight.cal = weight.cal/sdum
## Poids sur la valeur du rang de la temperature (pour avoir l'analogue le
## plus chaud ou le plus froid, selon la valeur de i.max)  
      d1=ana.d1
      Idum = match(ana.d1,XX.date)
      TN.d1 = XX[Idum] ## Temperatures analogues
      TN.d1.sort=sort(TN.d1,index.return=TRUE,decreasing=i.max,
                      na.last=TRUE,method="radix")
      weight.TN = exp(-alpha.TN*c(1:length(TN.d1)))
## Nombre d'analogues = length(TN.d1)
      sdum=sum(weight.TN,na.rm=TRUE)
      weight.TN[TN.d1.sort$ix] = weight.TN/sdum
## Produit des poids pour normalisation des probabilites
      weight.all=weight.cal * weight.TN
      weight.all[is.na(weight.all)]=0
      weight.all = weight.all/sum(weight.all,na.rm=TRUE)
      ndum=c(ndum,length(which(weight.all >= 0.05)))
## Echantillonnage aleatoire selon la combinaison des poids     
      d1.max = sample(d1,size=1,prob=weight.all)
      T.sim = c(T.sim,XX[XX.date == d1.max])
      t.sim = c(t.sim, d1.max)
      t0=ifelse(idyn0,d1.max,t1)
    }
    return(cbind(t.sim,T.sim))
}
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
## Simulation stochastique avec des poids sur la distance au jour calendaire
## a simuler, et le rang de la temperature.
## lanamax est le nombre de jours analogues autour du jour cible
## Version 2: poids sur les rangs des temperatures analogues (comme v1)
## mais on exclut les analogues de l'annee de t.start
## Prefered method
"simu.extrHW2" = function(t.start=20030601,Lsim=90,alpha.cal = 0.5,
                          alpha.TN = 0.1,lanamax=30,i.max=i.max)
{
    I0 = which(analo$date == t.start)
##    t0.cal = match(as.integer(substr(t.start,5,8)),moda)
## Modification par F. Pons (2022/21/04), correctly accounts for leap years 
    ## t0.cal <- yday(as.Date(paste(substr(t.start, start = 1, stop = 4),
    ##                              substr(t.start, start = 5, stop = 6), 
    ##                              substr(t.start, start = 7, stop = 8),
    ##                              sep = '-'))) 
    t0.cal=yday(as.Date(as.character(t.start), format='%Y%m%d')) 
    t0=t.start
    yy0=floor(t0/10000)
    T.sim=c(XX[XX.date == t0])
    t.sim=c(t0)
    ndum=c()
    for(i in 1:Lsim){
      I = which(analo$date == t0)
      I1 = I +1
      t1=analo$date[I1]
      ana.d1 = unlist(analo[I1,2:21])
      ana.d1 = intersect(ana.d1,XX.date)
## Poids relatif au jour calendaire pour respecter le cycle saisonnier  
## Modification par F. Pons (2022/21/04), correctly accounts for leap years    
      ana.d1.cal=yday(as.Date(as.character(ana.d1), format='%Y%m%d'))
      ## ana.d1.cal=yday(as.Date(paste(substr(ana.d1, start = 1, stop = 4),
      ##                               substr(ana.d1, start = 5, stop = 6), 
      ##                               substr(ana.d1, start = 7, stop = 8),
      ##                               sep = '-'))) 
##      ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
## Modification par C. Cadiou (2021/11/24), pour tenir compte du saut
## entre decembre et janvier
      diff.ana.cal.1 = abs(ana.d1.cal - ((t0.cal+i)%%366)) ## "-1" enlevé à la fin
      diff.ana.cal.2 = 365 - diff.ana.cal.1
      diff.ana.cal = pmin(diff.ana.cal.1,diff.ana.cal.2)
      weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
      sdum=sum(weight.cal,na.rm=TRUE)
      weight.cal = weight.cal/sdum
## Poids sur la valeur du rang de la temperature (pour avoir l'analogue le
## plus chaud ou le plus froid, selon la valeur de i.max)  
      d1=ana.d1
      Idum = match(ana.d1,XX.date)
      TN.d1 = XX[Idum] ## Temperatures analogues
      TN.d1.sort=sort(TN.d1,index.return=TRUE,decreasing=i.max,
                      na.last=TRUE,method="radix")
      weight.TN = exp(-alpha.TN*c(1:length(TN.d1)))
## Nombre d'analogues = length(TN.d1)
      sdum=sum(weight.TN,na.rm=TRUE)
      weight.TN[TN.d1.sort$ix] = weight.TN/sdum
## Produit des poids pour normalisation des probabilites
      weight.all=weight.cal * weight.TN
      weight.all[is.na(weight.all)]=0
## On met a 0 les poids des analogues qui sont dans l'événement observé
## Modification par C. Cadiou (09/06/2022)
      obs.dates = as.numeric(format(seq(as.Date(as.character(t.start),
                                                format="%Y%m%d"),
                                        by="day",length.out=Lsim),"%Y%m%d")) # vecteur des dates de l'événement observé
      weight.all[ana.d1 %in% obs.dates] = 0  # si une analogue est dans ce vecteur, on met son poids à 0
## Normalisation pour avoir une somme a 1      
      weight.all = weight.all/sum(weight.all,na.rm=TRUE)
      ndum=c(ndum,length(which(weight.all >= 0.05)))
## Echantillonnage aleatoire selon la combinaison des poids     
      d1.max = sample(d1,size=1,prob=weight.all)
      T.sim = c(T.sim,XX[XX.date == d1.max])
      t.sim = c(t.sim, d1.max)
      t0=ifelse(idyn0,d1.max,t1)
    }
    return(cbind(t.sim,T.sim))
}
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Definition de la fonction de simulation
fun.name = paste("simu.extrHW",weight.T,sep="")
print(paste("Applying",fun.name))
SIMU.FUNC = match.fun(fun.name)

## ------------------------------------------------------------------------
## Wrapper pour les simulations stochastiques
"wrap.extrHW" = function(k)
  {
      XX = SIMU.FUNC(t.start=t.start,Lsim=Lsim,
                       alpha.TN=alpha.TN,alpha.cal=alpha.cal,i.max=i.max)
    return(XX)
  }
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Wrapper pour le calcul des moyennes de temperature
"wrap.mean" = function(i)
  {
    mm = mean(Xdum[[i]][,2])
    return(mm)
  }
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Simulations pour la reanalyse
"simu.yy" = function(idyn=TRUE,yymin=1948,yymax=2018,
                     mo0=mo.start,day0=day.start,i.max=i.max)
{
    idyn0 <<- idyn
    l.X.mean.dyn=list()
    l.T.mean.dyn=list()
    l.X.dyn = list()
    for(yy in yymin:yymax){
        print(paste("Processing",yy))
        t.start <<- (yy*100+mo0)*100+day0
## Modification par C. Cadiou (26/04/2022) pour considérer l'hiver : si la simulation est à cheval sur deux ans, on l'associe à la deuxième année
##        overlap <<- (match(as.integer(substr(t.start,5,8)),moda)+Lsim>365)
## overlap <<- (yday(as.Date(paste(substr(t.start,start = 1, stop = 4),
##                                 substr(t.start, start = 5, stop = 6),
##                                 substr(t.start, start = 7, stop = 8),
##                                 sep = '-'))) + Lsim > 365)
        overlap <<- (yday(as.Date(as.character(t.start), format='%Y%m%d'))
            + Lsim > 365)
##      yy.sim <<- ifelse(overlap,yy,yy+1)
        yy.sim <<- ifelse(overlap,yy+1,yy)
        Xdum <<- mclapply(seq(1,nsim,by=1),wrap.extrHW,mc.cores=ncpus)
        l.X.dyn[[as.character(yy.sim)]]=Xdum
        X.mean = mclapply(seq(1,nsim,by=1),wrap.mean,mc.cores=ncpus)
        l.X.mean.dyn[[as.character(yy.sim)]]=unlist(X.mean)
        iTN=which(XX.date==t.start)
        l.T.mean.dyn[[as.character(yy.sim)]]=mean(XX[c(iTN:(iTN+Lsim))],
                                              na.rm=TRUE)
    }
    return(list(l.X.mean=l.X.mean.dyn,
                l.T.mean=l.T.mean.dyn,
                l.X=l.X.dyn,ymin=yymin,ymax=yymax))
}
## ------------------------------------------------------------------------
simu.dyn=simu.yy(idyn=TRUE,yymin=yymin,yymax=yymax,i.max=i.max)
simu.sta=simu.yy(idyn=FALSE,yymin=yymin,yymax=yymax,i.max=i.max)

setwd(OUTdir)
uncle.type=ifelse(i.max,"max","min")
fname=paste(varname,"-m",mo.start,"d",day.start,"L",Lsim,
            "_UNCLE-",uncle.type,"-animsa_cal",alpha.cal,"_",varname,
            alpha.TN,"meth",weight.T,"-",jobid,".Rdat",sep="")

save(file=fname,simu.dyn,alpha.cal,alpha.TN,simu.sta,args,weight.T,i.max)

q("no")
## END of simulation SCRIPT
## ------------------------------------------------------------------------


