## An extreme heatwave stochastic weather generator based on analogues
## and a simplified "importance sampling" algorithm
## Diagnostics on the simulations
## Pascal Yiou (LSCE), April 2018
SI=Sys.info()
if(SI[[1]] == "Darwin"){
##  Rsource="/Users/yiou/programmes/RStat/"
##  ANAdir="/Users/yiou/data/EUCLEIA/"
  OUTdir="/Users/yiou/data/EUCLEIA/"
}
if(SI[[1]] == "Linux"){
##  Rsource="/home/users/yiou/RStat/"
##  ANAdir="/home/estimr1/yiou/estimr1/NCEP/"
##  TNdir = "/home/estimr1/yiou/estimr1/WEGE/"
##  Tdir = "/home/estimr1/yiou/estimr1/ECAD/"
  OUTdir="/home/estimr2/yiou/IMPSAMP/" ## Needs to be user adjusted
}

library(parallel)

args=(commandArgs(TRUE))
print(args)
if(length(args)>0){
   fname=args[1]
}else{
    fname="TG-Orly-m6d1_extrHW-wege_cal0.3_TG0.5-276957.Rdat"
 }

setwd(OUTdir)

##save(file=fname,simu.dyn,alpha.cal,alpha.TN,simu.sta,args)
load(fname)
staname =args[1]
varname =args[2]
Lsim =as.integer(args[3])
mo.start =as.integer(args[4])
day.start =as.integer(args[5])
nsim =as.integer(args[6])
yymin=as.integer(args[7])
yymax=as.integer(args[8])
fileanalo=args[9]
jobid=args[10]

filout=paste(varname,"_m",mo.start,"d",day.start,
             "_extrHW_cal",alpha.cal,"_alphaX",alpha.TN,"_",
             jobid,".pdf",sep="")
if(!file.exists(filout)){
    pdf(filout)
    par(mar=c(4,5,1,1))
##    expression(paste("Temperature [",degree,"C]")),
    rangeplot=round(range(c(unlist(simu.sta$l.T.mean),
                            unlist(simu.dyn$l.X.mean)),na.rm=TRUE))
    ylab=expression(paste("TG (",degree,"C)"))
    plot(c(yymin,yymax),rangeplot,type="n",xlab="Years",
         ylab=ylab)
    boxplot(simu.sta$l.X.mean,at=c(yymin:yymax)-0.2,
            add=TRUE,axes=FALSE,col="blue")
    boxplot(simu.dyn$l.X.mean,at=c(yymin:yymax)+0.2,
            add=TRUE,axes=FALSE,col="red")
    lines(c(yymin:yymax),unlist(simu.sta$l.T.mean))
    legend("bottomright",lwd=c(1,5,5),col=c("black","blue","red"),
           legend=c("Obs.","Static","Dynamic"),bty="n")
    legend("bottomleft",
           legend=paste("t start = ",mo.start,"/",day.start),
           bty="n")
    legend("topleft","c",bty="n")
    dev.off()
}

q("no")
