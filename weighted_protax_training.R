#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message("Welcome to weighted_protax_training.R")
message("")
message("Step 1: Load necessary packages etc")
message("")
protaxdir<-args[1]
taxon<-args[2]
loci<-list.files(pattern=paste0(taxon,".database."))
for(i in 1:length(loci)){
  loci[i]<-sub(".fa","",loci[i])
  loci[i]<-sub(paste0(taxon,".database."),"",loci[i])
}

source(paste0(protaxdir,"/amcmc.rcode.txt"))
library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)
num.params=1+4
ind=1001:2000
message("")
message("Step 2: Run four iterations of training for each marker")
message("This will take some time...")
message("")

dirs = list.dirs(getwd(), recursive = F)

sel.dirs = NULL
for(i in c("12S","16S","Cytb")){
  sel.dir = sort(basename(dirs), decreasing=T)[grepl(i, sort(basename(dirs), decreasing=T), ignore.case = T)][1]
  sel.dirs = c(sel.dirs, sel.dir)
  sel.dirs = sel.dirs[!is.na(sel.dirs)]
}

for(locus in c("12S","16S","Cytb")){
  dir = sel.dirs[grepl(locus, sel.dirs, ignore.case = T)]
  if(length(dir) > 0){
    folder = paste0("./",sel.dirs[grepl(locus, sel.dirs, ignore.case = T)],"/")
    message(paste0("\nWorking on ",locus," in folder ",folder))

    for(level in c(1,2,3,4)){
      message(paste0("Working on level",level))
      dat=read.xdata(paste0(folder,"train.w",level,".xdat"))
      ppa=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
      initstate=initialize.adaptation(ppa$params[2000,])
      ppb=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
      initstate=initialize.adaptation(ppb$params[2000,])
      ppc=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
      initstate=initialize.adaptation(ppc$params[2000,])
      ppd=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
      pdf(paste0(folder,"weighted_training_plot_",locus,"_level",level,"a_MCMC.pdf"))
      traceplot.all(ppa,ind,num.levels=1, title="iter1")
      amcmc.diagnostic.plot(ppa)
      dev.off()
      pdf(paste0(folder,"weighted_training_plot_",locus,"_level",level,"b_MCMC.pdf"))
      traceplot.all(ppb,ind,num.levels=1, title="iter2")
      amcmc.diagnostic.plot(ppb)
      dev.off()
      pdf(paste0(folder,"weighted_training_plot_",locus,"_level",level,"c_MCMC.pdf"))
      traceplot.all(ppc,ind,num.levels=1, title="iter3")
      amcmc.diagnostic.plot(ppc)
      dev.off()
      pdf(paste0(folder,"weighted_training_plot_",locus,"_level",level,"d_MCMC.pdf"))
      traceplot.all(ppd,ind,num.levels=1, title="iter4")
      amcmc.diagnostic.plot(ppd)
      dev.off()
      k=which.max(ppa$postli[ind])
      write.postparams(ppa,paste0(folder,"w_mcmc",level,"a"),ind[k])
      k=which.max(ppb$postli[ind])
      write.postparams(ppb,paste0(folder,"w_mcmc",level,"b"),ind[k])
      k=which.max(ppc$postli[ind])
      write.postparams(ppc,paste0(folder,"w_mcmc",level,"c"),ind[k])
      k=which.max(ppd$postli[ind])
      write.postparams(ppd,paste0(folder,"w_mcmc",level,"d"),ind[k])
    }
  } 
  else {
    message(paste0("\nno results for ", locus)) 
    next
  }
}
message("")
message("End of weighted_protax_training.R")
q(save="no")
