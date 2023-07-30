rm(list=ls())

library(ggplot2)
library(gridExtra)
library(ggthemes)
library(xtable)
library(parallel)
library(VGAM)
library(BB)

ncores=5

dir.out='~/'
source(paste0(dir.out,"Design_Functions_bound_R1_1.R"))
source(paste0(dir.out,"design_bound_R1_1.R"))


# Run Normal design for given prior ---------------------------------------

oc_normal=designF(index=1,outcome='normal',prior.pars=c(0.25,sqrt(1/50)),
            designGrid=matrix(c(0.25,sqrt(1/50)),1,2),th0=0,
            w.list=seq(0,1,by=0.01),thR=0.15,n.list=1:250,sigma=1,
            alpha.up=0.15,cores=ncores)
oc=oc_normal
save(oc, file=paste0(dir.out,"normal_bound_R1.RData"))


# Run sensitivity analysis for Normal design ------------------------------


oc_normal_sensitivity=mclapply(1:21, function(id) 
                                      designF(index=id,outcome='normal',
                                      prior.pars=c(0.25,sqrt(1/50)),
                                      designGrid=matrix(c(seq(-0.7,1.3,by=0.1),
                                                          rep(sqrt(1/50),21)),21,2),
                                      th0=0,thR=0.15,sigma=1,alpha.up=0.15,
                                      cores=1,n.list=1:250),
                                      mc.cores = ncores)
oc=oc_normal_sensitivity
save(oc, file=paste0(dir.out,"normal_sens_bound_R1.RData"))


# Run Binomial design for given prior -------------------------------------


oc_binomial=designF(index=1,outcome='binomial',prior.pars=c(20,20),
                    designGrid=matrix(c(20,20),1,2),
                    th0=0.3,thR=0.5, n.list=1:250,alpha.up=0.15, 
                    w.list=seq(0,1,by=0.01),cores=ncores)
oc=oc_binomial
save(oc, file=paste0(dir.out,"binomial_bound_R1.RData"))


# Run sensitivity analysis for binomial design ----------------------------


oc_binomial_sensitivity=mclapply(1:21, function(id) 
                                      designF(index=id,w.list=0.5,
                                      outcome='binomial',prior.pars=c(20,20),
                                      designGrid=matrix(c(seq(0,40,by=2),
                                                          (40-seq(0,40,by=2))),22,2),
                                      th0=0.3,n.list=1:250,thR=0.5,cores=1,
                                      alpha.up=0.15), mc.cores = ncores)
oc=oc_binomial_sensitivity
save(oc, file=paste0(dir.out,"binomial_sens_bound_R1.RData",sep=""))


# Plots -------------------------------------------------------------------

source(paste0(dir.out,"Graphs_R1.R"))

#Designs & Sensitivity analyses
plots.design(outcome='normal',dir.out=dir.out)
plots.design(outcome='binomial',dir.out=dir.out)

#Frequentist risk
plot.fr(dir.out=dir.out)

# Case study results and plots ---------------------------------------------

source(paste0(dir.out,"case_study.R"))
designRunDataEx(dir.out=dir.out,alpha.up=0.15,bound=TRUE)
designRunDataEx(dir.out=dir.out,alpha.up=1,bound=FALSE)

