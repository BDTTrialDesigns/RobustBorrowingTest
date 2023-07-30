#Calculate operating characteristics for the various approaches

designF = function(index,       #row of designGrid matrix, for sensitivity analyses
                   outcome,     #data outcome : normal, binomial
                   prior.pars,  #informative prior parameters: normal: c(mean,SD), binomial: c(a,b)
                   designGrid,  #matrix of sampling prior parameters. each row corresponds to one prior.
                   th0,         #point for type I error rate evaluation/upper boundary of null Hp
                   thR,         #point for conditional power evaluation
                   n.list=seq(1,250,by=1),  #sample size grid, for sample size selection
                   w.list=seq(0,1,by=0.05), #weights grid, for CD and Mixture approaches
                   sigma=NULL,  #data sd for one observation. normal outcomes only.
                   cores,       #num of cores for parallelisation
                   alpha.up=1,  #upper type I error rate boundary for CD-Adapt
                   plot=FALSE,  #plot of adaptive weight
                   S=8,         #number of standard deviations for integration grid (normal outcomes only)
                   case.study=FALSE, #case study example settings
                   bound=TRUE)  #EB power prior bounded?
{
  if (outcome=="normal")
  { 
    #priors
    true.pars=c(designGrid[index,1],designGrid[index,2])
    sampling = c(1, true.pars[1], true.pars[2])
    info <-c(1,prior.pars[1],prior.pars[2])
    base <-c(1,0,100)
    base.mix = c(prior.pars[1],sigma)
    base.mix2 = c(prior.pars[1],100)
    
  }
  
  if (outcome=="binomial")
  {
    #priors
    true.pars=designGrid[index,]
    sampling = c(1, 1+true.pars[1], 1+true.pars[2])
    info <-c(1,1+prior.pars[1], 1+prior.pars[2])
    n0=prior.pars[1]+prior.pars[2]
    base <-c(1,0.001,1)
    base.mix = c(1,1)
    base.mix2 = c(0.001,1)
  }
  
  #test error costs
  c.alpha = 0.975
  c.beta = 0.025
  alpha=c.beta/(c.beta+c.alpha)
  
  #OCs of approaches with varying weight
  Out_List= mclapply(w.list, function(wT) {
    
    if(!case.study) mix <-cbind(c(wT,info[2:3]), c(1-wT,base.mix))
    if(case.study) mix <-cbind(c(wT,info[2:3]), c(1-wT,c(0, sqrt(8.27)))) 
    mix2 <- cbind(c(wT,info[2:3]), c(1-wT, th0, 0))
    mix3 <- cbind(c(wT,info[2:3]), c(1-wT, base.mix2))
    prior.list = list(base, mix, mix2,mix3)
    names(prior.list)=c("vague","mix","mixRB","mixAlt")
    
    lp=lapply(1:length(prior.list),function(y) {
      
      sp=sapply(n.list,function(x) {
        
        if (names(prior.list)[y]=='vague')
        {
          #Compute Type I error rate under BD
          if (outcome=="normal")
          { 
            gamma=sigma^2/(x*info[3]^2)
            alpha1=1-pnorm(gamma*(th0-info[2])*sqrt(x)/sigma + qnorm(1-alpha) *sqrt(1+gamma))
          }
          if (outcome=="binomial")
          { 
            priorPars=prior.pars
            thrroot=function(g)
            {
              post.a=info[2]+g
              post.b=info[3]+(x-g)
              fun=(pbeta(th0,post.a,post.b))-alpha
              fun
            }
            
            xt=0:x
            pr=thrroot(xt)
            x_targ=xt[which((1-abs(pr))*sign(-pr)==max((1-abs(pr))*sign(-pr)))]
            
            post.a=base[2]+x_targ
            post.b=base[3]+(x-x_targ)
            alpha1=pbeta(th0,post.a,post.b)
          }
          #Set CD threshold
          alphaS=(1-wT)*alpha+wT*alpha1
        }
        else
        {
          alphaS=alpha
        }
        
        if(outcome=="normal")
        {
          oc=integrand_Err(prior=names(prior.list)[y],
                           n=x,alpha=alphaS,
                           dataPar=sigma,priorParsWeights=prior.list[[y]],
                           lik=outcome,sampPrPars=sampling[2:3],
                           th0=th0,S=S)
          type1=oc[1] #Expected type I error rate
          type2=oc[2] #Expected type II error rate
          risk=oc[3]  #Integrated risk
          
          typeIfr=integrand_Err(prior=names(prior.list)[y],
                                n=x,alpha=alphaS,
                                dataPar=sigma,priorParsWeights=prior.list[[y]],
                                lik=outcome,sampPrPars=c(th0,0),
                                th0=th0,S=S)[1] #Frequentist type I error rate at th0
          
          typeIIfr=integrand_Err(prior=names(prior.list)[y],
                                 n=x,alpha=alphaS,
                                 dataPar=sigma,priorParsWeights=prior.list[[y]],
                                 lik=outcome,sampPrPars=c(thR,0),
                                 th0=th0,S=S)[2] #Frequentist type II error rate at thR
        }
        
        if(outcome=="binomial")
        {
          oc=integrand_Err(prior=names(prior.list)[y],
                           n=x,n0=n0,alpha=alphaS,
                           priorParsWeights=prior.list[[y]],
                           lik=outcome,sampPrPars=sampling[2:3],
                           th0=th0)
          type1=oc[1] #Expected type I error rate
          type2=oc[2] #Expected type II error rate
          risk=oc[3]  #Integrated risk
          
          typeIfr=integrand_Err(prior=names(prior.list)[y],
                                n=x,n0=n0,alpha=alphaS,
                                priorParsWeights=prior.list[[y]],
                                lik=outcome,sampPrPars=c(th0,0),
                                th0=th0)[1]  #Frequentist type I error rate at th0
          
          typeIIfr=integrand_Err(prior=names(prior.list)[y],
                                 n=x,n0=n0,alpha=alphaS,
                                 priorParsWeights=prior.list[[y]],
                                 lik=outcome,sampPrPars=c(thR,0),
                                 th0=th0)[2] #Frequentist type II error rate at thR
        }
        
        
        
        
        out=c('type1fr'= typeIfr,
              'type2fr'=typeIIfr,
              'ErrSum'=risk,
              'type1'=type1,
              'type2'=type2)
        out
      })
      colnames(sp)=n.list
      sp
    })
    names(lp)=prior.list
    lp
  },mc.cores = cores)
  
  names(Out_List)= w.list
  
  
  #OCs of FD, BD, and EBPD
  prior.list = list(base, info, info)
  names(prior.list)=c("vague","info","power")
  
  Out_List_UW= lapply(1:length(prior.list),function(y) {
    
    sp=sapply(n.list,function(x) {
      
      if(outcome=="normal")
      {
        oc=integrand_Err(prior=names(prior.list)[y],
                         n=x,alpha=alpha,
                         dataPar=sigma,priorParsWeights=prior.list[[y]],
                         lik=outcome,sampPrPars=sampling[2:3],
                         th0=th0,S=S,bound=bound)
        type1=oc[1] #Expected type I error rate
        type2=oc[2] #Expected type II error rate
        risk=oc[3]  #Integrated risk
        
        typeIfr=integrand_Err(prior=names(prior.list)[y],
                              n=x,alpha=alpha,
                              dataPar=sigma,priorParsWeights=prior.list[[y]],
                              lik=outcome,sampPrPars=c(th0,0),
                              th0=th0,S=S,bound=bound)[1] #Frequentist type I error rate at th0
        
        typeIIfr=integrand_Err(prior=names(prior.list)[y],
                               n=x,alpha=alpha,
                               dataPar=sigma,priorParsWeights=prior.list[[y]],
                               lik=outcome,sampPrPars=c(thR,0),
                               th0=th0,S=S,bound=bound)[2] #Frequentist type II error rate at thR
      }
      
      
      if(outcome=="binomial")
      {
        oc=integrand_Err(prior=names(prior.list)[y],
                         n=x,n0=n0,alpha=alpha,
                         priorParsWeights=prior.list[[y]],
                         lik=outcome,sampPrPars=sampling[2:3],
                         th0=th0,bound=bound)
        type1=oc[1] #Expected type I error rate
        type2=oc[2] #Expected type II error rate
        risk=oc[3]  #Integrated risk
        
        typeIfr=integrand_Err(prior=names(prior.list)[y],
                              n=x,n0=n0,alpha=alpha,
                              priorParsWeights=prior.list[[y]],
                              lik=outcome,sampPrPars=c(th0,0),
                              th0=th0,bound=bound)[1] #Frequentist type I error rate at th0
        
        typeIIfr=integrand_Err(prior=names(prior.list)[y],
                               n=x,n0=n0,alpha=alpha,
                               priorParsWeights=prior.list[[y]],
                               lik=outcome,sampPrPars=c(thR,0),
                               th0=th0,bound=bound)[2] #Frequentist type II error rate at thR
      }
      
      
      
      
      out=c('type1fr'= typeIfr,
            'type2fr'=typeIIfr,
            'ErrSum'=risk,
            'type1'=type1,
            'type2'=type2)
      out
    })
    colnames(sp)=n.list
    sp
  })
  names(Out_List_UW)=c("Vague","Info","Power")
  
  #OCs of CD-Adapt
  Out_List_Wad= sapply(n.list,function(x) {
    
    if(outcome=="normal")
    {
      gamma=sigma^2/(x*prior.pars[2]^2)
      alpha1=1-pnorm(gamma*(th0-prior.pars[1])*sqrt(x)/sigma + qnorm(1-alpha) *sqrt(1+gamma))
      
      
      oc=integrand_Err(prior='vague',
                       n=x,alpha=alpha,priorParsInfo=info,
                       dataPar=sigma,priorParsWeights=base,
                       lik=outcome,sampPrPars=sampling[2:3],
                       th0=th0,w_adapt=TRUE,alpha1=alpha1,alpha.up=alpha.up,plot=plot,S=S)
      type1=oc[1] #Expected type I error rate
      type2=oc[2] #Expected type II error rate
      risk=oc[3]  #Integrated risk
      
      typeIfr=integrand_Err(prior='vague',
                            n=x,alpha=alpha,priorParsInfo=info,
                            dataPar=sigma,priorParsWeights=base,
                            lik=outcome,sampPrPars=c(th0,0),
                            th0=th0,w_adapt=TRUE,alpha1=alpha1,alpha.up=alpha.up,S=S)[1] #Frequentist type I error rate at th0
      typeIIfr=integrand_Err(prior='vague',
                             n=x,alpha=alpha,priorParsInfo=info,
                             dataPar=sigma,priorParsWeights=base,
                             lik=outcome,sampPrPars=c(thR,0),
                             th0=th0,w_adapt=TRUE,alpha1=alpha1,alpha.up=alpha.up,S=S)[2] #Frequentist type II error rate at thR
      
    }
    
    if(outcome=="binomial")
    {
      thrroot=function(g)
      {
        post.a=info[2]+g
        post.b=info[3]+(x-g)
        fun=(pbeta(th0,post.a,post.b))-alpha
        fun
      }
      
      xt=0:x
      pr=thrroot(xt)
      x_targ=xt[which((1-abs(pr))*sign(-pr)==max((1-abs(pr))*sign(-pr)))]
      
      post.a=base[2]+x_targ
      post.b=base[3]+(x-x_targ)
      alpha1=pbeta(th0,post.a,post.b)
      
      oc=integrand_Err(prior='vague',
                       n=x,n0=n0,alpha=alpha,priorParsInfo=info,
                       priorParsWeights=base,
                       lik=outcome,sampPrPars=sampling[2:3],
                       th0=th0,w_adapt=TRUE,alpha1=alpha1,alpha.up=alpha.up,plot=plot)
      type1=oc[1] #Expected type I error rate
      type2=oc[2] #Expected type II error rate
      risk=oc[3]  #Integrated risk
      
      typeIfr=integrand_Err(prior='vague',
                            n=x,n0=n0,alpha=alpha,priorParsInfo=info,
                            priorParsWeights=base,
                            lik=outcome,sampPrPars=c(th0,0),
                            th0=th0,w_adapt=TRUE,alpha1=alpha1,alpha.up=alpha.up)[1] #Frequentist type I error rate at th0
      
      typeIIfr=integrand_Err(prior='vague',
                             n=x,n0=n0,alpha=alpha,priorParsInfo=info,
                             priorParsWeights=base,
                             lik=outcome,sampPrPars=c(thR,0), 
                             th0=th0,w_adapt=TRUE,alpha1=alpha1,alpha.up=alpha.up)[2] #Frequentist type II error rate rate at thR
    }
    
    out=c('type1fr'= typeIfr,
          'type2fr'=typeIIfr,
          'ErrSum'= risk,
          'type1'=type1,
          'type2'=type2)
    out
  })
  colnames(Out_List_Wad)=n.list
  Out_List_Wad=list('Compromise-Ad'=Out_List_Wad)
  
  return(c(Out_List,Out_List_UW,Out_List_Wad,list('Sampling_prior_params'=true.pars)))
  
}

