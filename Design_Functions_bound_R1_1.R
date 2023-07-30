# Power prior weight ------------------------------------------------------

#adapted from studyPrior package
power_par = function(dat,n,dataPar,priorPars,lik,bound=TRUE) {
  if(lik=="normal")
  {
    sD = dataPar^2/n
    d = priorPars[2]^2 / (pmax((dat - priorPars[1]) ^ 2, sD + priorPars[2]^2) - sD)
    if(bound) d = pmin(d,priorPars[2]^2/(dataPar^2/n))
  }
  if(lik=="binomial")
  {
    optfn=function(datC)
    {
      lik.d = function(d) VGAM::dbetabinom.ab(x=datC, size=n, d * (priorPars[1] -1) +1 ,  d *(priorPars[2]-1) +1) #uniform historical prior
      
      opd = BB::spg(par = .005,
                    fn = lik.d,
                    lower=0,
                    upper=1,
                    control=list(maximize=TRUE,
                                 trace=FALSE))
      
      if(opd$convergence!=0) print(opd)
      
      ds = opd$par
      
      return(ds)
    }
    d=sapply(dat,optfn)
    if(bound) d=pmin(d,n/(priorPars[1]+priorPars[2]))
  }
  return(d)
}


# Obtain test decision ----------------------------------------------------


testDecision= function(dat,prior,alpha,
                       n,n0=NULL,dataPar,
                       priorPars,weights,
                       lik,a0,th0,w_adapt=FALSE,
                       alpha1,alpha.up=NULL,priorParsInfo=NULL) 
{
  wT=1
  
  if (prior=="vague" | prior=="info" | prior=="sampling"){
    if(lik=="normal"){
      
      post.var = 1/(1/priorPars[2]^2+n/dataPar^2)
      post.mean= post.var*(priorPars[1]/(priorPars[2]^2)+dat/(dataPar^2/n))

      if(w_adapt){
        post.var.vague = 1/(1/priorPars[2]^2+n/dataPar^2)
        post.var.info = 1/(1/priorParsInfo[3]^2+n/dataPar^2)
        post.mean.info = post.var.info*(priorParsInfo[2]/(priorParsInfo[3]^2)+dat/(dataPar^2/n))
        wT=1-abs(pnorm(th0,dat,sqrt(post.var.info))-pnorm(th0,post.mean.info,sqrt(post.var.info)))
        alpha.w=pmin(rep(alpha.up,length(dat)),(1-wT)*alpha+wT*alpha1)
        wT=(alpha.w-rep(alpha,length(dat)))/(rep(alpha1,length(dat))-rep(alpha,length(dat)))
      } else {
        alpha.w=alpha
      }
      
      p.post=pnorm(th0,post.mean,sqrt(post.var))
      dec = (p.post<=alpha.w)
    }
    if(lik=="binomial") {
      post.a=priorPars[1]+dat
      post.b=priorPars[2]+(n-dat)
      
      if(w_adapt){
        post.a.info=priorParsInfo[2]+dat
        post.b.info=priorParsInfo[3]+(n-dat)
        wT=1-abs(pbeta(th0,post.a.info,post.b.info)-pbeta(th0,1+dat+dat*((n0/n)),1+n+n0-dat-dat*((n0/n))))
        alpha.w=pmin(rep(alpha.up,length(dat)),(1-wT)*alpha+wT*alpha1)
        wT=(alpha.w-rep(alpha,length(dat)))/(rep(alpha1,length(dat))-rep(alpha,length(dat)))
      }
      else{
        alpha.w=alpha
      }
      
      p.post=pbeta(th0,post.a,post.b)
      dec= (p.post<=alpha.w)
    }
  }
  if(grepl("mix",prior))  {#adapted from RBesT package
    if(lik=="normal"){
      post.var = 1/(1/priorPars[2,]^2+n/dataPar^2)
      post.mean = matrix(NA,length(dat),2)
      gamma=(dataPar^2/n)/((dataPar^2/n)+priorPars[2,]^2)
      post.mean[,1]=gamma[1]*priorPars[1,1]+(1-gamma[1])*dat
      post.mean[,2]=gamma[2]*priorPars[1,2] + (1-gamma[2])*dat
      dataParPred=matrix(NA,length(dat),2)
      dataParPred[,1] = sqrt(priorPars[2,1]^2 + dataPar^2/n)
      dataParPred[,2] = sqrt(priorPars[2,2]^2 + dataPar^2/n)
      
      margT= exp(log(weights[1]) + dnorm(dat, priorPars[1,1], dataParPred[,1], log=TRUE)) + exp(log(weights[2]) + dnorm(dat, priorPars[1,2], dataParPred[,2], log=TRUE))
      marg1= exp(log(weights[1]) + dnorm(dat, priorPars[1,1], dataParPred[,1], log=TRUE))
      post.weights=cbind(exp(log(marg1)-log(margT)),1-exp(log(marg1)-log(margT)))
      wT=post.weights[,1]
      p.post= (1 - post.weights[,2]) * pnorm(th0, post.mean[,1], rep(sqrt(post.var[1]),length(dat))) + post.weights[,2] * pnorm(th0, post.mean[,2], rep(sqrt(post.var[2]),length(dat)))
      dec=(p.post<=alpha)
    }
    if(lik=="binomial") {
      if(priorPars[2,2]!=0){
        post.a=matrix(NA,length(dat),2)
        post.a[,1]=priorPars[1,1]+dat
        post.a[,2]=priorPars[1,2]+dat
        post.b=matrix(NA,length(dat),2)
        post.b[,1]=priorPars[2,1]+(n-dat)
        post.b[,2]=priorPars[2,2]+(n-dat)
        
        marg1=log(weights[1]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,1], shape2=priorPars[2,1], log=TRUE)
        marg2=log(weights[2]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,2], shape2=priorPars[2,2], log=TRUE)
        
        margT= exp(marg1) + exp(marg2)
        post.weights=cbind(exp(marg1-log(margT)),1-exp(marg1-log(margT)))
        p.post= (1 - post.weights[,2]) * pbeta(th0, post.a[,1], post.b[,1]) + post.weights[,2] * pbeta(th0, post.a[,2], post.b[,2])
        wT=post.weights[,1]
        dec=(p.post<=alpha)
      }
      else {
        post.a=matrix(NA,length(dat),2)
        post.a[,1]=priorPars[1,1]+dat
        post.b=matrix(NA,length(dat),2)
        post.b[,1]=priorPars[2,1]+(n-dat)
        
        marg1=log(weights[1]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,1], shape2=priorPars[2,1], log=TRUE)
        marg2=log(weights[2]) + dbinom(x=dat, size=n, prob=th0, log=TRUE)
        
        margT= exp(marg1) + exp(marg2)
        post.weights=cbind(exp(marg1-log(margT)),1-exp(marg1-log(margT)))
        wT=post.weights[,1]
        p.post= (1 - post.weights[,2]) * pbeta(th0, post.a[,1], post.b[,1]) + post.weights[,2] * 1
        dec=(p.post<=alpha)
      }
    }
  }
  
  if (prior=="power"){
    if(lik=="normal") {
      prior.sd.EB=priorPars[2]/sqrt(a0)
      post.var=  1/(1/prior.sd.EB^2+n/dataPar^2)
      post.mean= post.var*(priorPars[1]/(prior.sd.EB^2)+dat/(dataPar^2/n))
      p.post=pnorm(th0,post.mean,sqrt(post.var))
      wT=a0
      dec=(p.post<=alpha)
    }
    if(lik=="binomial") {
      post.a =  a0 * (priorPars[1] -1) +1 + dat
      post.b =  a0 *(priorPars[2]-1) +1 +(n-dat)
      p.post=pbeta(th0,post.a,post.b)
      wT=a0
      dec= (p.post<=alpha)
    }
  }
  return(cbind(dec,p.post,wT))
}



# Integrated risk and error rates -----------------------------------------


integrand_Err = function(prior,n,n0=NULL,alpha,dataPar=NULL,
                         priorParsWeights,lik=c('binomial','normal'),sampPrPars,
                         th0,w_adapt=FALSE,alpha1=NULL,priorParsInfo=NULL,
                         alphas=0.025,alpha.up=NULL,plot=FALSE,
                         S=6,Fr=400,bound=NULL){
  if(grepl("mix",prior)){
    priorPars=priorParsWeights[2:3,]
    weights=priorParsWeights[1,] 
  } else {
    priorPars=c(priorParsWeights[2],priorParsWeights[3])
    weights=priorParsWeights[1]
  }
  
  #data outcomes grid
  if(lik=="normal"){
    #data integration grid
    x1=seq(sampPrPars[1]-S*sqrt(sampPrPars[2]^2+dataPar^2/n),sampPrPars[1]+S*sqrt(sampPrPars[2]^2+dataPar^2/n),by=sqrt(sampPrPars[2]^2+dataPar^2/n)/Fr)
    wi=dnorm(x1,sampPrPars[1],sqrt(sampPrPars[2]^2+dataPar^2/n))
    dy=t(wi/sum(wi))
    
    post.var.s = 1/(1/sampPrPars[2]^2+n/dataPar^2)
    if(sampPrPars[2]==0) {post.mean.s=sampPrPars[1]} else {post.mean.s=post.var.s*(sampPrPars[1]/(sampPrPars[2]^2)+x1/(dataPar^2/n))}
    ppos= pnorm(th0,post.mean.s,sqrt(post.var.s))
    pprs= pnorm(th0,sampPrPars[1],sampPrPars[2])
    
  }
  
  if(lik=="binomial"){
    x1=0:n
    if(sampPrPars[2]==0){
      ths=sampPrPars[1]
      dy <-dbinom(x1, size = n, ths)
      
      ppos= ths<=th0
      pprs= ths<=th0      
    } else{
      
      dy <- dbetabinom.ab(x1, 
                          size = n, 
                          shape1 = sampPrPars[1], 
                          shape2 = sampPrPars[2])  
      
      post.a.s=sampPrPars[1]+x1
      post.b.s=sampPrPars[2]+(n-x1)
      ppos= pbeta(th0,post.a.s,post.b.s)
      pprs= pbeta(th0,sampPrPars[1],sampPrPars[2])}
  }
  
  #power parameter
  a0=rep(1,length(x1))
  if (prior=="power") {
    if(lik=="binomial") {
      a0=power_par(dat=x1,n=n,dataPar=dataPar,priorPars=priorPars,lik=lik,bound=bound)
    }
    if(lik=="normal") {
      a0=power_par(dat=x1,n=n,dataPar=dataPar,priorPars=priorPars,lik=lik,bound=bound)
    }
  }
  

  #test decision & error rates
  if(lik=="binomial") {
    dec.w=testDecision(dat=x1,prior=prior,alpha=alpha,n=n,n0=n0,dataPar=dataPar,
                     priorPars=priorPars,weights=weights,lik=lik,a0=a0,th0=th0,
                     w_adapt=w_adapt,alpha1=alpha1,priorParsInfo=priorParsInfo,alpha.up=alpha.up)
    dec=dec.w[,1]
    w.out=dec.w[,3]
    
    t1e=dy%*%((ppos*dec)/pprs) 
    t2e=dy%*%(((1-ppos)*(1-dec))/(1-pprs))
    ir=dy%*%(ppos*dec+(alphas/(1-alphas))*(1-ppos)*(1-dec)) 
  }
  if(lik=="normal") {
    dec.w= testDecision(dat=x1,prior=prior,alpha=alpha,n=n,dataPar=dataPar,
                      priorPars=priorPars,weights=weights,lik=lik,a0=a0,th0=th0,
                      w_adapt=w_adapt,alpha1=alpha1,priorParsInfo=priorParsInfo,alpha.up=alpha.up)
    dec=dec.w[,1]
    w.out=dec.w[,3]
    

    t1e=dy%*%((ppos*dec)/pprs) #(expected) type I error
    t2e=dy%*%(((1-ppos)*(1-dec))/(1-pprs))  #(expected) type II error
    ir=dy%*%(ppos*dec+(alphas/(1-alphas))*(1-ppos)*(1-dec)) #integrated risk
  }
  
  #plot of CD-Adapt weight
  if(w_adapt & plot)
  {par(mar = c(5, 4, 4, 4) + 0.25)
    plot(x1[order(x1)],w.out[order(x1)],xlab=expression(bar(y)),ylab=expression(hat(w)),main='Compromise - Adaptive weight',type='l',lwd=2)
    axis(4)
    mtext(expression(f(bar(y)~";"~theta[0])), side = 4, line = 3, col = 2)
    plot(function(x) dnorm(x,th0,sqrt(dataPar/n)),add=TRUE,type='l',col=2,lwd=2,xlim=c(min(x1),max(x1)))
    }

  return(c(t1e,t2e,ir))
}