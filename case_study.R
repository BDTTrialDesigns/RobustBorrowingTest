#Case study analysis of Section 5

designRunDataEx = function(dir.out,alpha.up=0.15,bound=TRUE)
{
  outcome='normal'

  #data
  data=0.40 
  sigma=0.44 
  n=1
  
  #hypothesis threshold
  th0=0

  #priors
  n0=1
  info = c(1,0.51, sqrt(0.016))
  vague = c(1,0,100)
  mix = cbind(c(0.55, 0.51, sqrt(0.016)), c(1-0.55, 0, sqrt(8.27)))
  power = c(1,0.51, sqrt(0.016)) 
  prior.list = list("vague"=vague,"CD-Adapt"=vague,"mix"=mix,"info"=info,"power"=power)
  alpha=0.025
  
  #type I error rate
  gamma=sigma^2/(n*info[3]^2)
  alpha1=1-pnorm(gamma*(th0-info[2])*sqrt(n)/sigma + qnorm(1-alpha) *sqrt(1+gamma))
  w.grid=seq(0,1,0.05)
  alpha.grid=(1-w.grid)*alpha+w.grid*alpha1
  
  #analysis
  lp=sapply(1:length(prior.list),function(y) {
    
    w_adapt=FALSE
    prior=names(prior.list)[y]
    
   if (names(prior.list)[y]=='CD-Adapt')
  {
      prior='vague'
      w_adapt=TRUE
   }
  
  priorParsWeights=prior.list[[y]]
  
  
  if(grepl("mix",prior)){
    priorPars=priorParsWeights[2:3,]
    weights=priorParsWeights[1,] 
  } else {
    priorPars=c(priorParsWeights[2],priorParsWeights[3])
    weights=priorParsWeights[1]
  }
  
  a0=1
  if (prior=="power") {
    a0=power_par(dat=data,n=1,dataPar=sigma,priorPars=priorPars,lik=outcome,bound=bound)
  }

  out=testDecision(dat=data,prior=prior,alpha=alpha,n=1,dataPar=sigma,
                         priorPars=priorPars,weights=weights,
                         lik=outcome,a0=a0,th0=th0,w_adapt=w_adapt,
                         alpha1=alpha1,alpha.up=alpha.up,priorParsInfo=info) 
  
  alpha.v=ifelse(w_adapt,(1-out[3])*alpha+out[3]*alpha1,alpha)
  c(out,alpha.v)
})
  colnames(lp)=c("FD","CD-Adapt","RMD-Unit (w=0.55)","BD","EBPD")
  lp[3,1]=0
  rownames(lp)=c('Test decision',expression(Pr(H0|bar(y))),"Adaptive (post.) weight",'Threshold')
  lp=t(lp)
  
  #Table 1
  print(xtable(lp,digits=3))
  
  #Adaptive weight plot & OCs
  pdf(file = paste0(dir.out,"CaseStudy_Weight",alpha.up,".pdf"), width=5,height=5, pointsize=10)
  
  oc=mclapply(1, function(id) designF(index=id,outcome='normal',prior.pars=c(0.51, sqrt(0.016)),
                                              designGrid=matrix(c(0.51,sqrt(0.016)),1,2),th0=0,
                                              thR=0.15,sigma=sigma,cores=1,
                                              plot=TRUE,n.list=1,alpha.up=alpha.up,S=10,
                                              case.study=TRUE,bound=bound),mc.cores = 1)
  dev.off()

  #OCs plot
  scale_color_colorblind7 = function(.ColorList = 2L:9L){
    scale_color_discrete(type = colorblind_pal()(8)[.ColorList])
  }
  
  n.list=1  
  w.list=seq(0,1,by=0.05)
  
  Out_List=lapply(oc[[1]][1:21],function(y) {names(y)=c("Comp","Mix","TI RBD","Mix Vague");y})
  Out_List_UW=oc[[1]][22:24]
  Out_List_Bound=oc[[1]][25]
  
  data.oc= lapply(c("Comp","Mix","TI RBD","Mix Vague"), function(y) 
  {
    sapply(Out_List,function(x)
    {
      rls=(x[[y]]['ErrSum','1']-Out_List_UW$Info['ErrSum','1'])/(Out_List_UW$Vague['ErrSum','1']-Out_List_UW$Info['ErrSum','1'])
      typeIfr=x[[y]]['type1fr','1']
      type2=x[[y]]['type2','1']
      risk=x[[y]]['ErrSum','1']
      out=c('RSL'=rls,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
    })
  })
  
  w.grid=seq(0,1,0.05)
  rc="FD"; rmc="MD Vague"; rmu='RMD Unit'
  outComp=data.frame('RSL'=c((do.call(cbind,data.oc)[1,]),rep(1,length(w.grid)),rep(0,length(w.grid)),
                             rep((Out_List_Bound[[1]]['ErrSum','1']-Out_List_UW$Info['ErrSum','1'])/(Out_List_UW$Vague['ErrSum','1']-Out_List_UW$Info['ErrSum','1']),length(w.grid)),
                             rep((Out_List_UW$Power['ErrSum','1']-Out_List_UW$Info['ErrSum','1'])/(Out_List_UW$Vague['ErrSum','1']-Out_List_UW$Info['ErrSum','1']),length(w.grid))),
                     'Freq type I'=c((do.call(cbind,data.oc)[2,]),rep(Out_List_UW$Vague['type1fr','1'],each=length(w.grid)),
                                     rep(Out_List_UW$Info['type1fr','1'],each=length(w.grid)),
                                     rep(Out_List_Bound[[1]]['type1fr','1'],each=length(w.grid)),
                                     rep(Out_List_UW$Power['type1fr','1'],each=length(w.grid))),
                     'Avg.Power'=c((do.call(cbind,data.oc)[3,]),
                                   rep(1-Out_List_UW$Vague['type2','1'],each=length(w.grid)),
                                   rep(1-Out_List_UW$Info['type2','1'],each=length(w.grid)),
                                   rep(1-Out_List_Bound[[1]]['type2','1'],each=length(w.grid)),
                                   rep(1-Out_List_UW$Power['type2','1'],each=length(w.grid))),
                     'w'=rep(w.grid,8),
                     'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','CD-Adapt','EBPD'),each=length(w.grid)))
  
  outComp$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  theme=ggplot()+theme_light() +
    scale_color_colorblind7() +
    xlab('')
  
  
  g1 <- theme+ geom_line(aes(w, RSL,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("RSL")+ labs(color = "Decision")  + scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1,1.2,1.4), limits= c(0,1.4))
  
  g2 <- theme+ geom_line(aes(w, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("Type I error rate")+ labs(color = "Decision")
  
  g3 <- theme+ geom_line(aes(w, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("Exp. Power")+ labs(color = "Decision") +
    xlab('w')
  
  a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.text=element_text(size=10), legend.title=element_blank()) ,
                      g2,g3,ncol = 1,nrow = 3,common.legend = T,align="v")
  
  pdf(file = paste0(dir.out,"CaseStudy_OCs",alpha.up,".pdf"), width=7,height=9, pointsize=10)
  print(a)
  dev.off()
}
