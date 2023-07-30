plots.design = function(outcome, dir.out) #outcome can be 'binomial' or 'normal' 
{
  #color palette
  scale_color_colorblind7 = function(.ColorList = 2L:9L){
    scale_color_discrete(type = colorblind_pal()(8)[.ColorList])
  }
  
  if(outcome=="binomial")
  {rc="FD"; rmc="MD PM"; rmu='RMD Unif'}
  if(outcome=="normal")
  {rc="FD"; rmc="MD Vague"; rmu='RMD Unit'}
  
  # Fixed sample size -------------------------------------------------------
  
  load(paste0(dir.out,outcome,"_bound_R1.RData"))
  
  n.list=c(1:250)  
  w.list=seq(0,1,by=0.01)
  
  Out_List=lapply(oc[1:length(w.list)],function(y) {names(y)=c("Comp","Mix","TI RBD","Mix Vague");y})
  Out_List_UW=oc[102:104]
  Out_List_Bound=oc[105]
  
  #Sample size = 100
  
  data.oc= lapply(c("Comp","Mix","TI RBD","Mix Vague"), function(y) 
  {
    sapply(Out_List,function(x)
    {
      rls=(x[[y]]['ErrSum','100']-Out_List_UW$Info['ErrSum','100'])/(Out_List_UW$Vague['ErrSum','100']-Out_List_UW$Info['ErrSum','100'])
      typeIfr=x[[y]]['type1fr','100']
      type2=x[[y]]['type2','100']
      risk=x[[y]]['ErrSum','100']
      out=c('RSL'=rls,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
    })
  })
  
  
  outComp=data.frame('RSL'=c((do.call(cbind,data.oc)[1,]),rep(1,length(w.list)),rep(0,length(w.list)),
                             rep((Out_List_Bound[[1]]['ErrSum','100']-Out_List_UW$Info['ErrSum','100'])/(Out_List_UW$Vague['ErrSum','100']-Out_List_UW$Info['ErrSum','100']),length(w.list)),
                             rep((Out_List_UW$Power['ErrSum','100']-Out_List_UW$Info['ErrSum','100'])/(Out_List_UW$Vague['ErrSum','100']-Out_List_UW$Info['ErrSum','100']),length(w.list))),
                     'Freq type I'=c((do.call(cbind,data.oc)[2,]),rep(Out_List_UW$Vague['type1fr','100'],each=length(w.list)),
                                     rep(Out_List_UW$Info['type1fr','100'],each=length(w.list)),
                                     rep(Out_List_Bound[[1]]['type1fr','100'],each=length(w.list)),
                                     rep(Out_List_UW$Power['type1fr','100'],each=length(w.list))),
                     'Avg.Power'=c((do.call(cbind,data.oc)[3,]),
                                   rep(1-Out_List_UW$Vague['type2','100'],each=length(w.list)),
                                   rep(1-Out_List_UW$Info['type2','100'],each=length(w.list)),
                                   rep(1-Out_List_Bound[[1]]['type2','100'],each=length(w.list)),
                                   rep(1-Out_List_UW$Power['type2','100'],each=length(w.list))),
                     'w'=rep(w.list,8),
                     'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','CD-Adapt','EBPD'),each=length(w.list)))
  
  outComp$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  theme=ggplot()+theme_light() +
    scale_color_colorblind7() +
    xlab('')
  
  
  g1 <- theme+ geom_line(aes(w, RSL,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("n=100")+ labs(color = "Decision")  + scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1,1.2,1.4), limits= c(0,1.4))
  
  g2 <- theme+ geom_line(aes(w, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("n=100")+ labs(color = "Decision")
  
  g3 <- theme+ geom_line(aes(w, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("n=100")+ labs(color = "Decision") +
    xlab('w')
  
  
  #Sample size = 20
  
  data.oc20=
    lapply(c("Comp","Mix","TI RBD","Mix Vague"), function(y) 
    {
      sapply(Out_List,function(x)
      {
        rls=(x[[y]]['ErrSum','20']-Out_List_UW$Info['ErrSum','20'])/(Out_List_UW$Vague['ErrSum','20']-Out_List_UW$Info['ErrSum','20'])
        typeIfr=x[[y]]['type1fr','20']
        type2=x[[y]]['type2','20']
        out=c('RSL'=rls,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
      })
      
    })
  
  
  outComp20=data.frame('RSL'=c((do.call(cbind,data.oc20)[1,]),rep(1,length(w.list)),rep(0,length(w.list)),
                               rep((Out_List_Bound[[1]]['ErrSum','20']-Out_List_UW$Info['ErrSum','20'])/(Out_List_UW$Vague['ErrSum','20']-Out_List_UW$Info['ErrSum','20']),length(w.list)),
                               rep((Out_List_UW$Power['ErrSum','20']-Out_List_UW$Info['ErrSum','20'])/(Out_List_UW$Vague['ErrSum','20']-Out_List_UW$Info['ErrSum','20']),length(w.list))),
                       'Freq type I'=c((do.call(cbind,data.oc20)[2,]),rep(Out_List_UW$Vague['type1fr','20'],each=length(w.list)),
                                       rep(Out_List_UW$Info['type1fr','20'],each=length(w.list)),
                                       rep(Out_List_Bound[[1]]['type1fr','20'],each=length(w.list)),
                                       rep(Out_List_UW$Power['type1fr','20'],each=length(w.list))),
                       'Avg.Power'=c((do.call(cbind,data.oc20)[3,]),
                                     rep(1-Out_List_UW$Vague['type2','20'],each=length(w.list)),
                                     rep(1-Out_List_UW$Info['type2','20'],each=length(w.list)),
                                     rep(1-Out_List_Bound[[1]]['type2','20'],each=length(w.list)),
                                     rep(1-Out_List_UW$Power['type2','20'],each=length(w.list))),
                       'w'=rep(w.list,8),
                       'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','CD-Adapt','EBPD'),each=length(w.list)))
  outComp20$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  
  g4 <- theme+ geom_line(aes(w, RSL,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp20) + ylab("RSL")+
    ggtitle("n=20")+ labs(color = "Decision") + scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1,1.2,1.4), limits= c(0,1.4))
  
  g5 <- theme+ geom_line(aes(w, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp20) + ylab("Type I error rate")+
    ggtitle("n=20")+ labs(color = "Decision")
  
  g6 <- theme+ geom_line(aes(w, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp20) + ylab("Exp. power")+
    ggtitle("n=20")+ labs(color = "Decision") +
    xlab('w')
  
  
  a=ggpubr::ggarrange(g4+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.text=element_text(size=10), legend.title=element_blank()) ,
                      g1,g5,g2,g6,g3,ncol = 2,nrow = 3,common.legend = T,align="v")
  
  pdf(file = paste0(dir.out,outcome,".pdf"), width=7,height=7, pointsize=10)
  grid.arrange(a,nrow = 1, ncol=1)
  graphics.off()
  
  
  
  # Sample size optimization ------------------------------------------------
  
  #Robust priors
  data.oc_ss=
    lapply(c("Comp","Mix","TI RBD","Mix Vague"), function(y) 
    {
      sapply(Out_List,function(x)
      {
        type2=x[[y]]['type2',]
        n=max(ifelse(max(type2)<0.2,0,min(which(type2==type2[type2>0.2][length(type2[type2>0.2])])+1,length(n.list))),1)
        type2=x[[y]]['type2',n]
        typeIfr=x[[y]]['type1fr',n]
        n.out=n.list[n]
        out=c('n'=n.out,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
      })
    })
  
  #Vague, informative, and power priors
  data.oc_ssUW=lapply(Out_List_UW, function(x){
    type2=x['type2',]
    n=max(ifelse(max(type2)<0.2,0,min(which(type2==type2[type2>0.2][length(type2[type2>0.2])])+1,length(n.list))),1)
    type2=x['type2',n]
    typeIfr=x['type1fr',n]
    n.out=n.list[n]
    out=c('n'=n.out,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
    out})
  
  #adapt prior
  data.oc_ssBound=lapply(Out_List_Bound, function(x){
    type2=x['type2',]
    n=max(ifelse(max(type2)<0.2,0,min(which(type2==type2[type2>0.2][length(type2[type2>0.2])])+1,length(n.list))),1)
    type2=x['type2',n]
    typeIfr=x['type1fr',n]
    n.out=n.list[n]
    out=c('n'=n.out,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
    out})
  
  
  
  ##Fixed sampling prior - varying weight
  
  outComp=data.frame('n'=c((do.call(cbind,data.oc_ss)[1,]),rep(data.oc_ssUW$Vague['n'],each=length(w.list)),rep(data.oc_ssUW$Info['n'],each=length(w.list)),
                           rep(data.oc_ssUW$Power['n'],each=length(w.list)),rep(data.oc_ssBound[[1]]['n'],each=length(w.list))),
                     'Freq type I'=c((do.call(cbind,data.oc_ss)[2,]),rep(data.oc_ssUW$Vague['Type I cond.'],each=length(w.list)),
                                     rep(data.oc_ssUW$Info['Type I cond.'],each=length(w.list)),rep(data.oc_ssUW$Power['Type I cond.'],each=length(w.list)),
                                     rep(data.oc_ssBound[[1]]['Type I cond.'],each=length(w.list))),
                     'Avg.Power'=c((do.call(cbind,data.oc_ss)[3,]),rep(data.oc_ssUW$Vague['Exp. power'],each=length(w.list)),
                                   rep(data.oc_ssUW$Info['Exp. power'],each=length(w.list)),rep(data.oc_ssUW$Power['Exp. power'],each=length(w.list)),
                                   rep(data.oc_ssBound[[1]]['Exp. power'],each=length(w.list))),
                     'w'=rep(w.list,8),
                     'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','EBPD','CD-Adapt'),each=length(w.list)))
  outComp$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  theme=ggplot()+theme_light() +
    scale_color_colorblind7()+
    xlab('') 
  
  if (outcome=="normal") {a1=0; a2=250; b1=-0.00001; b2=0.85}
  if (outcome=="binomial") {a1=0; a2=250; b1=-0.00001; b2=1}
  
  g1 <- theme+ geom_line(aes(w, n,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("Min. sample size")+
    labs(color = "Decision") + ylim(c(a1,a2)) #ylim(c(80,220))
  
  g2 <- theme+ geom_line(aes(w, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("Type I error rate")+
    labs(color = "Decision")
  
  g3 <- theme+ geom_line(aes(w, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("Exp. power")+
    labs(color = "Decision") + ylim(c(b1,b2))+
    xlab('w') 
  
  
  ##Fixed weight - varying sampling prior
  
  load(paste0(dir.out,outcome,"_sens_bound_R1.RData"))
  
  n.list=c(1:250)  
  if (outcome=='binomial')
  {
    designGrid=matrix(c(seq(0,40,by=2),(40-seq(0,40,by=2))),21,2)
    mean_vec=designGrid[,1]/rowSums(designGrid)
    mus=0.5
  }
  if(outcome=='normal')
  {
    designGrid=matrix(c(seq(-0.7,1.3,by=0.1),rep(1/50,21)),21,2)
    mean_vec=designGrid[,1]
    mus=0.25
  }
  
  Out_List=lapply(oc[1:length(mean_vec)],function(y) 
  {names(y$'0.5')=c("Comp","Mix","TI RBD","Mix Vague")
  out=list(y$'0.5'$"Comp",y$'0.5'$"Mix",y$'0.5'$"TI RBD",y$'0.5'$"Mix Vague",y$'Vague',y$'Info',y$'Power',y$'Compromise-Ad')
  names(out)=c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt')
  out})
  
  
  data.oc_ss=
    lapply(c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt'), function(y) 
    {
      sapply(Out_List,function(x)
      {
        type2=x[[y]]['type2',]
        n=max(ifelse(max(type2)<0.2,0,min(max(which(type2==type2[type2>0.2][length(type2[type2>0.2])]))+1,length(n.list))),1)
        type2=x[[y]]['type2',n]
        typeIfr=x[[y]]['type1fr',n]
        n.out=n.list[n]
        out=c('n'=n.out,'Type I cond.'=typeIfr, 'Exp. power'=1-type2)
      })
      
    })
  names(data.oc_ss)=c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt')
  
  
  outComp=data.frame('n'=c((do.call(cbind,data.oc_ss)[1,])),
                     'Freq type I'=c((do.call(cbind,data.oc_ss)[2,])),
                     'Avg.Power'=c((do.call(cbind,data.oc_ss)[3,])),
                     'theta'=rep(mean_vec,8),
                     'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','EBPD','CD-Adapt'),each=length(mean_vec)))
  outComp$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  theme=ggplot()+theme_light() +
    scale_color_colorblind7()+
    xlab('') + geom_vline(xintercept = mus, linetype="dashed", 
                          color = "gray", size=0.5)
  
  
  g4 <- theme+ geom_line(aes(theta, n,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    labs(color = "Decision") + ylab("")
  
  g5 <- theme+ geom_line(aes(theta, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    labs(color = "Decision")
  
  g6 <- theme+ geom_line(aes(theta, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    labs(color = "Decision")+
    xlab('Sampling prior mean')
  
  
  
  a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.text=element_text(size=10),legend.title=element_blank()),
                      g4,g2,g5,g3,g6,ncol = 2,nrow = 3,common.legend = T,align="v")
  
  
  pdf(file = paste0(dir.out,outcome,"_SampleSize.pdf"), width=7,height=7, pointsize=10)
  grid.arrange(a,nrow = 1, ncol=1)
  graphics.off()
  
  
  # Sensitivity analyses ----------------------------------------------------
  
  
  ##Fixed weight - varying sampling prior 
  
  load(paste0(dir.out,outcome,"_sens_bound_R1.RData"))
  
  n.list=c(1:250)  
  if (outcome=='binomial')
  {
    designGrid=matrix(c(seq(0,40,by=2),(40-seq(0,40,by=2))),21,2)
    mean_vec=designGrid[,1]/rowSums(designGrid)
  }
  if(outcome=='normal')
  {
    designGrid=matrix(c(seq(-0.7,1.3,by=0.1),rep(1/50,21)),21,2)
    mean_vec=designGrid[,1]
  }
  
  Out_List=lapply(oc[1:length(mean_vec)],function(y) 
  {names(y$'0.5')=c("Comp","Mix","TI RBD","Mix Vague")
  out=list(y$'0.5'$"Comp",y$'0.5'$"Mix",y$'0.5'$"TI RBD",y$'0.5'$"Mix Vague",y$'Vague',y$'Info',y$'Power',y$'Compromise-Ad')
  names(out)=c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt')
  out})
  
  data.oc_ss=
    lapply(c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt'), function(y) 
    {
      sapply(Out_List,function(x)
      {
        n=20
        type2=x[[y]]['type2',n]
        typeIfr=x[[y]]['type1fr',n]
        n.out=n.list[n]
        risk=x[[y]]['ErrSum',n]
        out=c('n'=n.out,'Type I cond.'=typeIfr, 'Exp. power'=1-type2,'Int.risk'=risk)
      })
      
    })
  names(data.oc_ss)=c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt')
  
  
  outComp=data.frame('Int.Risk'=c((do.call(cbind,data.oc_ss)[4,])),
                     'Freq type I'=c((do.call(cbind,data.oc_ss)[2,])),
                     'Avg.Power'=c((do.call(cbind,data.oc_ss)[3,])),
                     'theta'=rep(mean_vec,8),
                     'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','EBPD','CD-Adapt'),each=length(mean_vec)))
  outComp$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  theme=ggplot()+theme_light() +
    scale_color_colorblind7()+
    xlab('')+ geom_vline(xintercept = mus, linetype="dashed", 
                         color = "gray", size=0.5)
  
  
  g1 <- theme+ geom_line(aes(theta, Int.Risk,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("Integrated risk")+
    ggtitle("n=20")+labs(color = "Decision")
  
  g2 <- theme+ geom_line(aes(theta, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("Type I error rate")+
    ggtitle("n=20")+labs(color = "Decision")
  
  g3 <- theme+ geom_line(aes(theta, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("Exp. power")+
    ggtitle("n=20")+labs(color = "Decision")+
    xlab('Sampling prior mean')
  
  
  data.oc_ss=
    lapply(c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt'), function(y) 
    {
      sapply(Out_List,function(x)
      {
        n=100
        type2=x[[y]]['type2',n]
        typeIfr=x[[y]]['type1fr',n]
        n.out=n.list[n]
        risk=x[[y]]['ErrSum',n]
        out=c('n'=n.out,'Type I cond.'=typeIfr, 'Exp. power'=1-type2,'Int.risk'=risk)
      })
      
    })
  names(data.oc_ss)=c("Comp","Mix","TI RBD","Mix Vague",'Vague','Info','EB Power','CD-Adapt')
  
  
  outComp=data.frame('Int.Risk'=c((do.call(cbind,data.oc_ss)[4,])),
                     'Freq type I'=c((do.call(cbind,data.oc_ss)[2,])),
                     'Avg.Power'=c((do.call(cbind,data.oc_ss)[3,])),
                     'theta'=rep(mean_vec,8),
                     'Decision'=rep(c('CD',rmu,'TI RBD',rmc,rc,'BD','EBPD','CD-Adapt'),each=length(mean_vec)))
  outComp$Decision=factor(outComp$Decision,levels=c('CD','CD-Adapt',rmu,'BD',rc,'EBPD',rmc,'TI RBD'))
  
  theme=ggplot()+theme_light() +
    scale_color_colorblind7()+
    xlab('')+ geom_vline(xintercept = mus, linetype="dashed", 
                         color = "gray", size=0.5)
  
  
  g4 <- theme+ geom_line(aes(theta, Int.Risk,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("n=100")+labs(color = "Decision") + ylab("")
  
  g5 <- theme+ geom_line(aes(theta, Freq.type.I,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("n=100")+labs(color = "Decision")
  
  g6 <- theme+ geom_line(aes(theta, Avg.Power,  color=Decision, linetype=Decision), 
                         size=0.8, data=outComp) + ylab("")+
    ggtitle("n=100")+labs(color = "Decision")+
    xlab('Sampling prior mean')
  
  
  a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.text=element_text(size=10),legend.title=element_blank()),
                      g4,g2,g5,g3,g6,ncol = 2,nrow = 3,common.legend = T,align="v")
  
  
  pdf(file = paste0(dir.out,outcome,"_Sensitivity.pdf"), width=7,height=7, pointsize=10)
  grid.arrange(a,nrow = 1, ncol=1)
  graphics.off()
}

# Frequentist risk plot ---------------------------------------------------

plot.fr = function(dir.out) #outcome can be 'binomial' or 'normal' 
{
  #data parameters
  n=100
  sigma=1
  n0=50
  
  #prior & grid
  prior.pars=c(-0.5,0.5,sqrt(sigma/n0))
  grid=round(seq(-0.7,1.3,by=0.01),2)
  
  #test error costs & thresholds
  c.alpha = 0.975
  c.beta = 0.025
  alpha=c.beta/(c.beta+c.alpha)
  
  th0=0
  gamma=sigma^2/(n*prior.pars[3]^2)
  alpha1=1-pnorm(gamma*(th0-prior.pars[2])*sqrt(n)/sigma + qnorm(1-alpha) *sqrt(1+gamma))
  alphaS=0.5*alpha+0.5*alpha1
  alpha_v=c(rep(alpha,3),alphaS)
  
  #priors
  info0 <-c(1,prior.pars[1],prior.pars[3])
  info1 <-c(1,prior.pars[2],prior.pars[3])
  base <-c(1,0,100)
  prior.list = list(base, info0, info1, base)
  names(prior.list)=c("vague","info","info","vague")
  
  #compute frequentist risk
  lp=sapply(1:length(prior.list),function(y) {
    
    fr.risk= mclapply(grid, function(i)
    {
      integrand_Err(prior=names(prior.list)[y],
                    n=n,n0=n0,alpha=alpha_v[y],
                    dataPar=sigma,priorParsWeights=prior.list[[y]],
                    lik='normal',sampPrPars=c(i,0),
                    th0=th0,S=6)[3]
    },mc.cores=1)
    
    return(do.call(c,fr.risk))
  })
  
  data.risk=data.frame('theta'=grid,'risk'=as.vector(lp),
                       'Decision'=rep(c('Frequentist','Bayes - Null','Bayes - Alt.', 'Compromise with Bayes - Alt.'),each=nrow(lp))
                       )
  
  theme=ggplot()+theme_light() +
    scale_color_viridis_d(begin = 0.1, end = 0.9)
  
  
  g1 <- theme+ geom_line(aes(theta, risk, group=Decision, color=Decision), 
                         size=0.5, data=subset(data.risk,theta<=0.00))+ 
    geom_line(aes(theta, risk, group=Decision, color=Decision), 
    size=0.5, data=subset(data.risk,theta>0.00)) +
    ylab(expression(paste("Frequentist risk ",R(theta)))) + xlab(expression(theta)) +
    geom_vline(xintercept = th0, linetype="dashed", 
               color = "gray", size=0.5) +
    geom_point(aes(x=theta, y=risk, color=Decision), 
               data=subset(data.risk,theta==0.00))

  pdf(file = paste0(dir.out,"risk_ex.pdf"), width=7,height=5, pointsize=10)
  grid.arrange(g1,nrow = 1, ncol=1)
  graphics.off()
  
}
