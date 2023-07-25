rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

###### RUN PARAMS
Nsim = 8000; Nspag = 500 ; textsize = 15
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
starthist = 1500
endhist = 1815
starthistartif = 1816
endhistartif = 1969

#### Run or plot only ?
RUN = F

#### PLOTS EN
dir.plots = "C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/PlotsEN/"

###### DATA LOADING
###Histrhone data 
Call = read.csv2(paste0(dir.data,"CX_All.csv"))
### 1816-2020 AMAX discharges with uncertainty and spags
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
Q = read.table(
  "C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]

###### MODELS AND PARAMS
#### 4 samples
samples =  c("C3","C4","Artif1","Artif2")
Qconts =   list(Q, Q, Q[which(Q$an>endhistartif),], Q[which(Q$an>endhistartif),])
Spagsss =  list(Spags[, (1:Nspag) ], Spags[, (1:Nspag) ],
                Spags[which(Q$an>endhistartif), (1:Nspag) ],
                Spags[which(Q$an>endhistartif), (1:Nspag) ])
Occs =     list(Call[which(Call$An>=starthist & Call$An<=endhist), 1],
                Call[which(Call$An>=starthist & Call$An<=endhist & Call$Cat == "C4"), 1],
                Q$an[which(Q$mp > 7000 & Q$an <= endhistartif)],
                Q$an[which(Q$mp > 9000 & Q$an <= endhistartif)]  )
Nhists = list(endhist-starthist+1,endhist-starthist+1,
             endhistartif-starthistartif+1,endhistartif-starthistartif+1)
#### Threshold priors
Thres = c(7000,9000,7000,9000)
sdThres = 2000
#### Nban priors (when Nban is draw in Unif dist.) : end of histo. period - date of first OT flood 
NbanS = c(endhist-Call$An[which(Call$An > 1500)][1]+1,
          endhist-Call$An[which(Call$An > 1500 & Call$Cat == "C4")][1]+1,
          endhistartif-Q$an[which(Q$mp > 7000)][1]+1,
          endhistartif-Q$an[which(Q$mp > 9000)][1]+1)

StartH = c(starthist,starthist,starthistartif,starthistartif)

MaxUnif_Nban = 500
#### 4 models
models =  c("A","B","C","D")
shortm = c(": u0",": uS",": ut*",": uS&t*")
distsTs = c("FIX","Gaussian","FIX","Gaussian")
distsNb = c("FIX","FIX","Uniform","Uniform")
#### Baselines
Baselines =     c("1816-2020","1970-2020")
Qbaselines =    list(Q, Q[which(Q$an>endhistartif),])
SpagsBaseline = list(Spags[, (1:Nspag) ], Spags[which(Q$an>endhistartif), (1:Nspag) ])

#### create dirs
# for(s in 1:length(samples)){
  # ds=paste0(dir.res,samples[s]); dir.create(ds, showWarnings = F)
  # for(m in 1:length(models)){
    # dir.create(paste0(ds,"/model_",models[m]), showWarnings = F)}
    # dir.create(paste0(dir.plots,samples[s]),showWarnings = F)}
# for(b in 1:length(Baselines)){dir.create(paste0(dir.res,"/Baseline_",Baselines[b]),showWarnings = F)}

# verif
# for(s in 1:length(samples)){
#   for(m in 1:length(models)){
#       print(samples[s])
#       print(models[m])
#       print( length(Qconts[[s]][,1]))
#       print( dim(Spagsss[[s]]))
#       print( length(data.frame(an=Occs[[s]])[,1]))
#       print(list(distsTs[m],c(Thres[s],sdThres),Thres[s]))
#       print(list(distsNb[m],c(NbanS[s],MaxUnif_Nban),NbanS[s]))
#       
#   }
# }

# 
#### RUN MODELS
if(RUN == T){
  #### Run GEV_Binoms
  for(s in 1:length(samples)){
    for(m in 1:length(models)){
      dir.case = paste0(dir.res,samples[s],"/model_",models[m])
      #true number of years for FIX distribution vs first flood number for Unif
      if(distsNb[m] == "FIX"){ aNb = Nhists[[s]] } else { aNb = NbanS[s] }
      # Fun run
      GEV_Binom(Qcont = Qconts[[s]],
                SpagCont = Spagsss[[s]],
                Occ = data.frame(an=Occs[[s]]),
                Ts = list(distsTs[m],c(Thres[s],sdThres),Thres[s]),#dist,prior,init
                Nban = list(distsNb[m],c(aNb,MaxUnif_Nban),aNb),#dist,prior,init
                prob = prob,
                Nsim = Nsim,
                Nhist = Nhists[[s]],
                dir.res = dir.case)
    }
  }
### RUN Baselines
  for(b in 1:length(Baselines)){
    dir.case = paste0(dir.res,"/Baseline_",Baselines[b])
    # Fun run
    GEV_Spag(Q = Qbaselines[[b]],
             Spags = SpagsBaseline[[b]],
             prob = prob,Nsim = Nsim,dir.res = dir.case)
  }
}

rm(Spagsss); rm(SpagsBaseline)

ind100 = length(Pr) - 18
ind1000 = length(Pr) - 9

#### Results reading init
Quant100 = list(); Quant1000 = list(); ParamAll = list(); Mp_par = list(); QuantsAll = list()
palmod = brewer.pal(length(models),"RdYlBu")[c(2,1,3,4)] 
palbase = c("#4dac26","#b8e186")
colprior = "lightgrey"
  
for(s in 1:length(samples)){
  # s=4
  GGQuants = list()
  for(m in 1:length(models)){
    dir.case = paste0(dir.res,samples[s],"/model_",models[m])
    ## All quants
    QuantsAll[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[,(2:4)]
    QuantsAll[[m]]$model = paste0(models[m],shortm[m]); QuantsAll[[m]]$sample = samples[s]
    ##Q100 & Q1000
    Quant100[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[ind100,(2:4)]
    Quant100[[m]]$model = paste0(models[m],shortm[m]); Quant100[[m]]$sample = samples[s]
    Quant1000[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[ind1000,(2:4)]
    Quant1000[[m]]$model = paste0(models[m],shortm[m]); Quant1000[[m]]$sample = samples[s]
    ## posterior params
    ParamAll[[m]] = read.table(paste0(dir.case,"/Params.txt"), header = T)
    ParamAll[[m]]$model = paste0(models[m],shortm[m]); ParamAll[[m]]$sample = samples[s]
    ## maxpost params
    Mp_par[[m]] = read.table(paste0(dir.case,"/Maxpost_Par.txt"), header = T)
    Mp_par[[m]]$model = paste0(models[m],shortm[m]); Mp_par[[m]]$sample = samples[s]
  }
  for(b in 1:length(Baselines)){
    m = length(models) + b
    dir.case = paste0(dir.res,"Baseline_",Baselines[b])
    ## All quants
    QuantsAll[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[,(2:4)]
    QuantsAll[[m]]$model = paste0("GEV",Baselines[b]); QuantsAll[[m]]$sample = Baselines[b]
    ##Q100 & Q1000
    Quant100[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[ind100,(2:4)]
    Quant100[[m]]$model = paste0("GEV ",Baselines[b]); Quant100[[m]]$sample = Baselines[b]
    Quant1000[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[ind1000,(2:4)]
    Quant1000[[m]]$model = paste0("GEV ",Baselines[b]); Quant1000[[m]]$sample = Baselines[b]
    ## posterior params
    ParamAll[[m]] = read.table(paste0(dir.case,"/Params.txt"), header = T)
    ParamAll[[m]]$model = paste0("GEV ",Baselines[b]); ParamAll[[m]]$sample = Baselines[b]
    ## maxpost params
    Mp_par[[m]] = read.table(paste0(dir.case,"/Maxpost_Par.txt"), header = T)
    Mp_par[[m]]$model = paste0("GEV ",Baselines[b]); Mp_par[[m]]$sample = Baselines[b]
    if(s<3){break}
  }
  
  ### RESULTS TIDYING
  Params = melt.list(ParamAll) ; Params$model = as.factor(Params$model)
  Params$model = fct_relevel(Params$model, 
                             as.vector(unique(Params$model)[c(5:length(Mp_par),1,2,3,4)]) )
  Quants100 = cast(melt.list(Quant100) );  Quants1000 = cast(melt.list(Quant1000) ) 
  ParamsMp = melt.list(Mp_par) ; ParamsMp$model = as.factor(ParamsMp$model)
  ParamsMp$model = fct_relevel(ParamsMp$model, 
                             as.vector(unique(ParamsMp$model)[c(5:length(Mp_par),1,2,3,4)]) )
  Quants100$model = as.factor(Quants100$model);Quants1000$model = as.factor(Quants1000$model)
  Quants100$model =fct_relevel(Quants100$model,
                               as.vector(Quants100$model[c(5:length(Mp_par),1,2,3,4)]))
  Quants1000$model=fct_relevel(Quants1000$model,
                               as.vector(Quants1000$model[c(5:length(Mp_par),1,2,3,4)]))

  
  #### DF RESUME RESULTS
  ResAll = data.frame(model = c(models,paste0("GEV",Baselines[1:(-4+length(Quants100$model))])),
                      sd.S=NA, mp.S=NA, sd.t = NA, mp.t = NA,
                      sd.xi = NA, mp.xi = NA,  sd.Q100 = NA, mp.Q100 = NA, 
                      sd.Q1000 = NA, mp.Q1000 = NA) 
  
  endPhist = Qconts[[s]]$an[1]
  for(m in 1:length(models)){
    ### S
    ResAll[m,]$sd.S = round(sd(Params$value[which(Params$model == paste0(models[m],shortm[m]) & 
                                         Params$variable == "Thres") ] ) )
    ResAll[m,]$mp.S = round(ParamsMp$value[which(ParamsMp$model == paste0(models[m],shortm[m]) & 
                                             ParamsMp$variable == "Seuil") ] )
    ### t*
    ResAll[m,]$sd.t = round(sd(Params$value[which(Params$model == paste0(models[m],shortm[m]) & 
                                              Params$variable == "Nban") ] ) )
    ResAll[m,]$mp.t = round(endPhist - ParamsMp$value[which(ParamsMp$model == paste0(models[m],shortm[m]) & 
                                             ParamsMp$variable == "NBan") ] )
    ### xi
    ResAll[m,]$sd.xi = round(sd(Params$value[which(Params$model == paste0(models[m],shortm[m]) & 
                                              Params$variable == "Shape") ] ),3)
    ResAll[m,]$mp.xi = round(ParamsMp$value[which(ParamsMp$model == paste0(models[m],shortm[m]) & 
                                             ParamsMp$variable == "Form") ],3)
    ### Q100
    ResAll[m,]$sd.Q100 = round( ((Quants100$Q_9 - Quants100$Mp)/2)[m] )
    ResAll[m,]$mp.Q100 = round( Quants100$Mp[m] ) 
    #Q1000
    ResAll[m,]$sd.Q1000 = round( ((Quants1000$Q_9 - Quants1000$Mp)/2)[m] )
    ResAll[m,]$mp.Q1000 = round( Quants1000$Mp[m] )
  } 
  ### BASELINES
  for(b in 1:(-4+length(Quants100$model))){
    ### xi
    ResAll[(4+b),]$sd.xi = round(sd(Params$value[which(
                                          Params$model == paste0("GEV ", Baselines[b]) & 
                                          Params$variable == "Shape") ] ),3)
    ResAll[(4+b),]$mp.xi = round(ParamsMp$value[which(
                                          ParamsMp$model == paste0("GEV ", Baselines[b]) & 
                                          ParamsMp$variable == "Form") ],3)
    ### Q100
    ResAll[(4+b),]$sd.Q100 = round( ((Quants100$Q_9 - Quants100$Mp)/2)[(4+b)] )
    ResAll[(4+b),]$mp.Q100 = round( Quants100$Mp[(4+b)] ) 
    #Q1000
    ResAll[(4+b),]$sd.Q1000 = round( ((Quants1000$Q_9 - Quants1000$Mp)/2)[(4+b)] )
    ResAll[(4+b),]$mp.Q1000 = round( Quants1000$Mp[(4+b)] )
  }
  ### write resume results
  write.csv2(ResAll, paste0(dir.res, samples[s],"/ResultsAll",samples[s],".csv"), row.names = F)
  
  ############ BARPLOT QUANTILES ############ 
  ### Q100
  BarQ100 = ggplot(data=Quants100, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
    geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
    geom_errorbar(width=0.3, colour="black", alpha=0.6, size=0.6)+
    scale_fill_manual(values = c(palbase[1:(length(Quants100$model)-4)],palmod))+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    # ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
    theme_light(base_size = textsize)+
    ggtitle(paste0("Q100 - Sample n°", s))+
    # ggtitle(paste0("Q100 - Echantillon n°", s))+
    scale_x_discrete(labels=c(rep("GEV",length(Quants100$model)-4),models))
    # scale_x_discrete(labels=c(rep("GEV",length(Quants100$model)-4),c("u0","uS","uN","uSN")))
  
  BarQ1000 = ggplot(data=Quants1000, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
    geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
    geom_errorbar(width=0.3, colour="black", alpha=0.6, size=0.6)+
    scale_fill_manual(values = c(palbase[1:(length(Quants100$model)-4)],palmod))+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    # ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
    theme_light(base_size = textsize)+
    ggtitle(paste0("Q1000 - Sample n°", s) )+
    # ggtitle(paste0("Q1000 - Echantillon n°", s) )+
    scale_x_discrete(labels=c(rep("GEV",length(Quants100$model)-4),models))
    
  ggarrange(BarQ100+coord_cartesian(ylim = c(5000,max(Quants1000$Q_9)) )+
              theme(axis.title.x = element_blank(),legend.title = element_blank()),
            BarQ1000+coord_cartesian(ylim = c(5000,max(Quants1000$Q_9)) )+
              theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                    legend.title = element_blank()),
            common.legend = T, legend = "right")
  
  ggsave(path = paste0(dir.plots,samples[s]), filename = paste0("Barplots_QX_",samples[s],".pdf"),
         width = 13,height = 8)
  
  
  ############ PLOT POSTERIOR PARAMETERS ############ 
  
  ####  SHAPE
  GGShape = ggplot(Params[which(Params$variable=="Shape"),], 
                 aes(y = model, x = value, group = model, fill = model)) +
    geom_density_ridges(alpha=0.7,stat = "binline",bins = 60, scale = 1.2,color = NA)+
    geom_segment(data = ParamsMp[
                  which(ParamsMp$variable == "Form"),][rev(c(5:length(Quants100$model),1,2,3,4)),],
                 aes(x = value, xend = value,
                     y=(1:length(Quants100$model)),yend=(2:(length(Quants100$model)+1)),color = model)
                 , lwd = 1)+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_discrete(expand = c(0, 0),
                     limits = rev(as.vector(unique(ParamsMp$model)[c(5:length(Mp_par),1,2,3,4)]) ),
                     labels = rev(c(rep("GEV",(length(Quants100$model)-4)),(models))) )+
    theme_light(base_size = textsize)+
    scale_fill_manual(values = c(palbase[1:(length(Quants100$model)-4)],palmod) )+
    scale_color_manual(values = (c(palbase[1:(length(Quants100$model)-4)],palmod) ) )+
    xlab("[-]")+
    theme(axis.title.y = element_blank(),legend.title = element_blank())+
    # ggtitle("GEV Shape parameter")+
    # ggtitle("Paramètre de forme \u03bE")+
    ggtitle("Shape parameter \u03bE")+
    coord_cartesian(xlim=c(0-3*sd(Params$value[which(Params$variable=="Shape")]),
                           0+4*sd(Params$value[which(Params$variable=="Shape")])))
  
  
  ggsave(path = paste0(dir.plots,samples[s]), filename = paste0("Shape_",samples[s],".pdf"),
         device = cairo_pdf, width = 10,height = 8)
  
  #### Threshold
  priorthres = Params[which(Params$variable=="Thres" & Params$model==paste0(models[1],shortm[1])
                            & Params$sample==samples[s]),]
  priorthres$value = rnorm(length(priorthres$model),Thres[s],sdThres)
  priorthres$model = "1-Prior"
  ParamsMp$value[which(ParamsMp$variable == "Seuil")][which(distsTs == "FIX")] = 0
  for(i in 1:length(models[which(distsTs == "FIX")])){Params$value[which(Params$variable=="Thres" &
                       Params$model== paste0(models,shortm)[which(distsTs == "FIX")][i])] = 0}
  GGThres = ggplot( rbind(priorthres,Params[which(Params$variable == "Thres" ),]), 
                  aes(x = value, fill = model)) +
    
    geom_histogram(aes(y=..density..),position = "identity", bins = 300, alpha = 0.7, color = NA)+
    geom_density(data=priorthres, aes(x = value),fill = NA, lty = 2)+
    geom_vline(data = ParamsMp[which(ParamsMp$variable == "Seuil"),],
               aes( xintercept = value, color = model), lwd = 1, alpha = 1,show.legend = F)+
    scale_fill_manual(values = c(colprior, palmod) )+
    scale_color_manual(values = palmod)+
    coord_cartesian(ylim=c(0,1.5*max(density(Params$value[which(Params$variable=="Thres" &
        Params$sample==samples[s] & Params$model=="B: uS")])$y)),
        xlim=c(Thres[s]-2*sdThres, Thres[s]+2*sdThres))+
    theme_light(base_size = textsize)+
    geom_vline(xintercept = Thres[s], lwd = 1, lty = 2)+
    xlab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    ylab("Density")+
    ggtitle("Perception threshold")+
    # xlab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
    # ylab("Densité")+
    # ggtitle("Seuil de perception S")+
    theme(legend.title = element_blank())

  #### Historical period length
  priornban = Params[which(Params$variable=="Nban" & Params$model==paste0(models[1],shortm[1])
                            & Params$sample==samples[s]),]
  priornban$value = runif(length(priornban$model),NbanS[s],MaxUnif_Nban)
  priornban$model = "1-Prior"
  for(i in 1:length(models[which(distsNb == "FIX")])){Params$value[which(Params$variable=="Nban" &
                    Params$model==paste0(models,shortm)[which(distsNb == "FIX")][i])] = 0  }
  ParamsMp$value[which(ParamsMp$variable=="NBan")][which(distsNb == "FIX")] = 0
  
  # endPhist = Qconts[[s]]$an[1]
  #### CALCULATION OF START DATE FROM ILARIA PRODOSCIMI
    ## number of floods histo
    k = length(Occs[[s]])
    ## differences between time of occurrence of historical events and start of systematic records
    tk = (endPhist) - Occs[[s]][1]
    T_Ilaria = (endPhist) - round(tk + tk/k - 1)
  ## Calcul TGEV(S) Michel
    T_RpThres = Occs[[s]][1] - round(1/(1-prob[which.min(abs(QuantsAll[[s]]$Mp - Thres[s]))]))
  
  GGNban = ggplot(rbind(priornban,Params[which(Params$variable == "Nban"),]), 
                aes(x = endPhist-value, fill = model)) +
    geom_histogram(aes(y=..density..),position = "identity", bins = 600, alpha = 0.7,color = NA)+
    geom_density(data=priornban, aes(x = endPhist-value), adjust = 0.4, fill = NA, lty = 2)+
    geom_vline(data = ParamsMp[which(ParamsMp$variable == "NBan"),],
               aes(xintercept = endPhist-value, color = model), lwd = 1, alpha = 1,
               show.legend = T, lty=1)+
    geom_vline(aes(xintercept = StartH[s], color = "Reference"), lwd = 1, lty = 2)+
    geom_vline(aes(xintercept = T_Ilaria, color = "Prosdocimi"), lwd = 1, lty = 2)+
    geom_vline(aes(xintercept = T_RpThres, color = "T_GEV(S)"), lwd = 1, lty = 2)+
    scale_fill_manual(values = c(colprior, palmod) )+
    scale_color_manual(values = c(palmod,"#31a354","black","#e7298a"))+
    coord_cartesian(ylim = c(0, 1.1*max(density(Params$value[which(Params$variable=="Nban" &
                              Params$sample==samples[s] & Params$model==paste0(models,shortm)[3])])$y)),
                    # xlim = c( endPhist-max(ParamsMp$value[which(ParamsMp$variable=="NBan")])-30,
                    #           endPhist-NbanS[s]+2 )  )+
                    xlim = c( endPhist-(max(priornban$value)+2),
                              endPhist-NbanS[s]+2 )  )+
    theme_light(base_size = textsize)+
    xlab("Year")+
    ylab("Density")+
    ggtitle("Start date of historical period")+
    # xlab("Année")+
    # ylab("Densité")+
    # ggtitle("Date de début de la période historique t*")+
    theme(legend.title = element_blank())
    
  ## ARRANGE
  ggarrange(GGThres,GGNban+theme(axis.title.y = element_blank()),
            ncol = 2, common.legend = T, legend = "bottom", align = "hv")
  ggsave(path = paste0(dir.plots,samples[s]), filename = paste0("Params_",samples[s],".pdf"),
         width = 15,height = 7)
  
  ######## QUANTILES
  ymax = max( max(QuantsAll[[2]]$Q_9), max(QuantsAll[[length(QuantsAll)]]$Q_9) )
  for(b in 1:(length(Quants1000$model)-4)){
    Freq = Qbaselines[[b]][order(Qbaselines[[b]]$mp),]
    Freq$Fr = (seq(1:length(Qbaselines[[b]]$an))-0.5)/length(Qbaselines[[b]]$an)
    Freq$Pr = 1/(1-Freq$Fr)
    
    GGQuants[[b]] = ggplot()+
      geom_ribbon(data=QuantsAll[[length(models)+b]],aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                             fill="95% uncertainty interval"),alpha=0.8)+
                                             # fill="Intervalle de conf. 95%"),alpha=0.8)+
      geom_line(data=QuantsAll[[length(models)+b]],aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
      geom_point(data = Freq, aes(x=Pr, y = mp))+
      geom_errorbar(data=Freq, aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
      scale_x_continuous(trans="log10")+
      xlab("Return period [years]")+
      ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
      # xlab("Période de retour [années]")+
      # ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
      scale_fill_manual(name = element_blank(),values = c("#67a9cf"))+
      scale_color_manual(values = c("yellow"))+ 
      theme_bw(base_size=textsize)+
      coord_cartesian(xlim=c(1,max(Pr)),ylim=c(2000, ymax) )+
      theme(legend.title=element_blank(),
            plot.title = element_text(hjust = 0.01, vjust = -7),
            legend.position = c(0.8,0.2))+
      ggtitle(paste0("GEV",Baselines[b]))
  }
  
  for(m in 1:length(models)){

    #### plot positions histo
    alpha = 0.5
    # Durée periode continue
    NC = length(Qconts[[s]]$an)
    # Durée periode histo
    NH = Nhists[[s]]
    #Nombre de crues histo sup. au seuil
    NS_H = length(Occs[[s]])
    #Nombre de crues de la periode continue sup. au seuil
    NS_C = length(which(Qconts[[s]]$mp > Thres[s]))
    #Nombre tot. de crues sup. au seuil 
    NS = NS_H + NS_C
    #Nombre d'annees total
    NAn = NC + NH
    ############ Sup seuil
    SupSeuil = data.frame(rang=1:NS)
    SupSeuil$flag = "Hist"
    # On tire exactement NS_C crues "cont"
    set.seed(2020)
    SupSeuil$flag[sample(1:NS,NS_C,replace = F)] = "Cont"
    SupSeuil$mp = NA
    SupSeuil$tot2.5 = NA
    SupSeuil$tot97.5 = NA
    SupSeuil$mp[which(SupSeuil$flag == "Cont")] = Qconts[[s]][,2][order(Qconts[[s]][,2],
                                                                        decreasing = T)]
    SupSeuil$tot2.5[which(SupSeuil$flag == "Cont")] = Qconts[[s]][,3][order(Qconts[[s]][,2],
                                                                            decreasing = T)]
    SupSeuil$tot97.5[which(SupSeuil$flag == "Cont")] = Qconts[[s]][,4][order(Qconts[[s]][,2],
                                                                       decreasing = T)]
    # Freq HIRSCH 
    SupSeuil$Fr = (NS/NAn) * ((SupSeuil$rang-alpha)/(NS+1-2*alpha))
    ############ Sous seuil
    SousSeuil = Qconts[[s]][which(Qconts[[s]][,2]<=Thres[[s]]),(2:4)]
    SousSeuil = SousSeuil[order(SousSeuil$mp, decreasing = T),]
    # Freq HIRSCH
    SousSeuil$rang = (NS+1):(NS_H+NC)
    SousSeuil$flag = "Cont"
    # Freq HIRSCH
    SousSeuil$Fr = (NS/NAn) + (1-NS/NAn) * ( ((SousSeuil$rang-NS)-alpha)/(NC-NS_C+1-2*alpha) )
    ## Tab freq empirique
    Freq = merge.data.frame(SupSeuil,SousSeuil, all = T)
    Freq$tot2.5[which(Freq$flag == "Hist")] = Thres[s]
    Freq$tot97.5[which(Freq$flag == "Hist")] = 50000
    Freq$Pr = round(1/Freq$Fr,4)
    rhist = which(Freq$flag=="Hist")
    rcont = which(Freq$flag=="Cont")
    cptcont = 1
    for(i in 1:length(rhist)){
      while(rhist[i]-rcont[cptcont] > 0){cptcont=cptcont+1}
      Freq$tot2.5[rhist[i]] = Freq$mp[rcont[cptcont]] }

    ### PLOT QUANTS
    GGQuants[[m+(length(Quants1000$model)-4)]] = ggplot()+
      geom_ribbon(data=QuantsAll[[m]],aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                             fill="95% uncertainty interval"),alpha=0.8)+
                                          # fill="Intervalle de conf. 95%"),alpha=0.8)+
      geom_line(data=QuantsAll[[m]],aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
      geom_errorbar(data=Freq[which(Freq$flag=="Cont"),],
                    aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
      geom_point(data = Freq[which(Freq$flag=="Cont"),],
                 aes(x=Pr, y = mp))+
      geom_errorbar(data=Freq[which(Freq$flag=="Hist"),],
                    aes(x=Pr,ymin = tot2.5, ymax = tot97.5),lty=2, width=0,
                    col="lightgrey")+
      geom_point(data = Freq[which(Freq$flag=="Hist"),], aes(x=Pr,y=tot2.5),shape = 17,
                 col="lightgrey")+
      scale_x_continuous(trans="log10")+
      xlab("Return period [years]")+
      ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
      # xlab("Période de retour [années]")+
      # ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
      scale_fill_manual(name = element_blank(),
                        values = c("#67a9cf"))+
      scale_color_manual(values = c("yellow"))+
      theme_bw(base_size=textsize)+
      coord_cartesian(xlim=c(1,max(Pr)),ylim = c(2000, ymax ))+
      theme(legend.title=element_blank(),
            plot.title = element_text(hjust = 0.01, vjust = -7),
            legend.position = c(0.8,0.2))+
      ggtitle(paste0("Model ",models[m]))
      # ggtitle(paste0("Modèle ",models[m]," ",shortm[m]))
    
  }
  
  if(length(GGQuants) == 5){
    ggarrange(GGQuants[[1]],
              GGQuants[[2]],
              GGQuants[[3]],
              GGQuants[[4]],
              GGQuants[[5]],
              ncol = 2, nrow = 3,
              common.legend = T, legend = "bottom",align = "hv")
    
  } else {
    ggarrange(GGQuants[[1]],
              GGQuants[[2]],
              GGQuants[[3]],
              GGQuants[[4]],
              GGQuants[[5]],
              GGQuants[[6]],
              ncol = 2, nrow = 3,
              common.legend = T, legend = "bottom",align = "hv")
    
  }
  ggsave(path = paste0(dir.plots,samples[s],"/"), filename = paste0("Quantiles_",samples[s],".pdf"),
         width = 18, height = 14)
    
  
  ## SCATTERPLOTS
  npoints = length(ParamAll[[m]][,1])
  # ggarrange(ggscatmat(ParamAll[[1]][ sample(x = 1:npoints, size = 10000, replace = F) , (1:3)]),
  #           ggscatmat(ParamAll[[2]][ sample(x = 1:npoints, size = 10000, replace = F) , (1:3)]),
  #           ggscatmat(ParamAll[[3]][ sample(x = 1:npoints, size = 10000, replace = F) , (1:3)]),
  #           ggscatmat(ParamAll[[4]][ sample(x = 1:npoints, size = 10000, replace = F) , (1:3)]),
  #           ncol = 2, nrow = 2,align = "hv", 
  #           labels = models)
  cairo_pdf(paste0(dir.plots,samples[s],"/","ScatterD_",samples[s],".pdf"), height = 6, width=6)
    pairs.panels(ParamAll[[4]][ sample(x = 1:npoints, size = 10000, replace = F) , (1:3)],
                 hist.col = palmod[4], density = T,lm = F,smooth = F, labels = c("\u03bE","S","n"))
  dev.off()
  # ggsave(filename = paste0("ScatterD_",samples[s],".pdf"),path = paste0(dir.plots,samples[s],"/"),
  #        height = 10, width = 10)
  
  
}
beep()


