rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

###### RUN PARAMS
Nsim = 8000; Nspag = 500 ; textsize = 15
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
starthist = 1500
endhist = 1815

#### Run or plot only ?
RUN = F

## DIRS
case = "C4short"
dir.res.case = paste0(dir.res,"/",case) ; dir.create(dir.res.case)
dir.plot.case = paste0(dir.plots,"/",case) ; dir.create(dir.plot.case)

###### DATA LOADING
###Histrhone data 
Call = read.csv2(paste0(dir.data,"CX_All.csv"))
### 1816-2020 AMAX discharges with uncertainty and spags
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
Q = read.table(
  "C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]

samples =  "C4"
Occs =  Call[which(Call$An>=starthist & Call$An<=endhist & Call$Cat == "C4"), 1]
Nhists = endhist-starthist+1
#### Threshold priors
Thres = 9000
sdThres = 500
#### Nban priors (when Nban is draw in Unif dist.) : end of histo. period - date of first OT flood 
NbanS = endhist-Call$An[which(Call$An > 1500 & Call$Cat == "C4")][1]+1

StartH = starthist

MaxUnif_Nban = (Occs[1] - starthist)*2
#### 4 models
models =  "D"
shortm = ": uS&t*"
distsTs = "Gaussian"
distsNb = "Uniform"
#### Baselines
# Baselines =     c("1816-2020","1970-2020")
# Qbaselines =    list(Q, Q[which(Q$an>endhistartif),])
# SpagsBaseline = list(Spags[, (1:Nspag) ], Spags[which(Q$an>endhistartif), (1:Nspag) ])

#### RUN MODELS
if(RUN == T){
  #### Run GEV_Binom
  # Fun run
  GEV_Binom(Qcont = Q,
            SpagCont = Spags,
            Occ = data.frame(an=Occs),
            Ts = list(distsTs,c(Thres,sdThres),Thres),#dist,prior,init
            Nban = list(distsNb,c(NbanS,NbanS+MaxUnif_Nban),NbanS),#dist,prior,init
            prob = prob,
            Nsim = Nsim,
            Nhist = Nhists,
            dir.res = dir.res.case)

}



ind100 = length(Pr) - 18
ind1000 = length(Pr) - 9

#### Results reading init
Quant100 = list(); Quant1000 = list(); ParamAll = list(); Mp_par = list(); QuantsAll = list()
palmod = brewer.pal(4,"RdYlBu")[4] 
palbase = "#4dac26"
colprior = "lightgrey"
  

######## RES GEV BINOM
## All quants
QuantsAll[[1]] = read.table(paste0(dir.res.case,"/Quants.txt"), header = T)[,(2:4)]
QuantsAll[[1]]$model = paste0(models,shortm); QuantsAll[[1]]$sample = samples
##Q100 & Q1000
Quant100[[1]] = read.table(paste0(dir.res.case,"/Quants.txt"), header = T)[ind100,(2:4)]
Quant100[[1]]$model = paste0(models,shortm); Quant100[[1]]$sample = samples
Quant1000[[1]] = read.table(paste0(dir.res.case,"/Quants.txt"), header = T)[ind1000,(2:4)]
Quant1000[[1]]$model = paste0(models,shortm); Quant1000[[1]]$sample = samples
## posterior params
ParamAll[[1]] = read.table(paste0(dir.res.case,"/Params.txt"), header = T)
ParamAll[[1]]$model = paste0(models,shortm); ParamAll[[1]]$sample = samples
## maxpost params
Mp_par[[1]] = read.table(paste0(dir.res.case,"/Maxpost_Par.txt"), header = T)
Mp_par[[1]]$model = paste0(models,shortm); Mp_par[[1]]$sample = samples


######## RES GEV BASELINE
dir.base = paste0(dir.res,"Baseline_1816-2020")
## All quants
# QuantsAll[[m]] = read.table(paste0(dir.case,"/Quants.txt"), header = T)[,(2:4)]
# QuantsAll[[m]]$model = paste0("GEV",Baselines[b]); QuantsAll[[m]]$sample = Baselines[b]
##Q100 & Q1000
Quant100[[2]] = read.table(paste0(dir.base,"/Quants.txt"), header = T)[ind100,(2:4)]
Quant100[[2]]$model = paste0("GEV_1816-2020"); Quant100[[2]]$sample = "1816-2020"
Quant1000[[2]] = read.table(paste0(dir.base,"/Quants.txt"), header = T)[ind1000,(2:4)]
Quant1000[[2]]$model = paste0("GEV_1816-2020"); Quant1000[[2]]$sample = "1816-2020"
## posterior params
ParamAll[[2]] = read.table(paste0(dir.base,"/Params.txt"), header = T)
ParamAll[[2]]$model = paste0("GEV_1816-2020"); ParamAll[[2]]$sample = "1816-2020"
## maxpost params
Mp_par[[2]] = read.table(paste0(dir.base,"/Maxpost_Par.txt"), header = T)
Mp_par[[2]]$model = paste0("GEV_1816-2020"); Mp_par[[2]]$sample ="1816-2020"

## RESULTS TIDYING
Params = melt.list(ParamAll) ; Params$model = as.factor(Params$model)
Quants100 = cast(melt.list(Quant100) );  Quants1000 = cast(melt.list(Quant1000) )
ParamsMp = melt.list(Mp_par)
Quants100$model = as.factor(Quants100$model);Quants1000$model = as.factor(Quants1000$model)
Quants100$model =fct_relevel(Quants100$model,
                             as.vector(Quants100$model[c(2,1)]) )
Quants1000$model=fct_relevel(Quants1000$model,
                             as.vector(Quants1000$model[c(2,1)]) )

 
#### DF RESUME RESULTS
ResAll = data.frame(model = c(models,"GEV_1816-2020"),
                    sd.S=NA, mp.S=NA, sd.t = NA, mp.t = NA,
                    sd.xi = NA, mp.xi = NA,  sd.Q100 = NA, mp.Q100 = NA,
                    sd.Q1000 = NA, mp.Q1000 = NA)

endPhist = Q$an[1]
m=1
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

### BASELINES
b=2
  ### xi
  ResAll[b,]$sd.xi = round(sd(Params$value[which( Params$model == "GEV_1816-2020" &
      Params$variable == "Shape") ] ),3)
  ResAll[b,]$mp.xi = round(ParamsMp$value[which( ParamsMp$model == "GEV_1816-2020" &
      ParamsMp$variable == "Form") ],3)
  ### Q100
  ResAll[b,]$sd.Q100 = round( ((Quants100$Q_9 - Quants100$Mp)/2)[b] )
  ResAll[b,]$mp.Q100 = round( Quants100$Mp[b] )
  #Q1000
  ResAll[b,]$sd.Q1000 = round( ((Quants1000$Q_9 - Quants1000$Mp)/2)[b] )
  ResAll[b,]$mp.Q1000 = round( Quants1000$Mp[b] )

### write resume results
write.csv2(ResAll, paste0(dir.res.case,"/ResultsAll",case,".csv"), row.names = F)


############ BARPLOT QUANTILES ############ 
### Q100
BarQ100 = ggplot(data=Quants100, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.3, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = c(palbase,palmod))+
  # ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
  theme_light(base_size = textsize)+
  ggtitle(paste0("Q100 - Echantillon n° 2"))+
  scale_x_discrete(labels=c('GEV',models))

BarQ1000 = ggplot(data=Quants1000, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.3, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = c(palbase,palmod))+
  # ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
  theme_light(base_size = textsize)+
  ggtitle(paste0("Q1000 - Echantillon n° 2"))+
  scale_x_discrete(labels=c('GEV',models))

ggarrange(BarQ100+coord_cartesian(ylim = c(5000,max(Quants1000$Q_9)) )+
            theme(axis.title.x = element_blank(),legend.title = element_blank()),
          BarQ1000+coord_cartesian(ylim = c(5000,max(Quants1000$Q_9)) )+
            theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                  legend.title = element_blank()),
          common.legend = T, legend = "right")

ggsave(path = paste0(dir.plot.case), filename = paste0("Barplots_QX_",case,".pdf"),
       width = 13,height = 8)


############ PLOT POSTERIOR PARAMETERS ############ 
####  SHAPE
GGShape = ggplot(Params[which(Params$variable=="Shape"),], 
                 aes(y = model, x = value, group = model, fill = model)) +
  geom_density_ridges(alpha=0.7,stat = "binline",bins = 60, scale = 1.2,color = NA)+
  geom_segment(data = ParamsMp[which(ParamsMp$variable == "Form"),],
               aes(x = value, xend = value,
                   y=(1:length(Quants100$model)),yend=(2:(length(Quants100$model)+1)),color = model)
               , lwd = 1)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0)) +
  theme_light(base_size = textsize)+
  scale_fill_manual(values = c(palmod,palbase) )+
  scale_color_manual(values = c(palmod,palbase) )+
  xlab("[-]")+
  theme(axis.title.y = element_blank(),legend.title = element_blank())+
  # ggtitle("GEV Shape parameter")+
  ggtitle("Paramètre de forme \u03bE")+
  coord_cartesian(xlim=c(0-2*sd(Params$value[which(Params$variable=="Shape")]),
                         0+5*sd(Params$value[which(Params$variable=="Shape")])))

ggsave(path = paste0(dir.plot.case), filename = paste0("Shape_",case,".pdf"),
       device = cairo_pdf, width = 10,height = 8)


#### Threshold
priorthres = Params[which(Params$variable=="Thres" & Params$model==paste0(models[1],shortm[1])
                          & Params$sample==samples),]
priorthres$value = rnorm(length(priorthres$model),Thres,sdThres)
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
                                            Params$sample==samples & Params$model=="D: uS&t*")])$y)),
                  xlim=c(Thres-2.5*sdThres, Thres+3*sdThres))+
  theme_light(base_size = textsize)+
  geom_vline(xintercept = Thres, lwd = 1, lty = 2)+
  # xlab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  # ylab("Density")+
  # ggtitle("Perception threshold")+
  xlab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
  ylab("Densité")+
  ggtitle("Seuil de perception S")+
  theme(legend.title = element_blank())



#### Historical period length
priornban = Params[which(Params$variable=="Nban" & Params$model==paste0(models[1],shortm[1])
                         & Params$sample==samples),]
priornban$value = runif(length(priornban$model),NbanS,NbanS+MaxUnif_Nban)
priornban$model = "1-Prior"
for(i in 1:length(models[which(distsNb == "FIX")])){Params$value[which(Params$variable=="Nban" &
                               Params$model==paste0(models,shortm)[which(distsNb == "FIX")][i])] = 0  }
ParamsMp$value[which(ParamsMp$variable=="NBan")][which(distsNb == "FIX")] = 0

endPhist = Q$an[1]
#### CALCULATION OF START DATE FROM ILARIA PRODOSCIMI
## number of floods histo
k = length(Occs)
## differences between time of occurrence of historical events and start of systematic records
tk = (endPhist) - Occs[1]
T_Ilaria = (endPhist) - round(tk + tk/k - 1)
## Calcul TGEV(S) Michel
T_RpThres = Occs[1] - round(1/(1-prob[which.min(abs(QuantsAll[[1]]$Mp - Thres))]))

GGNban = ggplot(rbind(priornban,Params[which(Params$variable == "Nban"),]), 
                aes(x = endPhist-value, fill = model)) +
  geom_histogram(aes(y=..density..),position = "identity", bins = 600, alpha = 0.7,color = NA)+
  geom_density(data=priornban, aes(x = endPhist-value), adjust = 0.4, fill = NA, lty = 2)+
  geom_vline(data = ParamsMp[which(ParamsMp$variable == "NBan"),],
             aes(xintercept = endPhist-value, color = model), lwd = 1, alpha = 1,
             show.legend = T, lty=1)+
  geom_vline(aes(xintercept = StartH, color = "Référence"), lwd = 1, lty = 2)+
  geom_vline(aes(xintercept = T_Ilaria, color = "Prosdocimi"), lwd = 1, lty = 2)+
  geom_vline(aes(xintercept = T_RpThres, color = "T_GEV(S)"), lwd = 1, lty = 2)+
  scale_fill_manual(values = c(colprior, palmod) )+
  scale_color_manual(values = c(palmod,"#31a354","black","#e7298a"))+
  coord_cartesian(ylim = c(0, 1.1*max(density(Params$value[which(Params$variable=="Nban" &
                                    Params$sample==samples & Params$model==paste0(models,shortm))])$y)),
                  # xlim = c( endPhist-max(ParamsMp$value[which(ParamsMp$variable=="NBan")])-30,
                  #           endPhist-NbanS[s]+2 )  )+
                  xlim = c( endPhist-(max(priornban$value)+5),
                            endPhist-NbanS+5 )  )+
  theme_light(base_size = textsize)+
  # xlab("Year")+
  # ylab("Density")+
  # ggtitle("Start date of historical period")+
  xlab("Année")+
  ylab("Densité")+
  ggtitle("Date de début de la période historique t*")+
  theme(legend.title = element_blank())

## ARRANGE
ggarrange(GGThres,GGNban+theme(axis.title.y = element_blank()),
          ncol = 2, common.legend = T, legend = "bottom", align = "hv")

ggsave(path = dir.plot.case, filename = paste0("Params_",case,".pdf"),
       width = 15,height = 7)



