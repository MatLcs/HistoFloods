rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
#case folder
case = "HistoSliding"
dir.res.case = paste0(dir.res,case,"/"); dir.create(dir.res.case,showWarnings = F)
dir.plot.case = paste0(dir.plots,case,"/"); dir.create(dir.plot.case,showWarnings = F)
#### Run or plot only ?
RUN = T
#### Threshold
Ts = 8000
sd.Ts = 1000
#### PARAMS SIMULATION
# Ychange = 1970 ;
Start = 1891 ; End = 2000
Nsim = 2000 ; Nspag = 200 ; Andeb = 1816  #Anfin = 1970
# RP up top 10 000 years
# prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
prob = c(0.99,0.999) ; Pr = 1/(1-prob)
#### Color palette
pal = c("grey",brewer.pal(4,"RdYlBu")[c(2,1,3,4)])
#### Results reading init
Quant100 = list();Quant1000 = list();ParamAll = list();Mp_par = list();QuantsAll = list()
#GEV params
Pos = parameter(name='Pos',init = 6000) 
Ech =  parameter(name='Ech',init = 1000) 
Form = parameter(name='Form',init = 0.01,priorDist='Gaussian',priorPar=c(0,0.3))
#Data
SpagsAll = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
QAll = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
                  header = T)[,c(1,2,7,8)]
require(progress)
pb = txtProgressBar(min=0,max=End-Start,style=2)
#### RUN 
if(RUN == T){
    #while the number of historical floods is > to 2 floods
    for(Yhisto in Start:End){
      dir.y = paste0(dir.res.case,Yhisto,"/"); dir.create(dir.y, showWarnings = F)
      dir.case = dir.y
      #histofloods
      CX = QAll[which(QAll$mp > Ts & QAll$an <= Yhisto),]
      names(CX) = c("An","mp","tot2.5","tot97.5")
      Anfin = Yhisto
      #Systematic floods
      Q = QAll[which(QAll$an > Yhisto),]
      Spags = SpagsAll[which(QAll$an > Yhisto),]
      #stoods execution
      source(paste0(dir.codes,"B1_Sliding.R"))  
      beep()
      setTxtProgressBar(pb, Yhisto)
  }
}


