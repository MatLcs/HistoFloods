rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001))

###Histrhone data since 1500
Call = data.frame(an=read.csv2(paste0(dir.data,"CX_All.csv"))[,1])
### 1816-2020 AMAX discharges with uncertainty and spaghetti 
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
Q = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]

Ts = list("Gaussian",c(7000,2000),7000)
Nban = list("Uniform",c(1816-1300,1000),1816-1300)

# 
# Qcont = Q
# SpagCont = Spags[,c(1:10)]
# Occ = data.frame(an=Call$an)
# Nsim = 100
# dir.res = dir.test

GEV_Binom = function(Qcont,   #df with an, mp, tot2.5, tot97.5 of continuous period amax discharges
                     SpagCont,#matrix of hydrometric uncertainty of continuous period amax
                     Occ,     #dates of OT occurences of historic period
                     Ts,      #list with dist, prior, init
                     Nban,    #list with dist, prior, init
                     prob,    #quantiles to extract
                     Nsim,    #number of MCMC sim
                     dir.res){#dir to run and write results
  
  if(Ts[[1]] == "FIX"){
    Seuil = parameter(name = "Seuil", init = Ts[[3]], priorDist = "FIX")
  } else {Seuil = parameter(name = "Seuil", init = Ts[[3]], priorDist = Ts[[1]], priorPar = Ts[[2]])  }
  
  if(Nban[[1]] == "FIX"){
    NBan = parameter(name = "NBan", init = Nban[[3]], priorDist = "FIX")
  } else {
    NBan = parameter(name = "NBan", init = Nban[[3]], priorDist = Nban[[1]], priorPar = Nban[[2]])  }

  Pos = parameter(name='Pos',init = 6000)  
  Ech =  parameter(name='Ech',init = 1000) 
  Form = parameter(name='Form',init = 0.01,priorDist='Gaussian',priorPar=c(0,0.3))
  
  MegaSpag = data.frame()
  Nspag = ncol(SpagCont)
  for(spag in 1:Nspag){
    print(paste0("Spag = ",spag))
    Y = data.frame(  c(SpagCont[,spag], length(Occ[,1])) )
    Var = data.frame(factor(c(rep('Q',length(Qcont[,1])),"Occ")))
    dat <- dataset(Y = Y, var = Var)
    mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
                 par=list(Pos,Ech,Form,Seuil,NBan),
                 ### PARAM GEV
                 formula = c("Pos=Pos","Ech=Ech","Form = Form",
                             ### PARAM BINOMIAL
                             "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                             "NBan = NBan"))
    #### Create spag subfolder 
    dir.spag = file.path(dir.res,"Spag"); dir.create(dir.spag,showWarnings = F)
    STooDs(mod,workspace = dir.spag, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
    #### Read MCMC results for each spag
    MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
    MegaSpag = rbind(MegaSpag, MCMCres[,(1:5)])
  }
  MegaGev = matrix(nrow = nrow(MegaSpag),ncol = length(prob))
  #### Compute the quantiles for all the GeV realisations spags
  #### WARNING, HERE SHAPE HAS TO BE THE OPPOSITE THE SHAPE GIVEN BY STOODS (BY CONVENTION)
  for ( i in 1 : nrow(MegaSpag)) {MegaGev[i,] = qgev(p = prob,loc = MegaSpag$Pos[i],
                                                     scale = MegaSpag$Ech[i], 
                                                     shape = -1*MegaSpag$Form[i]) }
  PostForm = MegaSpag$Form
  PostThres = MegaSpag$Seuil
  PostNban = MegaSpag$NBan
  #### Compute the true maxpost quantiles : maxpost of hydro sample x maxpost of GeV estim
  Y = data.frame(c(Qcont[,2], length(Occ[,1])) )
  Var = data.frame(factor(c(rep('Q',length(Qcont[,2])),"Occ")))
  dat <- dataset(Y = Y, var = Var)
  #### Create MP sub-folder & run Stoods
  dir.mp = file.path(dir.res,"/Mp") ; dir.create(dir.mp,showWarnings = F)
  dat <- dataset(Y = Y, var = Var)
  mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
               par=list(Pos,Ech,Form,Seuil,NBan),
               ### PARAM GEV
               formula = c("Pos=Pos","Ech=Ech","Form = Form",
                           ### PARAM BINOMIAL
                           "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                           "NBan = NBan"))
  STooDs(mod,workspace = dir.mp, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
  MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
  Mp.GeV = MCMCres[which.max(MCMCres$post),(1:5)]
  Mp.Quant = qgev(p = prob,loc = Mp.GeV$Pos, scale = Mp.GeV$Ech,shape = -1*Mp.GeV$Form)
  quant = apply(MegaGev,MARGIN = 2,FUN = quantile, probs = c(0.025,0.975))
  Quants = data.frame(Pr=1/(1-prob), Mp=Mp.Quant, Q_2=quant[1,] ,Q_9=quant[2,])
  ### write results in res dir
  write.table(round(Quants,3), paste0(dir.res,"/Quants.txt"), row.names = F)
  write.table(round(data.frame(Shape = PostForm, Thres = PostThres, Nban = PostNban),3),
              paste0(dir.res,"/Params.txt"), row.names = F)
  write.table(Mp.GeV,  paste0(dir.res,"/Maxpost_Par.txt"), row.names = F)
  
  #### plot positions histo
  NC = length(Qcont[,1])
  NH = Nban[[3]]
  #Nombre de crues histo sup. au seuil
  NS_H = nrow(Occ)
  #Nombre de crues de la periode continue sup. au seuil
  NS_C = length(which(Qcont[,2] > Ts[[3]] ))
  #Nombre tot. de crues sup. au seuil 
  NS = NS_H + NS_C
  #Nombre de crues inf. au seuil
  Ninf = length(which(Qcont[,2] <= Ts[[3]]))
  #Nombre d'annees total
  NAn = NC + NH
  #proba crue continue
  PC = NS_C/NS
  #proba crue histo
  PH = NS_H/NS
  #data frame crues supseuil
  SupSeuil = data.frame(rang=1:NS)
  SupSeuil$flag = "Hist"
  set.seed(2020)
  SupSeuil$flag[sample(1:NS,NS_C,replace = F)]="Cont"
  SupSeuil$mp = NA
  SupSeuil$tot2.5 = NA
  SupSeuil$tot97.5 = NA
  SupSeuil$mp[which(SupSeuil$flag == "Cont")] = Qcont[,2][order(Qcont[,2],decreasing = T)]
  SupSeuil$tot2.5[which(SupSeuil$flag == "Cont")] = Qcont[,3][order(Qcont[,3],decreasing = T)]
  SupSeuil$tot97.5[which(SupSeuil$flag == "Cont")] = Qcont[,4][order(Qcont[,4],decreasing = T)]
  SousSeuil = Qcont[which(Qcont[,2]<=Ts[[3]]),(2:4)]
  SousSeuil = SousSeuil[order(SousSeuil$mp, decreasing = T),]
  SousSeuil$rang = (NS+1):(NS_H+NC)
  SousSeuil$flag = "Cont"
  Freq = merge.data.frame(SupSeuil,SousSeuil, all = T)
  Freq$tot2.5[which(Freq$flag == "Hist")] =  Ts[[3]]
  Freq$tot97.5[which(Freq$flag == "Hist")] = 50000
  Freq$Fr = (Freq$rang-0.5)/(NS+Ninf)
  Freq$Pr = 1/Freq$Fr
  
  
  ### PLOT QUANTS
  Quant = ggplot()+
    geom_ribbon(data=Quants,aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                fill="95% uncertainty interval"),alpha=0.8)+
    geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
    geom_errorbar(data=Freq[which(Freq$flag=="Hist"),],
                  aes(x=Pr,ymin = tot2.5, ymax = tot97.5),lty=2,col="lightgrey",lwd=0.5)+
    geom_point(data = Freq[which(Freq$flag=="Hist"),], aes(x=Pr,y=tot2.5),shape = 17,col="lightgrey")+
    geom_errorbar(data=Freq[which(Freq$flag=="Cont"),],
                  aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
    geom_point(data = Freq[which(Freq$flag=="Cont"),],
               aes(x=Pr, y = mp))+
    scale_x_continuous(trans="log10")+
    xlab("Return period [years]")+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    scale_fill_manual(name = element_blank(),
                      values = c("#67a9cf"))+
    scale_color_manual(values = c("yellow"))+ 
    theme_bw(base_size=15)+
    coord_cartesian(xlim=c(1,max(1/(1-prob))),ylim = c(min(Quants$Q_2),max(Quants$Q_9)))+
    theme(legend.title=element_blank(),
          plot.title = element_text(hjust = 0.01, vjust = -7),
          legend.position = c(0.9,0.1))
  
  ggsave(Quant, path = dir.res, filename = paste0("Quantiles.pdf"),
         width = 12, height = 8)

}

dir.test = paste0(dir.res,"/test") ; dir.create(dir.test,showWarnings = F)
GEV_Binom(Qcont = Q, SpagCont = Spags[,c(1:10)],Occ = data.frame(an=Call$an), 
          Ts = Ts, Nban =  Nban,prob = prob,Nsim = 100, dir.res = dir.test) 



