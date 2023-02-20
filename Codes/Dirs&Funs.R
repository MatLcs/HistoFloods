require(lubridate);require(ggplot2);require(trend);require(ggpubr); require(forcats)
require(cowplot);require(RSTooDs);require(evd);require(reshape);require(ggridges);
require(RColorBrewer);require(beepr);require(GGally);require(psych)

dir.all = "C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/"
dir.data = paste0(dir.all,"Data/") 
dir.codes = paste0(dir.all,"Codes/") 
dir.plots = paste0(dir.all,"Plots/") 
dir.res = paste0(dir.all,"Results/")



#### fonction test poisson michel

#Fonction calcul proba poisson
Wkcond = function(N,K,J){
  Wkcond = choose(N,K) * (J/N)^K * (1-(J/N))^(N-K)
  return(Wkcond)} 

### Fonction test de poisson simple
#### data -> un vecteur des annees correspondant a des crues
#### Andeb -> date de debut echantillon, possiblement anterieure a la date de la premiere crue
#### Anfin -> date de fin echantillon, possiblement posterieur a la date de la derniere crue
#### alpha -> defaut 0.05 pour intervalle confiance a 90%
#### plot -> defaut TRUE, si on veut plotter ou seulement recuperer les donnees
#### Return = un data frame avec donnees d entree et donnees pour plotter les intervalles de conf
TestPois = function(data, Andeb, Anfin, alpha=0.05, plot = T){
  require(binom)
  #nombre de crues
  NC = length(data)
  #Flood number indice
  Fnb = seq_along(data)
  #data avec debut et fin
  Yflood = c(Andeb,data,Anfin)
  #dates experimentales
  Yexp = c(seq(Andeb,Anfin,(Anfin-Andeb)/NC),Anfin)
  #demi longueur d intervalle pour calcul ellipses 
  N1 = NC%/%2 + 1
  if (NC%%2 == 0){N1 = NC/2}
  N2 = NC - N1 - 1
  #initialisation epsilon
  EpsInf = rep(0,NC)
  EpsSup = EpsInf
  #pour la premiere moitie de l ellipse
  #Borne inf
  for (j in 1:N1){
    SumP = 0 
    K    = -1
    #On cherche le nombre de crues K pour lequel somme des WKCOND est egale a alpha
    while(SumP <= alpha){
      K = K+1
      SumP = SumP + Wkcond(NC,K,j)   
      if (K == 0){
        EpsInf[j] = -j}
      else {
        #interpolation entre les 2 derniers indices
        EpsInf[j] = (K-j) - (SumP - alpha)/Wkcond(NC,K,j)}
    }
  }
  #Borne sup
  for (j in 1:N1){
    SumP = 0 
    K    = -1
    while(SumP <= (1-alpha)){
      K = K+1
      SumP = SumP + Wkcond(NC,K,j)   
      EpsSup[j] = 1 + (K-j) - (SumP - (1-alpha))/Wkcond(NC,K,j)
    }
  }
  
  #seconde moitie par symetrie de l ellipse
  for (j in 1:N2){
    EpsSup[NC-j] = -EpsInf[j]
    EpsInf[NC-j] = -EpsSup[j]}
  
  Fnb = c(0,Fnb,NC)
  EpsIN = c(0,EpsInf,0)
  EpsSN = c(0,EpsSup,0)
  
  if (plot == T){
    #### PLOT
    plot(x = Yflood, y = Fnb, type='b', pch = 19, col ="blue",
         ylab = "Cumulated flood number", xlab = "Years")
    points(x = Yflood[1], y = Fnb[1], pch = 19, col = "indianred", cex = 1.5)
    points(x = tail(Yflood,1), y = tail(Fnb,1), pch = 19, col = "indianred", cex = 1.5)
    lines(x = Yexp, y = Fnb + EpsIN, lty = 2)
    lines(x = Yexp, y = Fnb + EpsSN, lty = 2)}
  
    legend(x = "topleft",legend = c("Flood event",
                                  "Startpoint/Changepoint",
                                  "95% confidence interval"),
         col = c("blue","indianred","black"),lty = c(1,0,2), pch = c(20,20,NA),
         pt.cex = c(1,2,1), cex = 0.8)
  #### data return
  DataPois = data.frame(AnneeCrue = Yflood,
                        FloodNumber = Fnb,
                        AnExp = Yexp,
                        EllipInf = Fnb+EpsIN,
                        EllipSup =Fnb+EpsSN)
  
  return(DataPois)
}


##### FONCTION TEST POISSON SUCCESSIFS (2 sous periodes)
### idem parametres TestPois et on ajoute Anrupt
### Anrupt -> date de decoupage en 2 sous echantillons, doit correspondre a une date de crue
### Return = liste des 2 sous echantillons avec les donnes pour plotter les intervalles de conf
### Une erreur apparait si un des 2 sous echantillons contient moins de 5 crues


PoissSucc = function(data,Andeb,Anfin,Anrupt=NULL,alpha=0.05){
  
  if (is.null(Anrupt) == T){return(TestPois(data,Andeb,Anfin,alpha))}
  else {
    if (length(which(data == Anrupt)) > 0){
      if (length(which(data <= Anrupt)) >= 5 & length(which(data >=  Anrupt)) > 5){
        
        P1 = TestPois(data[which(data <= Anrupt)],Andeb,Anrupt,alpha,plot = F)
        P2 = TestPois(data[which(data >  Anrupt)],Anrupt,Anfin,alpha, plot = F)
        Fnb2 = nrow(P1)-2
        #PLOT
        TestPois(data,Andeb,Anfin,alpha)
        lines(x = P1$AnExp, y = P1$EllipInf, lty = 2, col = "green")
        lines(x = P1$AnExp, y = P1$EllipSup, lty = 2, col = "green")
        lines(x = P2$AnExp, y = P2$EllipInf + Fnb2, lty = 2, col = "green")
        lines(x = P2$AnExp, y = P2$EllipSup + Fnb2, lty = 2, col = "green")
        return(list(P1,P2))
      }
      else{print("Error ! One of the sub-periods has less than 5 floods")}
    }
    else{print("Error ! The choosen changepoint date isn't corresponding to a flood date")}
  }
}


##### FONCTION TEST POISSONS SUCCESSIFS (n sous periodes)
##### memes parametres que precedemment, mais Anrupt devient une liste pouvant contenir n elements
##### Si Anrupt n est pas renseigne l invervalle est trace pour la periode entiere
##### Anrupt --> list(Date1,Date2,...,DateN), par defaut Anrupt = NULL
##### Return = Une liste composee des elements de chaque sous periode pour plotter les intervalles
##### Une erreur apparait si Anrupt n est pas une liste, ou si une des sous periodes contient 
##### moins de 5 crues
##### Sorties graphiques : un plot pour chaque sous periode et un plot global


PoissoN = function(data,Andeb,Anfin,Anrupt=list(NULL),alpha=0.05){
  
  #Anrupt est une liste ?
  if(is.list(Anrupt)==F) {print("Anrupt n'est pas une liste");return()}
  
  #Si pas de rupture
  if(length(Anrupt) == 1  &  is.null(Anrupt[[1]]) == T){
    return(TestPois(data = data, Andeb = Andeb, Anfin = Anfin, alpha = alpha, plot = T)) }
  
  #Si n ruptures
  if(length(Anrupt) >= 1){
    #On ajoute la date finale
    Anrupt[[(length(Anrupt)+1)]] = Anfin
    #On intialise la liste contenant les infos de chaque sous periode
    VecDat = vector("list",length = length(Anrupt))
    #On boucle sur toutes les sous periodes
    for (i in 1:length(Anrupt)){
      #pour la premiere periode
      if (i == 1){
        #plus de 5 crues dans la sous-periode ?
        if( (length(which(data <= Anrupt[[i]]))) <= 5 ){
          print("Error ! One of the sub-periods has less than 5 floods");return()}
        #si oui on applique la fonction TestPois a la sous periode
        VecDat[[i]] = TestPois(data = data[which(data <= Anrupt[[i]])],
                               Andeb = Andeb,
                               Anfin = Anrupt[[i]],
                               alpha = alpha, plot = T )}
      #pour les autres
      if (i > 1 ){
        #plus de 5 crues dans la sous-periode ?
        if( (length(which(data <= Anrupt[[i]] & data >= Anrupt[[(i-1)]]))) < 5 ){
          print("Error ! One of the sub-periods has less than 5 floods");return()}
        #si oui on applique la fonction TestPois a la sous periode
        VecDat[[i]] = TestPois(
          data = data[which(data <= Anrupt[[i]] & data >= Anrupt[[(i-1)]])],
          Andeb = Anrupt[[(i-1)]], Anfin = Anrupt[[(i)]], alpha = alpha, plot = T)}
    }
    
    #On applique TestPois a tout l echantillon pour faire le plot de base  
    TestPois(data, Andeb = Andeb, Anfin = Anfin, alpha = alpha)
    legend(x = "topleft",legend = c("Flood event",
                                    "Startpoint/Changepoint",
                                    "Global confidence interval",
                                    "Sub-period confidence interval"),
           col = c("blue","indianred","black","green"),lty = c(1,0,2,2), pch = c(20,20,NA,NA),
           pt.cex = c(1,2,1,1), cex = 0.8)
    #On boucle sur toutes les sous periodes pour tracer les IC
    for (i in 1:length(Anrupt)){
      #pour la premiere periode
      if (i == 1){
        lines(x = VecDat[[i]]$AnExp, y = VecDat[[i]]$EllipInf, lty = 2, col = "green")
        lines(x = VecDat[[i]]$AnExp, y = VecDat[[i]]$EllipSup, lty = 2, col = "green")}
      #pour les autres  
      if (i > 1){
        Fnb = length(which(data < Anrupt[[(i-1)]]))
        lines(x = VecDat[[i]]$AnExp, y = Fnb + VecDat[[i]]$EllipInf,
              lty = 2, col = "green")
        lines(x = VecDat[[i]]$AnExp, y = Fnb + VecDat[[i]]$EllipSup,
              lty = 2, col = "green")
        points(x = VecDat[[i]]$AnneeCrue[1],
               y = Fnb, pch = 19, col = "indianred", cex = 1.5)
      }
    }
    
    
    return(VecDat)}
}


#### POUR TESTER LES FONCTIONS, DE-COMMENTER LES LIGNES SUIVANTES ET MODIFIER LE CHEMIN DU CSV
# ###Chargement donnees
# Bale = read.csv2("C://Users/mathieu.lucas/Desktop/Segmentation/TestMichel/Bale1225_2017.csv")
# 
# data = Bale[,1]
# Andeb = 1225
# Anfin = 2017
# alpha = 0.05
# 
# 
# # TEST SIMPLE
# TestPois(data = data, Andeb = Andeb, Anfin = Anfin, alpha = alpha, plot = T)
# 
# 
# # TEST AVEC 2 SOUS ECHANTILLONS, PRECISER Anrupt
# PoissSucc(data,Andeb,Anfin,Anrupt = 1764)
# 
# # TEST AVEC N SOUS ECHANTILLONS, attention, Anrupt est une liste
# PoissoN(data,Andeb,Anfin,Anrupt = list(1600,1760))


GEV_Binom = function(Qcont,   #df with an, mp, tot2.5, tot97.5 of continuous period amax discharges
                     SpagCont,#matrix of hydrometric uncertainty of continuous period amax
                     Occ,     #dates of OT occurences of historic period
                     Ts,      #list with dist, prior, init
                     Nban,    #list with dist, prior, init
                     prob,    #quantiles to extract
                     Nsim,    #number of MCMC sim
                     Nhist,   #historical period length
                     dir.res){#dir to run and write results
  
  if(Ts[[1]] == "FIX"){
    Seuil = parameter(name = "Seuil", init = Ts[[3]], priorDist = "FIX")
  } else {Seuil = parameter(name = "Seuil", init = Ts[[3]], priorDist = Ts[[1]],
                            priorPar = Ts[[2]])  }
  
  if(Nban[[1]] == "FIX"){
    NBan = parameter(name = "NBan", init = Nban[[3]], priorDist = "FIX")
  } else {
    NBan = parameter(name = "NBan", init = Nban[[3]], priorDist = Nban[[1]],
                     priorPar = Nban[[2]])  }
  
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
  alpha = 0.5
  # Durée periode continue
  NC = nrow(Qcont)
  # Durée periode histo
  NH = Nhist
  #Nombre de crues histo sup. au seuil
  NS_H = nrow(Occ)
  #Nombre de crues de la periode continue sup. au seuil
  NS_C = length(which(Qcont[,2] > Ts[[3]]))
  #Nombre tot. de crues sup. au seuil 
  NS = NS_H + NS_C
  #Nombre d'annees total
  NAn = NC + NH
  ############ Sup seuil
  SupSeuil = data.frame(rang=1:NS)
  SupSeuil$flag = "Hist"
  # On tire exactement NS_C crues "cont"
  SupSeuil$flag[sample(1:NS,NS_C,replace = F)] = "Cont"
  SupSeuil$mp = NA
  SupSeuil$tot2.5 = NA
  SupSeuil$tot97.5 = NA
  SupSeuil$mp[which(SupSeuil$flag == "Cont")] = Qcont[,2][order(Qcont[,2],decreasing = T)]
  SupSeuil$tot2.5[which(SupSeuil$flag == "Cont")] = Qcont[,3][order(Qcont[,2],decreasing = T)]
  SupSeuil$tot97.5[which(SupSeuil$flag == "Cont")] = Qcont[,4][order(Qcont[,2],decreasing = T)]
  # Freq HIRSCH 
  SupSeuil$Fr = (NS/NAn) * ((SupSeuil$rang-alpha)/(NS+1-2*alpha))
  ############ Sous seuil
  SousSeuil = Qcont[which(Qcont[,2]<=Ts[[3]]),(2:4)]
  SousSeuil = SousSeuil[order(SousSeuil[,1], decreasing = T),]
  # Freq HIRSCH
  SousSeuil$rang = (NS+1):(NS_H+NC)
  SousSeuil$flag = "Cont"
  # Freq HIRSCH
  SousSeuil$Fr = (NS/NAn) + (1-NS/NAn) * ( ((SousSeuil$rang-NS)-alpha)/(NC-NS_C+1-2*alpha) )
  ## Tab freq empirique
  Freq = merge.data.frame(SupSeuil,SousSeuil, all = T)
  Freq$tot2.5[which(Freq$flag == "Hist")] = Ts[[3]]
  Freq$tot97.5[which(Freq$flag == "Hist")] = 50000
  Freq$Pr = round(1/Freq$Fr,4)
  rhist = which(Freq$flag=="Hist")
  rcont = which(Freq$flag=="Cont")
  cptcont = 1
  for(i in 1:length(rhist)){
    while(rhist[i]-rcont[cptcont] > 0){cptcont=cptcont+1}
    Freq$tot2.5[rhist[i]] = Freq$mp[rcont[cptcont]] }
  
  ### PLOT QUANTS
  Quant = ggplot()+
    geom_ribbon(data=Quants,aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                fill="95% uncertainty interval"),alpha=0.8)+
    geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
    geom_errorbar(data=Freq[which(Freq$flag=="Hist"),],
                  aes(x=Pr,ymin = tot2.5, ymax = tot97.5),lty=2,col="lightgrey",lwd=0.5,width=0)+
    geom_point(data = Freq[which(Freq$flag=="Hist"),], aes(x=Pr,y=tot2.5),
               shape = 17,col="lightgrey")+
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
          legend.position = c(0.85,0.1))
  
  ggsave(Quant, path = dir.res, filename = paste0("Quantiles.pdf"),
         width = 12, height = 8)
}


GEV_Spag = function(Q,       #df with an, mp, tot2.5, tot97.5 of continuous period amax discharges
                    Spags,   #matrix of hydrometric uncertainty of continuous period amax
                    prob,    #quantiles to extract
                    Nsim,    #number of MCMC sim
                    dir.res){#dir to run and write results
  
  Pos = parameter(name='Pos',init = 6000)  
  Ech =  parameter(name='Ech',init = 1000) 
  Form = parameter(name='Form',init = 0.01,priorDist='Gaussian',priorPar=c(0,0.3))
  
  MegaSpag = data.frame()
  Nspag = ncol(Spags)
  for(spag in 1:Nspag){
    print(paste0("Spag = ",spag))
    dat <- dataset(Y = data.frame(Spags[,spag]))
    mod <- model(dataset=dat, parentDist ='GEV', par=list(Pos,Ech,Form))
    #### Create spag subfolder 
    dir.spag = file.path(dir.res,"Spag"); dir.create(dir.spag)
    STooDs(mod,workspace = dir.spag, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
    #### Read MCMC results for each spag
    MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
    MegaSpag = rbind(MegaSpag, MCMCres[,(1:3)])
  }
  MegaGev = matrix(nrow = nrow(MegaSpag),ncol = length(prob))
  #### Compute the quantiles for all the GeV realisations spags
  #### WARNING, HERE SHAPE HAS TO BE THE OPPOSITE THE SHAPE GIVEN BY STOODS (BY CONVENTION)
  for ( i in 1 : nrow(MegaSpag)) {MegaGev[i,] = qgev(p = prob,loc = MegaSpag$Pos[i],
                                                     scale = MegaSpag$Ech[i], 
                                                     shape = -1*MegaSpag$Form[i]) }
  PostForm = MegaSpag$Form
  #### Compute the true maxpost quantiles : maxpost of hydro sample x maxpost of GeV estim
  dat <- dataset(Y = data.frame(Q[,2]))
  mod <- model(dataset = dat, parentDist ='GEV', par = list(Pos,Ech,Form))
  #### Create MP sub-folder & run Stoods
  dir.mp = file.path(dir.res,"/Mp") ; dir.create(dir.mp,showWarnings = F)
  STooDs(mod,workspace = dir.mp, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
  MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
  Mp.GeV = MCMCres[which.max(MCMCres$post),(1:3)]
  Mp.Quant = qgev(p = prob,loc = Mp.GeV$Pos, scale = Mp.GeV$Ech,shape = -1*Mp.GeV$Form)
  quant = apply(MegaGev,MARGIN = 2,FUN = quantile, probs = c(0.025,0.975))
  Quants = data.frame(Pr=1/(1-prob), Mp=Mp.Quant, Q_2=quant[1,] ,Q_9=quant[2,])
  ### write results in res dir
  write.table(round(Quants,3), paste0(dir.res,"/Quants.txt"), row.names = F)
  write.table(data.frame(Shape = round(PostForm,3), Thres = NA, Nban = NA),
              paste0(dir.res,"/Params.txt"), row.names = F)
  write.table(round(Mp.GeV,3),  paste0(dir.res,"/Maxpost_Par.txt"), row.names = F)
  
  Freq = Q[order(Q[,2]),]
  Freq$Fr = (seq(1:length(Q[,1]))-0.5)/length(Q[,1])
  Freq$Pr = 1/(1-Freq$Fr)
  
  ### PLOT QUANTS
  Quant = ggplot()+
    geom_ribbon(data=Quants,aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                fill="95% uncertainty interval"),alpha=0.8)+
    geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
    geom_point(data = Freq, aes(x=Pr, y = mp))+
    geom_errorbar(data=Freq, aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
    scale_x_continuous(trans="log10")+
    xlab("Return period [years]")+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    scale_fill_manual(name = element_blank(),
                      values = c("#67a9cf"))+
    scale_color_manual(values = c("yellow"))+ 
    theme_bw(base_size=15)+
    coord_cartesian(xlim=c(1,max(1/(1-prob))))+
    theme(legend.title=element_blank(),
          plot.title = element_text(hjust = 0.01, vjust = -7),
          legend.position = c(0.85,0.1))  
  
  ggsave(Quant, path = dir.res, filename = paste0("Quantiles.pdf"),
         width = 12, height = 8)
  
}
