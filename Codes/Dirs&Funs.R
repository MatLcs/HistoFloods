require(lubridate);require(ggplot2);require(trend);require(ggpubr); require(forcats)
require(cowplot);require(RSTooDs);require(evd);require(reshape);require(ggridges);
require(RColorBrewer)

dir.all = "C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/"
dir.data = paste0(dir.all,"Data/") 
dir.codes = paste0(dir.all,"Codes/") 
dir.plots = paste0(dir.all,"Plots/") 
dir.res = paste0(dir.all,"Results/")

#### PARAMS SIMULATION
Andeb = 1500
Anfin = 2000
Nsim = 5000
Nspag = 500
prob = seq(0.01,0.999,0.001) ; Pr = 1/(1-prob)

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

