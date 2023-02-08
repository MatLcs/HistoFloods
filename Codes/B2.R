# rm(list=ls())
# source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

#XXXX-2020 / unknown threshold
case = "B2"
dir.case = paste0(dir.res.ts,case,"_case/"); dir.create(dir.case, showWarnings = F)

# Thresh = read.table(paste0(dir.data,"Threshold",caseTs,".txt"), header = T)
# Ts = Thresh$mp
# sd.Ts = (Thresh$tot97.5 - Thresh$mp)/2

# sd.NbAn = 50

# Ts = 6740
# sd.Ts = 500
# Andeb = 1500 
# Anfin = 2000
# Nspag = 100
# Nsim = 2000

#### Data loading
# Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
# CX = read.csv2(paste0(dir.data,"CX_All.csv"))
# Q = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
#                header = T)[,c(1,2,7,8)]
#### SToods parameters, common GEV priors to all cases
# Pos = parameter(name='Pos',init = 6000) 
# Ech =  parameter(name='Ech',init = 1000) 
# Form = parameter(name='Form',init = 0.01,priorDist='Gaussian',priorPar=c(0,0.3))
Seuil = parameter(name = "Seuil", init = Ts, priorDist = "Gaussian", priorPar = c(Ts,sd.Ts))
NbAn = parameter(name = "NbAn", init = Anfin-Andeb, priorDist = "Uniform", 
                 priorPar = distNbAn)
#### Quantiles extracted up to Q1000 & return period associated
# prob = seq(0.01,0.999,0.001) ; Pr = 1/(1-prob)
# 
# CX = CX[which(CX$An > Andeb),]
# if(caseTs == "C4"){
#   CX = CX[which(CX$Cat == caseTs),] }
# 
# CX = Call[which(Call$An > Andeb),]
# if(caseTs == "C4"){
#   CX = CX[which(CX$Cat == caseTs),] }

#### Initializing DF for MCMC results
MegaSpag = data.frame()

for(spag in 1:Nspag){
  # spag = 1
  print(paste0("Spag = ",spag))
  #draw a random spag within the 500 Q AMAX spags
  Y = data.frame(  c(Spags[,spag],#[,sample(1:500,1)],
                    length(CX$An)) )
  Var = data.frame(factor(c(rep('Q',length(Q$mp)),"Occ")))
  dat <- dataset(Y = Y, var = Var)
  mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
               par=list(Pos,Ech,Form,Seuil,NbAn),
               ### PARAM GEV
               formula = c("Pos=Pos","Ech=Ech","Form = Form",
                           ### PARAM BINOMIAL
                           "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                           "NbAn = NbAn"))
  #### Create spag subfolder 
  dir.spag = file.path(dir.case,"Spag"); dir.create(dir.spag,showWarnings = F)
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
PostNban = MegaSpag$NbAn
#### Compute the true maxpost quantiles : maxpost of hydro sample x maxpost of GeV estim
Y = data.frame(c(Q$mp, length(CX$An)) )
Var = data.frame(factor(c(rep('Q',length(Q$mp)),"Occ")))
dat <- dataset(Y = Y, var = Var)
#### Create MP sub-folder & run Stoods
dir.mp = file.path(dir.case,"/Mp") ; dir.create(dir.mp,showWarnings = F)
dat <- dataset(Y = Y, var = Var)
mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
             par=list(Pos,Ech,Form,Seuil,NbAn),
             ### PARAM GEV
             formula = c("Pos=Pos","Ech=Ech","Form = Form",
                         ### PARAM BINOMIAL
                         "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                         "NbAn = NbAn"))
STooDs(mod,workspace = dir.mp, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
Mp.GeV = MCMCres[which.max(MCMCres$post),(1:5)]
Mp.Quant = qgev(p = prob,loc = Mp.GeV$Pos, scale = Mp.GeV$Ech,shape = -1*Mp.GeV$Form)

quant = apply(MegaGev,MARGIN = 2,FUN = quantile, probs = c(0.025,0.975))
Quants = data.frame(Pr=Pr, Mp=Mp.Quant, Q_2=quant[1,] ,Q_9=quant[2,])

write.table(round(Quants,3), paste0(dir.case,case,"_Quants.txt"), row.names = F)
write.table(round(data.frame(Shape = PostForm, Thres = PostThres, Nban = PostNban),3),
            paste0(dir.case,case,"_Params.txt"), row.names = F)
write.table(Mp.GeV,  paste0(dir.case,case,"_Maxpost_Par.txt"), row.names = F)


#### plot positions histo
NC = length(Q$an)
NH = Anfin-Andeb
#Nombre de crues histo sup. au seuil
NS_H = length(CX$An)
#Nombre de crues de la période continue sup. au seuil
NS_C = length(which(Q$mp>Ts))
#Nombre tot. de crues sup. au seuil 
NS = NS_H + NS_C
#Nombre de crues inf. au seuil
Ninf = length(which(Q$mp <= Ts))
#Nombre d'années total
NAn = NC + NH
#proba crue continue
PC = NS_C/NS
#proba crue histo
PH = NS_H/NS
#data frame crues supseuil
SupSeuil = data.frame(rang=1:NS)
SupSeuil$flag = "Hist"
SupSeuil$flag[sample(1:NS,NS_C,replace = F)]="Cont"
SupSeuil$mp = NA
SupSeuil$tot2.5 = NA
SupSeuil$tot97.5 = NA
SupSeuil$mp[which(SupSeuil$flag == "Cont")] = Q$mp[order(Q$mp,decreasing = T)]
SupSeuil$tot2.5[which(SupSeuil$flag == "Cont")] = Q$tot2.5[order(Q$mp,decreasing = T)]
SupSeuil$tot97.5[which(SupSeuil$flag == "Cont")] = Q$tot97.5[order(Q$mp,decreasing = T)]
SousSeuil = Q[which(Q$mp<Ts),(2:4)]
SousSeuil = SousSeuil[order(SousSeuil$mp, decreasing = T),]
SousSeuil$rang = (NS+1):(NS_H+NC)
SousSeuil$flag = "Cont"
Freq = merge.data.frame(SupSeuil,SousSeuil, all = T)
Freq$tot2.5[which(Freq$flag == "Hist")] = Ts
Freq$tot97.5[which(Freq$flag == "Hist")] = 50000
Freq$Fr = (Freq$rang-0.5)/(NS+Ninf)
Freq$Pr = 1/Freq$Fr


### PLOT QUANTS
Quant = ggplot()+
  geom_ribbon(data=Quants,aes(x=Pr, ymin = Q_2, ymax=Q_9,
                              fill="95% uncertainty interval"),alpha=0.8)+
  geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
  geom_errorbar(data=Freq[which(Freq$flag=="Cont"),],
                aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
  geom_point(data = Freq[which(Freq$flag=="Cont"),],
             aes(x=Pr, y = mp))+
  geom_errorbar(data=Freq[which(Freq$flag=="Hist"),],
                aes(x=Pr,ymin = tot2.5, ymax = tot97.5),lty=2)+
  geom_point(data = Freq[which(Freq$flag=="Hist"),], aes(x=Pr,y=tot2.5),shape = 17)+
  
  scale_x_continuous(trans="log10")+
  xlab("Return period [years]")+
  ylab("Discharge [m3/s]")+
  scale_fill_manual(name = element_blank(),
                    values = c("#67a9cf"))+
  scale_color_manual(values = c("yellow"))+ #"royalblue")+
  theme_bw(base_size=15)+
  # labs(title = paste0(head(Q$an,1)," - ",tail(Q.case$an,1)))+
  coord_cartesian(xlim=c(1,max(Pr)),ylim = c(min(Quants$Q_2),20000))+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.01, vjust = -7),
        legend.position = c(0.8,0.2))

ggsave(Quant, path = dir.case, filename = paste0("Quantiles_",case,caseTs,".pdf"),
       width = 10, height = 7)
# Quant = ggplot()+
#   geom_ribbon(data=Quants,aes(x=Pr, ymin = Q_2, ymax=Q_9,
#                               fill="95% uncertainty interval"),alpha=0.8)+
#   geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
#   geom_point(data = Freq, aes(x=Pr, y = mp))+
#   geom_errorbar(data=Freq, aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
#   scale_x_continuous(trans="log10")+
#   xlab("Return period [years]")+
#   ylab("Discharge [m3/s]")+
#   scale_fill_manual(name = element_blank(),
#                     values = c("#67a9cf"))+
#   scale_color_manual(values = c("yellow"))+ #"royalblue")+
#   theme_bw(base_size=15)+
#   # labs(title = paste0(head(Q$an,1)," - ",tail(Q.case$an,1)))+
#   coord_cartesian(xlim=c(1,max(Pr)))+
#   theme(legend.title=element_blank(),
#         plot.title = element_text(hjust = 0.01, vjust = -7),
#         legend.position = c(0.8,0.2))  
# 
# ggsave(Quant, path = dir.case, filename = paste0("Quantiles_",case,caseTs,".pdf"),width = 10, height = 7)
# 
# 
