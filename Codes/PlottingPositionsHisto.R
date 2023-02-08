rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
QAll = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
                  header = T)[,c(1,2,7,8)]

Call = read.csv2(paste0(dir.data,"CX_All.csv"))
QX = QAll
CX = data.frame(an = Call$An[which(Call$An > 1500 & Call$An < 1816)],
                mp = NA, tot2.5 = 8000, tot97.5 = 20000)
Seuil = 8000
alph = 0.5

#Longueur période continue
NC = 205
#Longueur période histo
NH = 1816-1500
#Nombre de crues histo sup. au seuil
NS_H = length(CX$an)
#Nombre de crues de la période continue sup. au seuil
NS_C = length(which(QX$mp>Seuil))
#Nombre tot. de crues sup. au seuil 
NS = NS_H + NS_C
#Nombre de crues inf. au seuil
Ninf = length(which(QX$mp<=Seuil))
#Nombre d'années total
NAn = NC + NH
#proba crue continue
PC = NS_C/NS
#proba crue histo
PH = NS_H/NS

# on tire exactement NS_C crues "cont"
SupSeuil = data.frame(rang=1:NS)
SupSeuil$flag = "Hist"
SupSeuil$flag[sample(1:NS,NS_C,replace = F)]="Cont"
SupSeuil$mp = NA
SupSeuil$tot2.5 = NA
SupSeuil$tot97.5 = NA
SupSeuil$mp[which(SupSeuil$flag == "Cont")] = QX$mp[order(QX$mp,decreasing = T)]
SupSeuil$tot2.5[which(SupSeuil$flag == "Cont")] = QX$tot2.5[order(QX$mp,decreasing = T)]
SupSeuil$tot97.5[which(SupSeuil$flag == "Cont")] = QX$tot97.5[order(QX$mp,decreasing = T)]
SousSeuil = QX[which(QX$mp<Seuil),(2:4)]
SousSeuil = SousSeuil[order(SousSeuil$mp, decreasing = T),]
SousSeuil$rang = (NS+1):(NS_H+NC)
SousSeuil$flag = "Cont"
Freq = merge.data.frame(SupSeuil,SousSeuil, all = T)
Freq$tot2.5[which(Freq$flag == "Hist")] = Seuil
Freq$tot97.5[which(Freq$flag == "Hist")] = 20000
Freq$Fr = (Freq$rang-0.5)/(NS+Ninf)
Freq$Pr = 1/Freq$Fr

ggplot(Freq)+
  geom_point(aes(x=Pr,y=mp))+
  geom_errorbar(aes(x=Pr,ymin=tot2.5,ymax=tot97.5,col=flag))+
  scale_x_log10()+
  coord_cartesian(xlim=c(1,10000),ylim = c(2500,15000))+
  theme_bw()


# #ancienne m�thode
# SupSeuil$flag[sample(1:NS,NS_C,replace = F)]="Cont"
# # SupSeuil$prob = runif(NS)
# # SupSeuil$flag[which(SupSeuil$prob < PC)] = "Cont"
# SupSeuil$mp = NA
# SupSeuil$tot2.5 = NA
# SupSeuil$tot97.5 = NA
# SupSeuil$flag[which(SupSeuil$rang > which(SupSeuil$flag == "Cont")[NS_C])]="Hist"


