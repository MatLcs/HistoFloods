source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

#### Data loading
CXall = read.csv2(paste0(dir.data,"CX_All.csv"))
Q = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
#### Homogeneity of systematic records
# plot(x=Q$an, y = Q$mp,type="l", ylab = "Discharge",xlab="Date")
# mktrend = mk.test(Q$mp)
# ptdate = pettitt.test(Q$mp)
# mktrend$p.value<0.05
# ptdate$p.value<0.05
# #### Homogeneity of flood occurences
Andeb = 1500
Anfin = 1816

CX = CXall[which(CXall$An > Andeb & CXall$An < Anfin),]
#### CALCULATION OF START DATE FROM ILARIA PRODOSCIMI
## number of floods after 1500
k = nrow(CX[which(CX$An > Andeb & CX$An < Anfin),])
## differences between time of occurrence of historical events and start of systematic records
tk = Anfin - CX$An[which(CX$An > Andeb)][1]
NbAnIlaria = tk + tk/k - 1

###Poisson process
PCX=PoissoN(CX$An[which(CX$An>Andeb)],Andeb=Andeb,Anfin = Anfin)
# PC3=PoissoN(CX$An[which(CX$An>Andeb & CX$Cat=="C4")],Andeb=Andeb,Anfin = Anfin)
PC4=PoissoN(CX$An[which(CX$An>Andeb & CX$Cat=="C3")],Andeb=Andeb,Anfin = Anfin)
# # PoissoN(CXall$An[which(CXall$An > Andeb)],Andeb=Andeb,Anfin = 2000,Anrupt = list(Anfin))
# pdf(file = paste0(dir.plots,"Poisson_C3-C4_FR.pdf"), width = 12, height = 7)
plot(x = PCX$AnneeCrue, y = PCX$FloodNumber, type='b', pch = 19, col ="blue",
     # ylab = "Cumulated flood number", xlab = "Years")
     ylab = "Nombre cumulé de crues", xlab = "Années")
#C3 only
# points(x = PC3$AnneeCrue, y = PC3$FloodNumber, type = 'b', pch = 19, col = "blue")
#C4 only
points(x = PC4$AnneeCrue, y = PC4$FloodNumber, type = 'b', pch = 19, col = "red")
#confidence intervals
lines(x = PCX$AnExp, y =  PCX$EllipInf, lty = 2,col = "blue")
lines(x = PCX$AnExp, y =  PCX$EllipSup, lty = 2,col="blue")
# lines(x = PC3$AnExp, y =  PC3$EllipInf, lty = 2,col = "blue")
# lines(x = PC3$AnExp, y =  PC3$EllipSup, lty = 2,col="blue")
lines(x = PC4$AnExp, y =  PC4$EllipInf, lty = 2,col = "red")
lines(x = PC4$AnExp, y =  PC4$EllipSup, lty = 2,col="red")
#starting points
points(x = PCX$AnneeCrue[1], y = PCX$FloodNumber[1], pch = 19,col = "indianred", cex = 1.5)
points(x = tail(PCX$AnneeCrue,1),y = tail(PCX$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
# points(x = tail(PC3$AnneeCrue,1),y = tail(PC3$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
points(x = tail(PC4$AnneeCrue,1),y = tail(PC4$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
#legend EN
# legend(x = "topleft",legend = c("C3 & C4 floods",
#                                 # "C3 floods",
#                                 "C4 floods",
#                                 "Start/End ",
#                                 "95% confidence intervals"),
# 
#        col = c(
#          # "purple",
#                  "blue","red","indianred","black"),lty = c(1,#1,
#                        1,0,2),
#        pch = c(20,#20,
#                20,20,NA),   pt.cex = c(1,#1,
#                                        1,2,1), cex = 1.2)
#Legend FR
legend(x = "topleft",legend = c("Catég. C3&C4",
                                "Catég. C4",
                                "Début/Fin",
                                "Intervalle de conf. 95%"),

       col = c( "blue","red","indianred","black"),lty = c(1,1,0,2),
       pch = c(20, 20,20,NA),   pt.cex = c(1,1,2,1), cex = 1.2)
# dev.off()


P7 = PoissoN(Q$an[which(Q$mp>7000 & Q$an < 1970)],Andeb=1816,Anfin = 1970)
P9 = PoissoN(Q$an[which(Q$mp>9000 & Q$an < 1970)],Andeb=1816,Anfin = 1970)

P8_b = PoissoN(Q$an[which(Q$mp>8000)],Andeb=1816,Anfin = 2020, Anrupt = list(1940))

# pdf(file = paste0(dir.plots,"Poisson_Qrecent_FR.pdf"), width = 12, height = 7)
  plot(x = P7$AnneeCrue, y = P7$FloodNumber, type='b', pch = 19, col ="blue",
       ylab = "Nombre cumulé de crues", xlab = "Années")
  points(x = P9$AnneeCrue, y = P9$FloodNumber, type = 'b', pch = 19, col = "red")
  #confidence intervals
  lines(x = P7$AnExp, y =  P7$EllipInf, lty = 2,col = "blue")
  lines(x = P7$AnExp, y =  P7$EllipSup, lty = 2,col="blue")
  lines(x = P9$AnExp, y =  P9$EllipInf, lty = 2,col = "red")
  lines(x = P9$AnExp, y =  P9$EllipSup, lty = 2,col="red")
  #starting/end points
  points(x = P7$AnneeCrue[1], y = P7$FloodNumber[1], pch = 19,col = "indianred", cex = 1.5)
  points(x = tail(P7$AnneeCrue,1),y = tail(P7$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
  points(x = tail(P9$AnneeCrue,1),y = tail(P9$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
  #legend expression(paste("Discharge [",m^3,".",s^-1,"]",sep=""))
  legend(x = "topleft",legend = c(expression(paste("Seuil de perc. 7000 ",m^3,"/",s,sep = "")),
                                  expression(paste("Seuil de perc. 9000 ",m^3,"/",s,sep = "")),
                                  "Début/Fin",
                                  "Intervalle de conf. 95%"),
         col = c("blue","red","indianred","black"),lty = c(1,1,0,2),pch = c(20,20,20,NA),
                 pt.cex = c(1,1,2,1), cex = 1.2)
# dev.off()

#### Perception threshold determination
# Q = Q[which(Q$an < 2000),]
# Q$type = "none"
## associate floods from systematic period to pichard flood type
# for(y in 1:length(CXall$An)) { Q$type[which(Q$an == CXall$An[y])] = CXall$Cat[y] }
## minimum discharge of C3 flood
# Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],]
# TsC3 = Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],]
# TsC4 =  Q[which(Q$type == "C4")[(which.min(Q$mp[which(Q$type == "C4")]))],]

# write.table(TsC3, paste0(dir.data,"ThresholdC3.txt"),row.names = F)
# write.table(TsC4,  paste0(dir.data,"ThresholdC4.txt"),row.names = F)

# Categ = 
# ggplot(data=Q)+
#   geom_errorbar(aes(x = an, ymin = tot2.5,ymax=tot97.5,color=type))+
#   geom_point(aes(x = an, y=mp ,color=type))+
#   geom_line(aes(x = an, y = rep(min(mp[which(type=="C3")]),length(an))),
#             col = "blue",alpha = 0.6,lty=2,lwd=1)+
#   geom_line(aes(x = an, y = rep(min(mp[which(type=="C4")]),length(an))),
#             col="red",alpha = 0.6,lty=2,lwd=1)+
#   geom_line(aes(x = an, y = rep(TsC3$tot2.5,length(an))),
#             col="black",alpha = 0.6,lty=2,lwd=1)+
#   scale_color_manual("Flood type",values=c("blue","red","grey"))+
#   theme_light()+
#   labs(x="Years",y="Discharge [m3/s]")+
#   coord_cartesian(xlim = c(1817,1998))

# ggsave(path = dir.plots, filename = "C3-C4_SystematicPeriod.pdf", width = 12, height = 8)

# ### Pettit and MK tests
# Pett = rep(NA,Nspag); Mkt = Pett
# for(spag in 1:Nspag){
#   Pett[spag] =  pettitt.test(Spags[,spag])[4][1]
#   Mkt[spag] = mk.test(Spags[,spag])[2][1] }
# length(which(Pett>0.05))/Nspag
# length(which(Mkt>0.05))/Nspag

  
  
### PLOT ECHANTILLON MIXTE
SpagsGG = melt(Spags)
GGcrues = ggplot()+
  geom_errorbar(aes(x = Q$an, ymin = Q$tot2.5,ymax=Q$tot97.5,color="Période\ncontinue"),width = 2)+
  geom_point(aes(x = Q$an, y=Q$mp,color="Période\ncontinue" ))+
  #C3
  geom_polygon(aes(x = c(1500,1816,1816,1500),y = c(5000,5000,9000,9000),fill="S3"),alpha=0.1)+
  geom_errorbar(aes(x=CX$An[which(CX$Cat == "C3")], ymin = 7000, ymax = 20000,color="Cat. C3"),
                lty=2,width=0)+
  geom_point(aes(x=CX$An[which(CX$Cat == "C3")], y = 7000,color="Cat. C3"),
             shape = 17, size = 2.5)+
  #C4
  geom_polygon(aes(x = c(1500,1816,1816,1500),y = c(7000,7000,11000,11000),fill="S4"),alpha=0.1)+
  geom_errorbar(aes(x=CX$An[which(CX$Cat == "C4")], ymin = 9000, ymax = 20000,color="Cat. C4"),
                lty=2,width=0)+
  geom_point(aes(x=CX$An[which(CX$Cat == "C4")], y = 9000,color="Cat. C4"),
             shape = 17,size=2.5)+
  scale_fill_manual(values=c("blue","red"))+
  scale_color_manual(values=c("blue","red","black"))+
  coord_cartesian(xlim = c(1520,2000),ylim=c(2400,15000))+
  theme_light(base_size = 25)+
  xlab("Années")+
  ylab(expression(paste("Débit [",m^3,"/",s,"]",sep = "")))+
  labs(fill = "", colour = "")
  
# ggsave(plot = GGcrues, filename = "EchMixteBcr.pdf",path = dir.plots, width = 15, height = 10)



