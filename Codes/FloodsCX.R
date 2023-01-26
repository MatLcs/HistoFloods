# source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

#### Data loading
CXall = read.csv2(paste0(dir.data,"CX_All.csv"))
# Q = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
#   header = T)[,c(1,2,7,8)]
#### Homogeneity of systematic records
# plot(x=Q$an, y = Q$mp,type="l", ylab = "Discharge",xlab="Date")
# mktrend = mk.test(Q$mp)
# ptdate = pettitt.test(Q$mp)
# mktrend$p.value<0.05
# ptdate$p.value<0.05
# #### Homogeneity of flood occurences
# Andeb = 1500
# Anfin = 1816

CX = CXall[which(CXall$An > Andeb & CXall$An < Anfin),]
#### CALCULATION OF START DATE FROM ILARIA PRODOSCIMI
## number of floods after 1500
k = nrow(CX[which(CX$An > Andeb & CX$An < Anfin),])
## differences between time of occurrence of historical events and start of systematic records
tk = Anfin - CX$An[which(CX$An > Andeb)][1]
NbAnIlaria = tk + tk/k - 1

###Poisson process
# PCX=PoissoN(CX$An[which(CX$An>Andeb)],Andeb=Andeb,Anfin = Anfin, alpha = 0.025)
# PC3=PoissoN(CX$An[which(CX$An>Andeb & CX$Cat=="C4")],Andeb=Andeb,Anfin = Anfin)
# PC4=PoissoN(CX$An[which(CX$An>Andeb & CX$Cat=="C3")],Andeb=Andeb,Anfin = Anfin)
# pdf(file = paste0(dir.plots,"Poisson_C3-C4.pdf"), width = 12, height = 7)
# plot(x = PCX$AnneeCrue, y = PCX$FloodNumber, type='b', pch = 19, col ="purple",
#      ylab = "Cumulated flood number", xlab = "Years")
# #C3 only
# points(x = PC3$AnneeCrue, y = PC3$FloodNumber, type = 'b', pch = 19, col = "blue")
# #C4 only
# points(x = PC4$AnneeCrue, y = PC4$FloodNumber, type = 'b', pch = 19, col = "red")
# #confidence intervals
# lines(x = PCX$AnExp, y =  PCX$EllipInf, lty = 2,col = "purple")
# lines(x = PCX$AnExp, y =  PCX$EllipSup, lty = 2,col="purple")
# lines(x = PC3$AnExp, y =  PC3$EllipInf, lty = 2,col = "blue")
# lines(x = PC3$AnExp, y =  PC3$EllipSup, lty = 2,col="blue")
# lines(x = PC4$AnExp, y =  PC4$EllipInf, lty = 2,col = "red")
# lines(x = PC4$AnExp, y =  PC4$EllipSup, lty = 2,col="red")
# #starting points
# points(x = PCX$AnneeCrue[1], y = PCX$FloodNumber[1], pch = 19,col = "indianred", cex = 1.5)
# points(x = tail(PCX$AnneeCrue,1),y = tail(PCX$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
# points(x = tail(PC3$AnneeCrue,1),y = tail(PC3$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
# points(x = tail(PC4$AnneeCrue,1),y = tail(PC4$FloodNumber,1),pch = 19, col = "indianred", cex = 1.5)
# #legend
# legend(x = "topleft",legend = c("C3 & C4 floods",
#                                 "C3 floods",
#                                 "C4 floods",
#                                 "Start/End ",
#                                 "95% confidence intervals"),
#        col = c("purple","blue","red","indianred","black"),lty = c(1,1,1,0,2),
#        pch = c(20,20,20,20,NA),   pt.cex = c(1,1,1,2,1), cex = 1.2)
# dev.off()

#### Perception threshold determination
Q = Q[which(Q$an < 2000),]
Q$type = "none"
## associate floods from systematic period to pichard flood type
for(y in 1:length(CXall$An)) { Q$type[which(Q$an == CXall$An[y])] = CXall$Cat[y] }
## minimum discharge of C3 flood
# Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],]
TsC3 = Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],]
TsC4 =  Q[which(Q$type == "C4")[(which.min(Q$mp[which(Q$type == "C4")]))],]

write.table(TsC3, paste0(dir.data,"ThresholdC3.txt"),row.names = F)
write.table(TsC4,  paste0(dir.data,"ThresholdC4.txt"),row.names = F)

# Categ = 
ggplot(data=Q)+
  geom_errorbar(aes(x = an, ymin = tot2.5,ymax=tot97.5,color=type))+
  geom_point(aes(x = an, y=mp ,color=type))+
  geom_line(aes(x = an, y = rep(min(mp[which(type=="C3")]),length(an))),
            col = "blue",alpha = 0.6,lty=2,lwd=1)+
  geom_line(aes(x = an, y = rep(min(mp[which(type=="C4")]),length(an))),
            col="red",alpha = 0.6,lty=2,lwd=1)+
  geom_line(aes(x = an, y = rep(TsC3$tot2.5,length(an))),
            col="black",alpha = 0.6,lty=2,lwd=1)+
  scale_color_manual("Flood type",values=c("blue","red","grey"))+
  theme_light()+
  labs(x="Years",y="Discharge [m3/s]")+
  coord_cartesian(xlim = c(1817,1998))

ggsave(path = dir.plots, filename = "C3-C4_SystematicPeriod.pdf", width = 12, height = 8)

# ### Pettit and MK tests
# Pett = rep(NA,Nspag); Mkt = Pett
# for(spag in 1:Nspag){
#   Pett[spag] =  pettitt.test(Spags[,spag])[4][1]
#   Mkt[spag] = mk.test(Spags[,spag])[2][1] }
# length(which(Pett>0.05))/Nspag
# length(which(Mkt>0.05))/Nspag


