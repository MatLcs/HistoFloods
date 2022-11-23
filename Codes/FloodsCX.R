# source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

#### Data loading
CX = read.csv2(paste0(dir.data,"CX_bcr.csv"))
Q = read.table("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]

# #### Homogeneity of flood occurences
# Andeb = 1500
# Anfin = 2000

#### CALCULATION OF START DATE FROM ILARIA PRODOSCIMI
## number of floods after 1500
k = nrow(CX[which(CX$An > Andeb),])
## differences between time of occurrence of historical events and start of systematic records
tk = 2000 - CX$An[which(CX$An > Andeb)][1]
NbAnIlaria = tk + tk/k - 1

pdf(file = paste0(dir.plots,"Poisson_C3+C4.pdf"), width = 12, height = 7)
### WITH START DATE FROM ILARIA
PoissoN(CX$An[which(CX$An>Andeb)],Andeb=1507,Anfin = Anfin, alpha = 0.025)
dev.off()
# PoissoN(CX$An[which(CX$An>Andeb & CX$Cat=="C4")],Andeb=Andeb,Anfin = Anfin)
# PoissoN(CX$An[which(CX$An>Andeb & CX$Cat=="C3")],Andeb=Andeb,Anfin = Anfin)

#### Perception threshold determination
Q = Q[which(Q$an < 2000),]
Q$type = "none"
for(y in 1:length(CX$An)) { Q$type[which(Q$an == CX$An[y])] = CX$Cat[y] }

Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],]
Ts = Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],3]

write.table(Q[which(Q$type != "none")[(which.min(Q$mp[which(Q$type != "none")]))],],
            paste0(dir.data,"Threshold.txt"),row.names = F)

# Categ = 
ggplot(data=Q)+
  geom_errorbar(aes(x = an, ymin = tot2.5,ymax=tot97.5,color=type))+
  geom_point(aes(x = an, y=mp ,color=type))+
  geom_line(aes(x = an, y = rep(min(mp[which(type=="C3")]),length(an))),
            col = "blue",alpha = 0.6,lty=2,lwd=1)+
  geom_line(aes(x = an, y = rep(min(mp[which(type=="C4")]),length(an))),
            col="red",alpha = 0.6,lty=2,lwd=1)+
  geom_line(aes(x = an, y = rep(Ts,length(an))),
            col="black",alpha = 0.6,lty=2,lwd=1)+
  scale_color_manual("Flood type",values=c("blue","red","grey"))+
  theme_light()+
  labs(x="Years",y="Discharge [m3/s]")+
  coord_cartesian(xlim = c(1817,1998))

ggsave(path = dir.plots, filename = "C3-C4_SystematicPeriod.pdf", width = 12, height = 8)



