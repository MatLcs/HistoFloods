rm(list=ls())
require(ggplot2)
require(forcats)
require(ggpubr)

dirgit = "C://Users/mathieu.lucas/Desktop/GitMat/"

YC4 = c(1827,1840,1843,1856,1872,1882,1886,1907,1910,1935,1951,1955,1960,1993,1994)
YC3 = c(1841,1842,1889,1890,1896,1900,1917,1919,1923,1924,1926,1928,1933,1937,1954,1958,1963,1976)

Q = read.table(paste0(dirgit,"PropagMaxAn/Results/Quantiles_Amax.txt"),header = T)
Q = Q[which(Q$an < 2000),c(1,2,7,8)]
Q$type="none"

for(y in 1:length(YC4)) { Q$type[which(Q$an == YC4[y])] = "C4" }
for(y in 1:length(YC3)) { Q$type[which(Q$an == YC3[y])] = "C3" }

Categ = ggplot(data=Q)+
  geom_errorbar(aes(x = an, ymin = tot2.5,ymax=tot97.5,color=type))+
  geom_point(aes(x = an, y=mp ,color=type))+
  geom_line(aes(x = an, y = rep(min(mp[which(type=="C3")]),184)),
            col="blue",alpha = 0.6,lty=2,lwd=1)+
  geom_line(aes(x = an, y = rep(min(mp[which(type=="C4")]),184)),
            col="red",alpha = 0.6,lty=2,lwd=1)+
  scale_color_manual("Flood type",values=c("blue","red","grey"))+
  theme_light()+
  labs(x="Years",y="Discharge [m3/s]")+
  coord_cartesian(xlim = c(1817,1998))

Q$type = as.factor(Q$type)
Q$type = fct_relevel(Q$type,c("none","C3","C4"))

Dens = ggplot(data=Q)+
  geom_density(aes(x=mp,color=type),lwd=1,position = "identity")+
  geom_vline(aes(xintercept=min(mp[which(type=="C3")])),col="blue",alpha=0.6,lwd=1,lty=2)+
  geom_vline(aes(xintercept=min(mp[which(type=="C4")])),col="red",alpha=0.6,lwd=1,lty=2)+
  xlab("Discharge [m3/s]")+
  ylab("Density")+
  xlim(2000,14000)+
  theme_light()+
  scale_color_manual("Flood type",values=c("grey","blue","red"))


ggarrange(Categ, Dens, nrow = 2, ncol = 1, align = "hv",common.legend = T,legend = "right")
# ggsave("FloodTypes.pdf",path = paste0(dirgit,"/CruesHisto"),
#       width = 10,height = 8)
  

# 
# plot(density(Q$mp[which(Q$type=="none")]), col="grey", lwd = 2,
#      ylim = c(0,5.5E-4), xlim=c(min(density(Q$mp)$x),max(density(Q$mp)$x)),
#      main = "AMAX discharges distribution @Beaucaire (1816-2000)", xlab = "Discharge [m3/s]")
#   lines(density(Q$mp[which(Q$type=="C3")]),col="blue", lwd=2)
#   lines(density(Q$mp[which(Q$type=="C4")]),col="red",lwd=2)
#   legend("topright",legend = c("none","C3","C4"), lty = c(1,1,1),
#          col = c("grey","blue","red"), title = "Category", cex=0.8)
#   abline(v=min(Q$mp[which(Q$type=="none")]),col="grey", lty = 2)
#   abline(v=min(Q$mp[which(Q$type=="C3")]),col="blue", lty = 2)
#   abline(v=min(Q$mp[which(Q$type=="C4")]),col="red", lty = 2)
