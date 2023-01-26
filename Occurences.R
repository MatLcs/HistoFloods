rm(list=ls())
require(RSTooDs)

dirgit = "C://Users/mathieu.lucas/Desktop/GitMat/"
dir.res = paste0(dirgit,"CruesHisto/Results/");dir.create(dir.res)
dir.case = paste0(dir.res,"OccC3_C4_NbAnFIX");dir.create(dir.case)

C3 = rev(c(
       1811,1808,1802,1801,1791,1790,1776,1774,1765,1763,1760,1758,1756,1747,1726,1725,
       1711,1710,1709,1708,1702,1697,1689,1685,1679,1660,1657,1651,1647,1637,1636,1633,
       1629,1616,1593,1591,1587,1586,1583,1578,1571,1566,1561,1551,1544,1542,1541,1530,
       1529,1525,1522,1515,1511,1507,1500 ))
C4 = rev(c( 1810,1801,1755,1745,1711,1706,1705,1694,1674,1580,1573,1570,1548,1529 ))
YC4 = c(1827,1840,1843,1856,1872,1882,1886,1907,1910,1935,1951,1955,1960,1993,1994)
YC3 = c(1841,1842,1889,1890,1896,1900,1917,1919,1923,1924,1926,1928,1933,1937,1954,1958,1963,1976)
AnOcc = 1816 - 1500
NoccC3 = length(C3)
NoccC4 = length(C4)

# write.table(c(C4,YC4),file = paste0(dirgit,"CruesHisto/Data/C4.txt"),row.names = F,col.names = F)
# write.table(c(C3,YC3),file = paste0(dirgit,"CruesHisto/Data/C3.txt"),row.names = F,col.names = F)

Q = read.table(paste0(dirgit,"PropagMaxAn/Results/Quantiles_Amax.txt"),header = T)[,c(1,2,7,8)]



Both = c(C4,YC4,C3,YC3)[order(c(C4,YC4,C3,YC3))]
plot(x=Both, y = seq_along(Both),ylab="Cumulative frequency", xlab = "Year")
fitboth=lm((seq_along(Both[which(Both<1816)]))~Both[which(Both<1816)])
abline(fitboth)
points(x=c(C3,YC3),y=seq_along(c(C3,YC3)),col="blue")
fitC3=lm((seq_along(C3))~C3)
abline(fitC3,col="blue")
points(x=c(C4,YC4),y=seq_along(c(C4,YC4)),col="red")
fitC4=lm((seq_along(C4))~C4)
abline(fitC4,col="red")
abline(v = 1816, col = "grey")
legend("topleft",legend = c("C3+C4","C3","C4"), col = c("black","blue","red"),pch = c(1,1,1))


Spags = read.table(paste0(dirgit,"PropagMaxAn/Results/Spags_uTot_Amax.txt"))
Nsim = 5000

Y = data.frame(c(Q$mp,NoccC4+NoccC4))
Var = data.frame(factor(c(rep('Q',length(Q$mp)),"Occ")))

# prend uniquement un df
dat <- dataset(Y=Y,var = Var)

# param de forme de la GEV, centre sur 0 qui revient a Gumbel
Form = parameter(name='Form',init = 0.01,priorDist='Gaussian',priorPar=c(0,0.2))
Ech = parameter(name='Ech',init = 1000)
Pos = parameter(name='Pos',init = 6000)
Seuil = parameter(name = "Seuil", init = 7000, priorDist = "Gaussian", priorPar = c(7000,1000))
# NbAn = parameter(name = "NbAn", init = AnOcc, priorDist = "Gaussian",priorPar = c(AnOcc,100))
NbAn = parameter(name = "NbAn", init = AnOcc, priorDist = "FIX")
# definition du modele
mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
             par=list(Pos,Ech,Form,Seuil,NbAn),
             #l ordre dans la formule est important
                         ### PARAM GEV
             formula = c("Pos=Pos","Ech=Ech","Form = Form",
                         ### PARAM BINOMIAL
                         "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                         "NbAn = NbAn"))

# creation dossier pour stocker les resultats
# dir.create(path = dirSt)
STooDs(mod,workspace = dir.case,mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))

#Report file pdf pour regarder la trace integrale
repMc = readMCMC(file = file.path(dir.case,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2,
                 reportFile = file.path(dir.case,"TraceAll1.pdf"))

setwd(dir.case)
MCMCreport(mcmc=repMc,model = mod)

repMc[which.max(repMc$post),]


