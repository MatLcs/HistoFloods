rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

###### RUN PARAMS
Nsim = 8000; Nspag = 500 ; textsize = 15
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
starthistartif = 1816
endhistartif = 1970
Thres = 9000
Nbhists = endhistartif - starthistartif
Run = T

###### DATA LOADING
###Histrhone data 
Call = read.csv2(paste0(dir.data,"CX_All.csv"))
### 1816-2020 AMAX discharges with uncertainty and spags
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
Q = read.table(
  "C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]

mask = which(Q$mp > Thres)
maskhist = which(Q$mp > Thres & Q$an < endhistartif)

Qconts = rbind(Q[maskhist,],Q[which(Q$an >= endhistartif),])
Spagsss = rbind(Spags[maskhist,],Spags[which(Q$an >= endhistartif),])
Occs = data.frame(an=NULL)

dir.case = paste0(dir.res,"Censure"); dir.create(dir.case,showWarnings = F)
dir.plot = paste0(dir.plots,"Censure"); dir.create(dir.plot, showWarnings = F)

if(Run == T){
  # Fun run
  GEV_Binom(Qcont = Qconts,
            SpagCont = Spagsss,
            Occ = Occs,
            Ts = list("FIX",c(Thres,1),Thres),#dist,prior,init
            Nban = list("FIX",
                        c(Nbhists - length(maskhist),1),
                        Nbhists - length(maskhist)),#dist,prior,init
            prob = prob,
            Nsim = Nsim,
            Nhist = Nbhists,
            dir.res = dir.case)
}

##### COMPARE WITH BASELINE & A
Censure.dir = paste0(dir.res,"Censure")
pal.Censure = "#7b3294"
Base.dir = paste0(dir.res,"Baseline_1816-2020")
pal.Base = "#4dac26"
A.dir = paste0(dir.res,"Artif2/model_A")
pal.A = brewer.pal(4,"RdYlBu")[c(2,1,3,4)][1]

models = c("Baseline_1816-2020","Modèle_A","Censure")
dirs = c(Base.dir,A.dir,Censure.dir)
pals = c(pal.Base,pal.A,pal.Censure) 

ind100 = length(Pr) - 18
ind1000 = length(Pr) - 9

Quant100 = list(); Quant1000 = list()

for(m in 1:length(models)){
  ##Q100 & Q1000
  Quant100[[m]] = read.table(paste0(dirs[[m]],"/Quants.txt"), header = T)[ind100,(2:4)]
  Quant100[[m]]$model = models[m]; #Quant100[[m]]$sample = samples[s]
  Quant1000[[m]] = read.table(paste0(dirs[[m]],"/Quants.txt"), header = T)[ind1000,(2:4)]
  Quant1000[[m]]$model = models[m];# Quant1000[[m]]$sample = samples[s]

}


Quants100 = cast(melt.list(Quant100) );  Quants1000 = cast(melt.list(Quant1000) ) 
Quants100$model = as.factor(Quants100$model);Quants1000$model = as.factor(Quants1000$model)
Quants100$model =fct_relevel(Quants100$model,
                             as.vector(Quants100$model[c(1,3,2)]))
Quants1000$model=fct_relevel(Quants1000$model,
                             as.vector(Quants100$model[c(1,3,2)]))

############ BARPLOT QUANTILES ############ 
### Q100
BarQ100 = ggplot(data=Quants100, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.3, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = pals)+
  # ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
  theme_light(base_size = textsize)+
  ggtitle("Q100")+
  scale_x_discrete(labels=c("GEV","A","Censure"))
# scale_x_discrete(labels=c(rep("GEV",length(Quants100$model)-4),c("u0","uS","uN","uSN")))

BarQ1000 = ggplot(data=Quants1000, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.3, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = pals)+
  # ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  ylab(expression(paste("Débit [",m^3,"/",s,"]",sep="")))+
  theme_light(base_size = textsize)+
  # ggtitle(paste0("Q1000 - sample ", s) )+
  ggtitle("Q1000")+
  scale_x_discrete(labels=c("GEV","A","Censure"))

ggarrange(BarQ100+coord_cartesian(ylim = c(5000,max(Quants1000$Q_9)) )+
            theme(axis.title.x = element_blank(),legend.title = element_blank()),
          BarQ1000+coord_cartesian(ylim = c(5000,max(Quants1000$Q_9)) )+
            theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                  legend.title = element_blank()),
          common.legend = T, legend = "right")

ggsave(path = dir.plot, filename = paste0("Barplots_QX_censure.pdf"),
       width = 13,height = 8)



