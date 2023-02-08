rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
#### PARAMS SIMULATION
Andeb = 1500; Anfin = 1816; Nsim = 6000; Nspag = 500
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
#### Color palette
pal = c("grey",brewer.pal(4,"RdYlBu")[c(2,1,3,4)])
Quant100 = list();Quant1000 = list();ParamAll = list(); Mp_par = list(); QuantsAll = list()
ind100 = length(Pr) - 18;ind1000 = length(Pr) - 9

models = c("Baseline","A1","A1")
samples = c("C3","C3","C4")
Barplots = T
Params = T
Quants = T

#### READ RESULTS
for(res in 1:(length(models))){
  ## All quants
  QuantsAll[[res]] = read.table(paste0(dir.res,"/",samples[res],"/",models[res],"_case/",
                                        models[res],"_Quants.txt"), header = T)[,(2:4)]
  QuantsAll[[res]]$model = models[res]; QuantsAll[[res]]$sample = samples[res]
  ##Q100 & Q1000
  Quant100[[res]] = read.table(paste0(dir.res,"/",samples[res],"/",models[res],"_case/",
                                      models[res],"_Quants.txt"), header = T)[ind100,(2:4)]
  Quant100[[res]]$model = models[res];  Quant100[[res]]$sample = samples[res]
  Quant1000[[res]] = read.table(paste0(dir.res,"/",samples[res],"/",models[res],"_case/",
                                       models[res],"_Quants.txt"), header = T)[ind1000,(2:4)]
  Quant1000[[res]]$model = models[res];  Quant1000[[res]]$sample = samples[res]
  ## posterior params
  ParamAll[[res]] = read.table(paste0(dir.res,"/",samples[res],"/",models[res],"_case/",
                                       models[res],"_Params.txt"), header = T)
  ParamAll[[res]]$model = models[res]; ParamAll[[res]]$sample = samples[res]
  ## maxpost params
  Mp_par[[res]] = read.table(paste0(dir.res,"/",samples[res],"/",models[res],"_case/",
                                    models[res],"_Maxpost_Par.txt"),header = T)[1,]
  Mp_par[[res]]$model = models[res]; Mp_par[[res]]$sample = samples[res]
}

### RESULTS TIDYING
Params = melt.list(ParamAll); Quants100 = cast(melt.list(Quant100)); 
Quants1000 = cast(melt.list(Quant1000)) ; Params$model = as.factor(Params$model)

# Params$sample = fct_relevel(Params$case, c("Baseline","A1","A2","B1","B2") )
Params$value[which(Params$variable == "Thres" & 
                     (Params$model == "A1" | Params$model == "A2"))] = NA
Params$value[which(Params$variable == "Nban" & Params$model != "B2")] = NA
ParamsMp = melt.list(Mp_par) 
# ParamsMp$case = NA
# for(i in 1:length(cases)){ParamsMp$case[which(ParamsMp$L1 == i)] = cases[i]}
# ParamsMp$value = round(ParamsMp$value,3) ; ParamsMp$case = as.factor(ParamsMp$case)
# ParamsMp$case = fct_relevel(ParamsMp$case, c("Baseline","A1","A2","B1","B2") )

############ BARPLOT QUANTILES ############ 
### Q100
Quants100$model = as.factor(Quants100$model)
Quants100$model = fct_relevel(Quants100$model, c("Baseline","A1","A2","B1","B2") )
BarQ100 = ggplot(data=Quants100, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = pal)+
  ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  theme_light()+
  ggtitle(paste0("Q100: sample", samples))

### Q1000
Quants1000$model = as.factor(Quants1000$model)
Quants1000$model = fct_relevel(Quants1000$model, c("Baseline","A1","A2","B1","B2") )
BarQ1000 = ggplot(data=Quants1000, aes(x = model, y = Mp, ymin = Q_2, ymax = Q_9, fill = model))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = pal)+
  ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
  theme_light()+
  ggtitle(paste0("Q1000: sample", samples))

#### ARRANGE
ggarrange(BarQ100+coord_cartesian(ylim = c(5000,17000))+
            theme(axis.title.x = element_blank()),
          BarQ1000+coord_cartesian(ylim = c(5000,17000))+
            theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
          common.legend = T, legend = "right")
# ggsave(path = dir.plots.ts, filename = paste0("Barplot_QX_",casesTs[t],".pdf"), width = 12,height = 8)


