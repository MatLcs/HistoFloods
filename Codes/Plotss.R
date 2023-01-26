rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")

cases = c("HistoC3","HistoC4","Artif8000","Artif10000")
dirC3 = paste0(dir.res,"C3/")
dirC4 = paste0(dir.res,"C4/")
dir8000 = paste0(dir.res,"HistoArtif70/8000/")
dir10000 = paste0(dir.res,"HistoArtif70/10000/")
scenarios = c("A1","A2","B1","B2")
dirs = c(dirC3,dirC4,dir8000,dir10000)
# RP up top 10 000 years
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)

Quant100 = list();Quant1000 = list()
ind100 = length(Pr) - 18
ind1000 = length(Pr) -9
Base100 = read.table(paste0(dirC3,"Baseline_case/Baseline_Quants.txt"),
           header = T)[ind100,(2:4)]
Base100$case = "Baseline";  Base100$scenario = "Baseline"
Base1000 = read.table(paste0(dirC3,"Baseline_case/Baseline_Quants.txt"),
                     header = T)[ind1000,(2:4)]
Base1000$case = "Baseline";  Base1000$scenario = "Baseline"
cpt = 0
for(case in 1:length(dirs)){
  for(scen in 1:length(scenarios)){
  cpt=cpt+1
  Quant100[[cpt]] = read.table(paste0(dirs[case],scenarios[scen],"_case/",scenarios[scen],"_Quants.txt"),
                                header = T)[ind100,(2:4)]
  Quant100[[cpt]]$case = cases[case];  Quant100[[cpt]]$scenario = scenarios[scen]
  Quant1000[[cpt]] = read.table(paste0(dirs[case],scenarios[scen],"_case/",scenarios[scen],"_Quants.txt"),
                             header = T)[ind1000,(2:4)]
  Quant1000[[cpt]]$case = cases[case]  ; Quant1000[[cpt]]$scenario = scenarios[scen]
  }
}
Quant100[[17]] = Base100
Quant1000[[17]] = Base1000

Quants100 = cast(melt.list(Quant100)); 
Quants1000 = cast(melt.list(Quant1000)) 
Quants100$case = as.factor(Quants100$case)
Quants100$case = fct_relevel(Quants100$case, c("Baseline",cases) )
Quants1000$case = as.factor(Quants1000$case)
Quants1000$case = fct_relevel(Quants1000$case, c("Baseline",cases) )

pals = list(brewer.pal(5,"YlOrBr")[2:5],brewer.pal(5,"Reds")[2:5],
            brewer.pal(5,"PuBu")[2:5],brewer.pal(5,"Blues")[2:5])

for(scen in 1:length(scenarios)){
  inds = which(Quants100$scenario==scenarios[scen] | Quants100$scenario=="Baseline")
  ### Q100
    BarQ100 = ggplot(data=Quants100[inds,]
                     , aes(x = case, y = Mp, ymin = Q_2, ymax = Q_9, fill = case))+
      geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
      geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
      scale_fill_manual(values = c("grey",pals[[scen]]))+
      ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
      theme_light()+
      ggtitle(paste0("case ",scenarios[scen],": Q100"))
    ### Q1000
    BarQ1000 = ggplot(data=Quants1000[inds,]
                     , aes(x = case, y = Mp, ymin = Q_2, ymax = Q_9, fill = case))+
      geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
      geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
      scale_fill_manual(values = c("grey",pals[[scen]]))+
      ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
      theme_light()+
      ggtitle(paste0("case ",scenarios[scen],": Q1000"))
    
    #### ARRANGE
    ggarrange(BarQ100+coord_cartesian(ylim = c(5000,25000))+
                theme(axis.title.x = element_blank()),
              BarQ1000+coord_cartesian(ylim = c(5000,25000))+
                theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
              common.legend = T, legend = "right")
    ggsave(path = dir.plots, filename = paste0("Barplot_",scenarios[scen],"_all",".pdf"),
           width = 12,height = 8)

}