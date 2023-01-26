rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
#### PARAMS SIMULATION
Andeb = 1500; Anfin = 2000; Nsim = 4000; Nspag = 500
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
#### Run or plot only ?
RUN = F
### Threshold hypotheses
casesTs = c("C3","C4")
#### 5 simulation cases
cases = c("Baseline","A1","A2","B1","B2")
### 1816-2020 AMAX discharges with uncertainty and spaghetti 
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
Q = read.table(
  "C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]
###Histrhone data since 1500
source(paste0(dir.codes,"FloodsCX.R"))
Call = read.csv2(paste0(dir.data,"CX_All.csv"))
#### Color palette
pal = c("grey",brewer.pal(4,"RdYlBu")[c(2,1,3,4)])
#### Results reading init
Quant100 = list(); Quant1000 = list(); ParamAll = list(); Mp_par = list(); QuantsAll = list()
### Prior GEV params, common to all sim
Pos = parameter(name='Pos',init = 6000) ; Ech =  parameter(name='Ech',init = 1000) 
Form = parameter(name='Form',init = 0.01,priorDist='Gaussian',priorPar=c(0,0.3))
#### RUN 
if(RUN == T){
  #for C3 and C4 thresholds
  for(t in 1:length(casesTs)){
    #for all cases
    for(i in 1:length(cases)){
      caseTs = casesTs[t]
      dir.plots.ts = paste0(dir.plots,caseTs,"/")
      dir.res.ts = paste0(dir.res,caseTs,"/")
      dir.create(dir.plots.ts,showWarnings = F)
      dir.create(dir.res.ts,showWarnings = F)
      Thresh = read.table(paste0(dir.data,"Threshold",caseTs,".txt"), header = T)
      #maxpost or tot2.5?
      Ts = Thresh$tot2.5
      sd.Ts = (Thresh$tot97.5 - Thresh$mp)/2
      #histoflood
      CX = Call[which(Call$An > Andeb),]
      if(caseTs == "C4"){ CX = CX[which(CX$Cat == caseTs),] }
      #nban
      sd.NbAn = 50
      
      #stoods execution
      source(paste0(dir.codes,cases[i],".R"))  
      beep()
      }  
    }
}
  
Freq = Q[order(Q$mp),]
Freq$Fr = (seq(1:length(Q$an))-0.5)/length(Q$an)
Freq$Pr = 1/(1-Freq$Fr)

ind100 = length(Pr) - 18
ind1000 = length(Pr) -9

#### reading results
for(t in 1:length(casesTs)){
  
  GGQuants = list()
  caseTs = casesTs[t]
  dir.plots.ts = paste0(dir.plots,caseTs,"/")
  dir.res.ts = paste0(dir.res,caseTs,"/")

  for(case in 1: length(cases)){
    QuantsAll[[case]] = read.table(paste0(dir.res.ts,cases[case],"_case/",cases[case],"_Quants.txt"),
                                   header = T)[,(2:4)]
    Quant100[[case]] = read.table(paste0(dir.res.ts,cases[case],"_case/",cases[case],"_Quants.txt"),
                                  header = T)[ind100,(2:4)]
    Quant100[[case]]$case = cases[case]
    Quant1000[[case]] = read.table(paste0(dir.res.ts,cases[case],"_case/",cases[case],"_Quants.txt"),
                                   header = T)[ind1000,(2:4)]
    Quant1000[[case]]$case = cases[case]
    
    ParamAll[[case]] = read.table(paste0(dir.res.ts,cases[case],"_case/",cases[case],"_Params.txt"),
                                  header = T)
    ParamAll[[case]]$case = cases[case]
    
    Mp_par[[case]] = read.table(paste0(dir.res.ts,cases[case],"_case/",cases[case],"_Maxpost_Par.txt"),
                               header = T)[1,]
  }
  
  ### RESULTS TIDYING
  Params = melt.list(ParamAll); Quants100 = cast(melt.list(Quant100)); 
  Quants1000 = cast(melt.list(Quant1000)) ; Params$case = as.factor(Params$case)
  Params$case = fct_relevel(Params$case, c("Baseline","A1","A2","B1","B2") )
  Params$value[which(Params$variable == "Thres" & (Params$case == "A1" | Params$case == "A2"))] = NA
  Params$value[which(Params$variable == "Nban" & Params$case != "B2")] = NA
  ParamsMp = melt.list(Mp_par) ; ParamsMp$case = NA
  for(i in 1:length(cases)){ParamsMp$case[which(ParamsMp$L1 == i)] = cases[i]}
  ParamsMp$value = round(ParamsMp$value,3) ; ParamsMp$case = as.factor(ParamsMp$case)
  ParamsMp$case = fct_relevel(ParamsMp$case, c("Baseline","A1","A2","B1","B2") )
  
  ############ PLOT POSTERIOR PARAMETERS ############ 
  ####  SHAPE
  Shape = ggplot(Params[which(Params$variable=="Shape"),], 
                 aes(y = case, x = value, group = case, fill = case)) +
    geom_density_ridges(alpha=0.7,stat = "binline",bins = 60, scale = 1.2,color = NA)+
    geom_segment(data = ParamsMp[which(ParamsMp$variable == "Form"),],
                 aes(x = value, xend = value,
                     y = as.numeric(case), yend = as.numeric(case) + 1,
                     color = case),
                     lwd = 0.8)+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_discrete(expand = c(0, 0)) +
    theme_light()+
    scale_fill_manual(values = pal)+
    scale_color_manual(values = pal)+
    xlab("[-]")+
    theme(axis.title.y = element_blank())+
    ggtitle("GEV Shape parameter")
  
  #### Threshold
  Thres = ggplot(Params[which(Params$variable == "Thres" ),], 
         aes(x = value, fill = case)) +
    geom_histogram(position = "identity", bins = 60, alpha = 0.7, color = NA)+
    geom_vline(data = ParamsMp[which(ParamsMp$variable == "Seuil"),],
               aes( xintercept = value, color = case), lwd = 0.8, alpha = 1)+
    theme_light()+
    scale_fill_manual(values = pal[4:5])+
    scale_color_manual(values = pal)+
    coord_cartesian(ylim = c(0, (Nspag*Nsim/ 60 * 0.5)+1000) )+
    xlab("Discharge [m3/s]")+
    theme(axis.title.y = element_blank())+
    ggtitle("Perception threshold")
  
  #### NB AN
  Nban = ggplot(Params[which(Params$variable == "Nban"),], 
         aes(x = value, fill = case)) +
    geom_histogram(position = "identity", bins = 60, alpha = 0.7,color = NA)+
    geom_vline(data = ParamsMp[which(ParamsMp$variable == "NbAn" & ParamsMp$L1 != 1),],
               aes( xintercept = value, color = case), lwd = 0.8, alpha = 1)+
    theme_light()+
    scale_fill_manual(values = pal[5])+
    scale_color_manual(values = pal[2:5])+
    coord_cartesian(ylim = c(0,Nspag*Nsim/ 60 * 0.5))+
    xlab("Number of years")+
    theme(axis.title.y = element_blank())+
    ggtitle("Historical series length")
  
  ##### ARRANGE
  ggarrange(Shape,Thres,Nban, ncol = 3, common.legend = T, legend = "right", align = "h")
  ggsave(path = dir.plots.ts, paste0("Post_params_",casesTs[t],".pdf"), width = 14,height = 5)
  
  ############ PLOT POSTERIOR QUANTILES ############ 
  ### Q100
  Quants100$case = as.factor(Quants100$case)
  Quants100$case = fct_relevel(Quants100$case, c("Baseline","A1","A2","B1","B2") )
  BarQ100 = ggplot(data=Quants100, aes(x = case, y = Mp, ymin = Q_2, ymax = Q_9, fill = case))+
    geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
    geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
    scale_fill_manual(values = pal)+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    theme_light()+
    ggtitle("Q100")
    
  ### Q1000
  Quants1000$case = as.factor(Quants1000$case)
  Quants1000$case = fct_relevel(Quants1000$case, c("Baseline","A1","A2","B1","B2") )
  BarQ1000 = ggplot(data=Quants1000, aes(x = case, y = Mp, ymin = Q_2, ymax = Q_9, fill = case))+
    geom_bar( stat="identity",width = 0.8, alpha = 0.7)+
    geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
    scale_fill_manual(values = pal)+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    theme_light()+
    ggtitle("Q1000")
  
  #### ARRANGE
  ggarrange(BarQ100+coord_cartesian(ylim = c(5000,17000))+
              theme(axis.title.x = element_blank()),
            BarQ1000+coord_cartesian(ylim = c(5000,17000))+
              theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
            common.legend = T, legend = "right")
  ggsave(path = dir.plots.ts, filename = paste0("Barplot_QX_",casesTs[t],".pdf"), width = 12,height = 8)
  
  
  #### Quants
  for(case in 1:length(cases)){

      GGQuants[[case]] = ggplot()+
        geom_ribbon(data=QuantsAll[[case]],aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                    fill="95% uncertainty interval"),alpha=0.8)+
        geom_line(data=QuantsAll[[case]],aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
        geom_point(data = Freq, aes(x=Pr, y = mp))+
        geom_errorbar(data=Freq, aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
        scale_x_continuous(trans="log10")+
        xlab("Return period [years]")+
        ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
        scale_fill_manual(name = element_blank(),values = c("#67a9cf"))+
        scale_color_manual(values = c("yellow"))+ 
        theme_bw(base_size=15)+
        coord_cartesian(xlim=c(1,max(Pr)))+
        theme(legend.title=element_blank(),
              plot.title = element_text(hjust = 0.01, vjust = -7),
              legend.position = c(0.8,0.2))+
        ggtitle(cases[case])+
        coord_cartesian(ylim = c(2000,20000))
  
  }
  
  ggarrange(GGQuants[[1]],#+theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
            GGQuants[[2]],#+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                                #axis.text.y = element_blank(),axis.text.x = element_blank()),
            GGQuants[[3]],#+theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
            GGQuants[[4]],#+theme(axis.title.y = element_blank(),axis.text.y = element_blank()),
            GGQuants[[5]],
            ncol = 2, nrow = 3,
            common.legend = T, legend = "bottom",align = "hv")
  ggsave(path = dir.plots.ts, filename = paste0("Quantiles_",casesTs[t],".pdf"),
         width = 15, height = 17)
  
  
  
}



