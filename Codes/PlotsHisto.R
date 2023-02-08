rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
#### PARAMS SIMULATION
Andeb = 1500; Anfin = 1816; Nsim = 6000; Nspag = 500
prob = c(seq(0.01,0.999,0.001),seq(0.9991,0.9999,0.0001)) ; Pr = 1/(1-prob)
#### Run or plot only ?
RUN = T
### Threshold hypotheses
casesTs = c("C3","C4")
Thresh = c(7000,9000)
sdThresh = c(2000,2000)
### Uniform distribution of NbAn c(lastyear - firstfloodyear , 1000)
distNbAn = c(1816-1511,1000)
#### 5 simulation cases
cases = c("Baseline","A1","A2","B1","B2")
###Histrhone data since 1500
# source(paste0(dir.codes,"FloodsCX.R"))
Call = read.csv2(paste0(dir.data,"CX_All.csv"))
### 1816-2020 AMAX discharges with uncertainty and spaghetti 
Spags = read.table(paste0(dir.data,"Spags_uTot_Amax.txt"))
Q = read.table(
  "C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Results/Quantiles_Amax.txt",
  header = T)[,c(1,2,7,8)]
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
      # Thresh = read.table(paste0(dir.data,"Threshold",caseTs,".txt"), header = T)
      # #maxpost or tot2.5?
      Ts = Thresh[t]
      sd.Ts = sdThresh[t]
      #histoflood
      CX = Call[which(Call$An > Andeb & Call$An < Anfin),]
      if(caseTs == "C4"){ CX = CX[which(CX$Cat == caseTs),] }
      #nban
      # sd.NbAn = 50
      
      #stoods execution
      source(paste0(dir.codes,cases[i],".R"))  
      beep()
      }  
    }
}
  
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
  
  #to add prior shape to plot
  # priorshape = Params[which(Params$variable=="Shape"),(1:3)]
  # priorshape$variable = "PriorShape"
  # priorshape$value = rep(
  #                       rnorm( 
  #                         (length(Params[which(Params$variable=="Shape"),1])/length(cases)),
  #                         Form$priorPar[1],Form$priorPar[2] ), length(cases))
  ############ PLOT POSTERIOR PARAMETERS ############ 
  ####  SHAPE
  Shape = ggplot(Params[which(Params$variable=="Shape"),], 
                 aes(y = case, x = value, group = case, fill = case)) +
    
    # geom_density_ridges(data=priorshape,aes(y = case, x = value, group = case), fill = "lightgrey",
    #                     alpha=0.5,stat = "binline",bins = 60, scale = 1.2,color = NA )+
  
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
    # +coord_cartesian(xlim = c(-0.2,0.2))

  
  #### Threshold
  priorthres = data.frame( case = "1-Prior", variable= "Thres", value= rnorm(
    (length(Params[which(Params$variable=="Thres"),1])/length(cases)),
    Thresh[t],sdThresh[t] ),  L1 = NA)
  
  Thres = ggplot( rbind(priorthres,Params[which(Params$variable == "Thres" ),]), 
         aes(x = value, fill = case)) +
    geom_histogram(position = "identity", bins = 200, alpha = 0.7, color = NA)+
    geom_vline(data = ParamsMp[which(ParamsMp$variable == "Seuil"),],
               aes( xintercept = value, color = case), lwd = 0.8, alpha = 1)+
    theme_light()+
    # scale_fill_manual(values = pal[4:5])+
    scale_fill_manual(values = c("lightgrey",pal[c(4,5)]))+
    scale_color_manual(values = pal)+
    coord_cartesian(xlim = c(Thresh[t]-1.5*sdThresh[t],Thresh[t]+1.5*sdThresh[t]) )+
                    # ylim = c(0, (Nspag*Nsim)/(200*0.35)) )+
    xlab("Discharge [m3/s]")+
    theme(axis.title.y = element_blank())+
    ggtitle("Perception threshold")
  
  #### NB AN
  priornban = data.frame( case = "1-Prior", variable= "Nban", value= runif(
    n=(length(Params[which(Params$variable=="Nban"),1])/length(cases)),
    min=distNbAn[1],max=distNbAn[2] ),  L1 = NA)
  
  Nban = ggplot(rbind(priornban,Params[which(Params$variable == "Nban"),]), 
         aes(x = value, fill = case)) +
    geom_histogram(position = "identity", bins = 100, alpha = 0.7,color = NA)+
    geom_vline(data = ParamsMp[which(ParamsMp$variable == "NbAn" & ParamsMp$L1 != 1),],
               aes( xintercept = value, color = case), lwd = 0.8, alpha = 1)+
    theme_light()+
    scale_fill_manual(values = c("lightgrey",pal[5]))+
    scale_color_manual(values = pal[2:5])+
    # coord_cartesian(ylim = c(0, (Nspag*Nsim)/(100*0.35)) )+
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
    
      
    if(case == 1){
      Freq = Q[order(Q$mp),]
      Freq$Fr = (seq(1:length(Q$an))-0.5)/length(Q$an)
      Freq$Pr = 1/(1-Freq$Fr)
      
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
    } else {
      Ts = Thresh[t]
      CX = Call[which(Call$An > Andeb & Call$An < Anfin),]
      if(caseTs == "C4"){ CX = CX[which(CX$Cat == caseTs),] }
      #### plot positions histo
      NC = length(Q$an)
      NH = Anfin-Andeb
      #Nombre de crues histo sup. au seuil
      NS_H = length(CX$An)
      #Nombre de crues de la période continue sup. au seuil
      NS_C = length(which(Q$mp>Ts))
      #Nombre tot. de crues sup. au seuil 
      NS = NS_H + NS_C
      #Nombre de crues inf. au seuil
      Ninf = length(which(Q$mp <= Ts))
      #Nombre d'années total
      NAn = NC + NH
      #proba crue continue
      PC = NS_C/NS
      #proba crue histo
      PH = NS_H/NS
      #data frame crues supseuil
      SupSeuil = data.frame(rang=1:NS)
      SupSeuil$flag = "Hist"
      set.seed(2020)
      SupSeuil$prob = runif(NS)
      SupSeuil$flag[which(SupSeuil$prob < PC)] = "Cont"
      SupSeuil$mp = NA
      SupSeuil$tot2.5 = NA
      SupSeuil$tot97.5 = NA
      #bidouille pour se limiter à NS_C crues "Cont"
      SupSeuil$flag[which(SupSeuil$rang > which(SupSeuil$flag == "Cont")[NS_C+1])]="Hist"
      SupSeuil$mp[which(SupSeuil$flag == "Cont")] = Q$mp[order(Q$mp,decreasing = T)]
      SupSeuil$tot2.5[which(SupSeuil$flag == "Cont")] = Q$tot2.5[order(Q$mp,decreasing = T)]
      SupSeuil$tot97.5[which(SupSeuil$flag == "Cont")] = Q$tot97.5[order(Q$mp,decreasing = T)]
      #df crues sous seuil
      SousSeuil = Q[which(Q$mp<Ts),(2:4)]
      SousSeuil = SousSeuil[order(SousSeuil$mp, decreasing = T),]
      SousSeuil$rang = (NS+1):(NS_H+NC)
      SousSeuil$flag = "Cont"
      #Freq
      Freq = merge.data.frame(SupSeuil[,c(1,2,4,5,6)],SousSeuil, all = T)
      Freq$tot2.5[which(Freq$flag == "Hist")] = Ts
      Freq$tot97.5[which(Freq$flag == "Hist")] = 50000
      Freq$Fr = (Freq$rang-0.5)/(NS+Ninf)
      Freq$Pr = 1/Freq$Fr
      
      ### PLOT QUANTS
      GGQuants[[case]] = ggplot()+
        geom_ribbon(data=QuantsAll[[case]],aes(x=Pr, ymin = Q_2, ymax=Q_9,
                                    fill="95% uncertainty interval"),alpha=0.8)+
        geom_line(data=QuantsAll[[case]],aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
        geom_errorbar(data=Freq[which(Freq$flag=="Cont"),],
                      aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
        geom_point(data = Freq[which(Freq$flag=="Cont"),],
                   aes(x=Pr, y = mp))+
        geom_errorbar(data=Freq[which(Freq$flag=="Hist"),],
                      aes(x=Pr,ymin = tot2.5, ymax = tot97.5),lty=2)+
        geom_point(data = Freq[which(Freq$flag=="Hist"),], aes(x=Pr,y=tot2.5),shape = 17)+
        
        scale_x_continuous(trans="log10")+
        xlab("Return period [years]")+
        ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
        scale_fill_manual(name = element_blank(),
                          values = c("#67a9cf"))+
        scale_color_manual(values = c("yellow"))+ #"royalblue")+
        theme_bw(base_size=15)+
        # labs(title = paste0(head(Q$an,1)," - ",tail(Q.case$an,1)))+
        coord_cartesian(xlim=c(1,max(Pr)),ylim = c(2000,20000))+
        theme(legend.title=element_blank(),
              plot.title = element_text(hjust = 0.01, vjust = -7),
              legend.position = c(0.8,0.2))
      
    }
  
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



