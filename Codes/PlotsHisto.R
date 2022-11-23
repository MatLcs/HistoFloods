rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/CruesHisto/Codes/Dirs&Funs.R")
source(paste0(dir.codes,"FloodsCX.R"))

#### 5 cases
cases = c("Baseline","A1","A2","B1","B2")
#### Color palette
pal = c("grey",brewer.pal(4,"RdYlBu"))
#### Results reading init
Quant100 = list()
Quant1000 = list()
ParamAll = list()
Mp_par = list()

#### Run or plot only ?
RUN = F
#### RUN 
if(RUN == T){
  for(i in 1:length(cases)){source(paste0(dir.codes,cases[i],".R"))  }  }

ind100 = length(Pr) - 10
ind1000 = length(Pr) 

#### reading results
for(case in 1: length(cases)){
  Quant100[[case]] = read.table(paste0(dir.res,cases[case],"_case/",cases[case],"_Quants.txt"),
                                header = T)[ind100,(2:4)]
  Quant100[[case]]$case = cases[case]
  Quant1000[[case]] = read.table(paste0(dir.res,cases[case],"_case/",cases[case],"_Quants.txt"),
                                 header = T)[ind1000,(2:4)]
  Quant1000[[case]]$case = cases[case]
  ParamAll[[case]] = read.table(paste0(dir.res,cases[case],"_case/",cases[case],"_Params.txt"),
                                header = T)
  ParamAll[[case]]$case = cases[case]
  
  Mp_par[[case]] = read.table(paste0(dir.res,cases[case],"_case/",cases[case],"_Maxpost_Par.txt"),
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
ggsave(path = dir.plots, "Post_params.pdf", width = 14,height = 5)

############ PLOT POSTERIOR QUANTILES ############ 
### Q100
Quants100$case = as.factor(Quants100$case)
Quants100$case = fct_relevel(Quants100$case, c("Baseline","A1","A2","B1","B2") )
BarQ100 = ggplot(data=Quants100, aes(x = case, y = Mp, ymin = Q_2, ymax = Q_9, fill = case))+
  geom_bar( stat="identity",width = 0.8,alpha = 0.7)+
  geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = pal)+
  ylab("Discharge [m3/s]")+
  theme_light()+
  ggtitle("Q100")
  
### Q1000
Quants1000$case = as.factor(Quants1000$case)
Quants1000$case = fct_relevel(Quants1000$case, c("Baseline","A1","A2","B1","B2") )
BarQ1000 = ggplot(data=Quants1000, aes(x = case, y = Mp, ymin = Q_2, ymax = Q_9, fill = case))+
  geom_bar( stat="identity",width = 0.8, alpha = 0.7)+
  geom_errorbar(width=0.2, colour="black", alpha=0.6, size=0.6)+
  scale_fill_manual(values = pal)+
  ylab("Discharge [m3/s]")+
  theme_light()+
  ggtitle("Q1000")

#### ARRANGE
ggarrange(BarQ100+coord_cartesian(ylim = c(5000,17000))+
            theme(axis.title.x = element_blank()),
          BarQ1000+coord_cartesian(ylim = c(5000,17000))+
            theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
          common.legend = T, legend = "right")
ggsave(path = dir.plots, filename = "Barplot_QX_cases.pdf", width = 12,height = 8)





