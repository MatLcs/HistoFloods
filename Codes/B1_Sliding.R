#known nban / unkknown threshold
# case = "B1"
# dir.case = paste0(dir.res.ts,case,"_case/"); dir.create(dir.case, showWarnings = F)

Seuil = parameter(name = "Seuil", init = Ts, priorDist = "Gaussian", priorPar = c(Ts,sd.Ts))
NbAn = parameter(name = "NbAn", init = Anfin-Andeb, priorDist = "FIX")

#### Initializing DF for MCMC results
MegaSpag = data.frame()

for(spag in 1:Nspag){
  # spag = 1
  print(paste0("Spag = ",spag))
  #draw a random spag within the 500 Q AMAX spags
  Y = data.frame(  c(Spags[,sample(1:500,1)],#[,spag],
                     length(CX$An)) )
  Var = data.frame(factor(c(rep('Q',length(Q$mp)),"Occ")))
  dat <- dataset(Y = Y, var = Var)
  mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
               par=list(Pos,Ech,Form,Seuil,NbAn),
               ### PARAM GEV
               formula = c("Pos=Pos","Ech=Ech","Form = Form",
                           ### PARAM BINOMIAL
                           "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                           "NbAn = NbAn"))
  #### Create spag subfolder 
  dir.spag = file.path(dir.case,"Spag"); dir.create(dir.spag,showWarnings = F)
  STooDs(mod,workspace = dir.spag, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
  #### Read MCMC results for each spag
  MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
  MegaSpag = rbind(MegaSpag, MCMCres[,(1:5)])
}

MegaGev = matrix(nrow = nrow(MegaSpag),ncol = length(prob))
#### Compute the quantiles for all the GeV realisations spags
#### WARNING, HERE SHAPE HAS TO BE THE OPPOSITE THE SHAPE GIVEN BY STOODS (BY CONVENTION)
for ( i in 1 : nrow(MegaSpag)) {MegaGev[i,] = qgev(p = prob,loc = MegaSpag$Pos[i],
                                                   scale = MegaSpag$Ech[i], 
                                                   shape = -1*MegaSpag$Form[i]) }

PostForm = MegaSpag$Form
PostThres = MegaSpag$Seuil
PostNban = MegaSpag$NbAn
#### Compute the true maxpost quantiles : maxpost of hydro sample x maxpost of GeV estim
Y = data.frame(c(Q$mp, length(CX$An)) )
Var = data.frame(factor(c(rep('Q',length(Q$mp)),"Occ")))
dat <- dataset(Y = Y, var = Var)
#### Create MP sub-folder & run Stoods
dir.mp = file.path(dir.case,"/Mp") ; dir.create(dir.mp,showWarnings = F)
dat <- dataset(Y = Y, var = Var)
mod <- model(dataset=dat, parentDist=c("GEV","Binomial"), varName = c("Q","Occ") ,
             par=list(Pos,Ech,Form,Seuil,NbAn),
             ### PARAM GEV
             formula = c("Pos=Pos","Ech=Ech","Form = Form",
                         ### PARAM BINOMIAL
                         "Prob = 1 - exp((-1)*(1-Form*((Seuil-Pos)/Ech))^(1/Form))",
                         "NbAn = NbAn"))
STooDs(mod,workspace = dir.mp, mcmcOptions = mcmc(Nsim = Nsim, Nslim = 2))
MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
Mp.GeV = MCMCres[which.max(MCMCres$post),(1:5)]
Mp.Quant = qgev(p = prob,loc = Mp.GeV$Pos, scale = Mp.GeV$Ech,shape = -1*Mp.GeV$Form)

quant = apply(MegaGev,MARGIN = 2,FUN = quantile, probs = c(0.025,0.975))
## floods PR
Freq = Q[order(Q$mp),]
Freq$Fr = (seq(1:length(Q$an))-0.5)/length(Q$an)
Freq$Pr = 1/(1-Freq$Fr)

Quants = data.frame(Pr=Pr, Mp=Mp.Quant, Q_2=quant[1,] ,Q_9=quant[2,])

write.table(round(Quants,3), paste0(dir.case,case,"_Quants.txt"), row.names = F)
write.table(round(data.frame(Shape = PostForm, Thres = PostThres, Nban = PostNban),3),
            paste0(dir.case,case,"_Params.txt"), row.names = F)
write.table(Mp.GeV,  paste0(dir.case,case,"_Maxpost_Par.txt"), row.names = F)

# 
# Quant = ggplot()+
#   geom_ribbon(data=Quants,aes(x=Pr, ymin = Q_2, ymax=Q_9,
#                               fill="95% uncertainty interval"),alpha=0.8)+
#   geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"),lwd=1)+
#   geom_point(data = Freq, aes(x=Pr, y = mp))+
#   geom_errorbar(data=Freq, aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
#   scale_x_continuous(trans="log10")+
#   xlab("Return period [years]")+
#   ylab("Discharge [m3/s]")+
#   scale_fill_manual(name = element_blank(),
#                     values = c("#67a9cf"))+
#   scale_color_manual(values = c("yellow"))+ #"royalblue")+
#   theme_bw(base_size=15)+
#   # labs(title = paste0(head(Q$an,1)," - ",tail(Q.case$an,1)))+
#   coord_cartesian(xlim=c(1,max(Pr)))+
#   theme(legend.title=element_blank(),
#         plot.title = element_text(hjust = 0.01, vjust = -7),
#         legend.position = c(0.8,0.2))  
# 
# ggsave(Quant, path = dir.case, filename = paste0("Quantiles_",case,caseTs,".pdf"),width = 10, height = 7)
# 
# 


