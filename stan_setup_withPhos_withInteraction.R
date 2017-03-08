#=======================
# Phylogeny is a powerful tool for predicting plant biomass responses to nitrogen enrichment. 
# Rachel C. Wooliver, Zachary H. Marion, Christopher R. Peterson, Brad M. Potts, John K. Senior, Joseph K. Bailey, and Jennifer A. Schweitzer
#=======================


library(ape)
library(shinystan)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(yarrr)
library(randomcoloR)



setwd("~/Desktop/Bayesian_Hierarchical_model")






#=======================
# Data
#=======================
eucal <- read.csv("EucsGH.csv", stringsAsFactors=FALSE)
eucal$n.g.m2.yr <- eucal$n.kg.ha.yr/10
eucal$p.g.m2.yr <- eucal$p.kg.ha.yr/10

tempDat <- eucal[,c("sp","dead","subgenus","n.g.m2.yr","p.g.m2.yr","sla.cm2.g","srl.cm.g","ssd.mg.mm3","total.biomass.g","ht.cm")]
o <- tempDat$sp!= "E. morrisbyi" 
tempDat1 <- tempDat[o,]
p <- tempDat1$dead!= "y"
tempDat2 <- tempDat1[p,] # 931 individuals after removing morrisbyi and dead


#=======================
# Phylogeny
#=======================
tree29 <- read.tree("Euc.final.dated.tre")
tree29$tip.label <- paste("E.", tree29$tip.label, sep=" ")
tree <- drop.tip(tree29, c("E. morrisbyi","E. nitida","E. vernicosa" ,"E. johnstonii","E. perriniana","E. nitens" ))
tree$tip.label[8] <- "E. tenuiramis"

tree.transformed <- tree 
tree.transformed $edge.length[c(1,22)] <- tree.transformed $edge.length[c(1,22)]-0.96

phylocor <- vcv.phylo(tree,corr=TRUE)
corNames <- order(colnames(phylocor))
phyCor <- phylocor[corNames, corNames]

#=======================
# Model
#=======================
tempBiom <- tempDat2[, c("sp", "n.g.m2.yr", "p.g.m2.yr", "total.biomass.g", "sla.cm2.g", "srl.cm.g", "ssd.mg.mm3")]

tempBiom2 <- na.omit(tempBiom) # 722 total individuals
naDrop = which(is.na(tempBiom[,4]))

Biom <- tempBiom2$total.biomass.g 

nitroBiom <- tempBiom2$n.g.m2.yr # nitrogen treatment
phosBiom <- tempBiom2$p.g.m2.yr 
phosBiomDiscrete <- phosBiom; phosBiomDiscrete[phosBiomDiscrete>0] <- 1  # phosphorus treatment

sppBiom <- as.integer(as.factor(tempBiom2 $sp)) # species identity as integers: in alph order
funTraits <- log(tempBiom2[,5:7])

indMatBiom <- cbind(rep(1, length(nitroBiom)), nitroBiom, phosBiomDiscrete, nitroBiom*phosBiomDiscrete)

groupVecBiom <- rep(1, length(unique(sppBiom)))

groupMatBiom <- matrix(rep(1, length(unique(sppBiom))), nrow=length(unique(sppBiom)), ncol=1)

UmatBiom <- cbind(groupVecBiom)
UmatFun <- rbind(groupVecBiom, groupVecBiom, groupVecBiom)


biomData <- list(
  nSpp = length(unique(sppBiom)), # 23 species
  KBiom = ncol(indMatBiom), # 4 individual level predictors
  nObsBiom = length(Biom), # 722 total individuals
  nFun = ncol(funTraits), # 3 functional traits
  UmatBiom = UmatBiom[,1], # species predictor matrix (biomass)
  sppBiom = sppBiom, # species identities
  x = indMatBiom, # individual predictor matrix
  funTraits = funTraits, # functional trait matrix
  ObsBiom = Biom, # individual biomass
  LKJParam = 2, 
  phyloCor = phyCor,
  d = diag(diag(phyCor)))

logBiom = log(Biom)


 fit_full_withInteraction <- stan(file.path(getwd(), "HBM_fit_full_withInteraction.stan"), 
                                  data= biomData, iter=2000, chains=3, seed=5, 
                                  control = list(adapt_delta = .9, max_treedepth = 13), cores = 3)  



####
# Extract data
####

nospp <- 23
species <- data.frame(no=unique(as.integer(as.factor(tempBiom2 $sp))), sp=unique(as.character(tempBiom2 $sp))); rownames(species) <-unique(as.character(tempBiom2 $sp)); species <- species[tree$tip.label,]
lambdaStat <- monitor(extract(fit_full_withInteraction, pars=c("lambda"),permuted=FALSE),print=FALSE)[,c("mean", "2.5%","50%", "97.5%","Rhat")]

overallBetaEsts <-as.data.frame(monitor(extract(fit_full_withInteraction, pars=c("beta"),permuted=FALSE),print=FALSE)[,c("mean", "2.5%","50%", "97.5%","Rhat")])
rawBeta <- as.matrix(fit_full_withInteraction, pars=c("beta"))


  
  ########################################################################
  ###### PLOT REACTION NORMS FOR SPECIES/LINEAGES/SUBGENERA IN LOW VS. HIGH P #####
  ##### Revised
  
  pdf("Ecology_Figure 2_revised.pdf", width=6, height=7)
  layout(matrix(c(1,2,3,4,5,6),3,2, byrow = T), widths=c(1,1,0.5))
  par(oma=c(4.5,4,2,6.5), mar=c(0,0,0,0))
  nospp <- 23
  
  ######## Subgenus responses to N in low P ########
  
  overallBetaEsts$subgenus <- rep(c("E", "S", "S", "E", "S", "S", "E", "S", "S",
                                    "E", "S", "E", "E", "E", "E", "E", "S", "S", "E",
                                    "S", "E", "S", "S"), 4)
  
  overallBetaEstsSubInt <- aggregate(mean ~ subgenus, data=overallBetaEsts[c(1:nospp),], FUN="mean")
  overallBetaEstsSubSlope1 <- aggregate(mean ~ subgenus, data=overallBetaEsts[c(1:nospp)+nospp,], FUN="mean")
  overallBetaEstsSubSlope2 <- aggregate(mean ~ subgenus, data=overallBetaEsts[c(1:nospp)+nospp*2,], FUN="mean")
  overallBetaEstsSubSlope3 <- aggregate(mean ~ subgenus, data=overallBetaEsts[c(1:nospp)+nospp*3,], FUN="mean")
  
  xs<-seq(0,10,length.out = 1000)
  getys<-function(intercept,slope){
    predictedy<-NA
    for(i in 1:1000){
      predictedy[i] <- exp(intercept + slope*xs[i])
    }
    return(predictedy)
  }
  
  plot(0,0, ylim=c(0,40), xlim=c(0,10.5), type="n", las=1, 
       mgp=c(2, 0.65, 0),
       xlab="", xaxt="n")
  mtext("Biomass (g)", 2, line=2.25, cex=0.9)
  mtext(expression(paste("0 g P m"^-2, " yr"^-1)), 3, line=0, cex=0.8)
  abline(h=seq(0,40,5), col=scales::alpha("gray",0.4), lty=3)
  
  cols <- c("purple", "darkorange1")
    betaIters <- extract(fit_full_withInteraction, "beta")$beta 
    lineageMCdat <- data.frame(Euc=rep(NA, 3000),
                               Symph=rep(NA, 3000))
    for(i in 1:3000){
      lineageMCdat$Euc[i] <- mean(betaIters[i,species$no[1:11],2])
      lineageMCdat$Symph[i] <- mean(betaIters[i,species$no[12:23],2])
    }
    (exp(quantile(lineageMCdat$Euc, probs = 0.1))-1)*100
    (exp(quantile(lineageMCdat$Symph, probs = 0.1))-1)*100
  ltys <- c(1,1)
  for(i in 1:2){
    lines(xs,getys(overallBetaEstsSubInt$mean[i], 
                   overallBetaEstsSubSlope1$mean[i]),type="l", col=cols[i], lty=ltys[i])
  }
  
  ######## Subgenus responses to N in high P ########
  
  
  getysInt<-function(intercept,slope1,slope2,interact){
    predictedy<-NA
    for(i in 1:1000){
      predictedy[i] <- exp(intercept + slope1*xs[i] + slope2 + interact*xs[i])
    }
    return(predictedy)
  }
  
  plot(0,0, ylim=c(0,40), xlim=c(0,10.5), type="n", las=1, 
       mgp=c(2.25, 0.65, 0),
       xlab="", xaxt="n", yaxt="n")
  mtext(expression(paste("1.2 g P m"^-2, " yr"^-1)), 3, line=0, cex=0.8)
  abline(h=seq(0,40,5), col=scales::alpha("gray",0.4), lty=3)
  
    betaIters <- extract(fit_full_withInteraction, "beta")$beta 
    lineageMCdat <- data.frame(Euc=rep(NA, 3000),
                              Symph=rep(NA, 3000))
   for(i in 1:3000){
     lineageMCdat$Euc[i] <- mean(betaIters[i,species$no[1:11],2] + betaIters[i,species$no[1:11],4])
      lineageMCdat$Symph[i] <- mean(betaIters[i,species$no[12:23],2] + betaIters[i,species$no[12:23],4])
    }
    (exp(quantile(lineageMCdat$Euc, probs = 0.1))-1)*100
    (exp(quantile(lineageMCdat$Symph, probs = 0.1))-1)*100
    ltys <- c(1,1)
  for(i in 1:2){
    lines(xs,getysInt(overallBetaEstsSubInt$mean[i], 
                      overallBetaEstsSubSlope1$mean[i],
                      overallBetaEstsSubSlope2$mean[i],
                      overallBetaEstsSubSlope3$mean[i]),type="l", col=cols[i], lty=ltys[i])
  }
  
    mtext(c("Subgenus:", "Symphyomyrtus", "Eucalyptus"), side=4, at=c(40,36,32),las=1,col=c("black",rev(cols)), font=c(1,3,3), cex = 0.6, line=0.5)
    
  
  ######## Lineage responses to N in low P ########
  
  overallBetaEsts$lineage <- rep(c("PA", "AlpWBYgum", "AlpWBYgum", "PA", "AlpWBYgum", "Wgum", "PA", NA, "AlpWBYgum",
                                   "PA", "AlpWBYgum", "PA", "PA", "PA", "PA", "PA", "AlpWBYgum", "Wgum", "PA",
                                   "AlpWBYgum", "PA", "AlpWBYgum", "Wgum"), 4)
  
  overallBetaEstsLinInt <- aggregate(mean ~ lineage, data=overallBetaEsts[c(1:nospp),], FUN="mean")
  overallBetaEstsLinSlope1 <- aggregate(mean ~ lineage, data=overallBetaEsts[c(1:nospp)+nospp,], FUN="mean")
  overallBetaEstsLinSlope2 <- aggregate(mean ~ lineage, data=overallBetaEsts[c(1:nospp)+nospp*2,], FUN="mean")
  overallBetaEstsLinSlope3 <- aggregate(mean ~ lineage, data=overallBetaEsts[c(1:nospp)+nospp*3,], FUN="mean")
  
  overallBetaEstsLinSlope1$meanBT <- (exp(overallBetaEstsLinSlope1$mean)-1)*100
  
  xs<-seq(0,10,length.out = 1000)
  getys<-function(intercept,slope){
    predictedy<-NA
    for(i in 1:1000){
      predictedy[i] <- exp(intercept + slope*xs[i])
    }
    return(predictedy)
  }
  
  plot(0,0, ylim=c(0,40), xlim=c(0,10.5), type="n", las=1, 
       mgp=c(2, 0.65, 0),
       xlab="", xaxt="n")
  mtext("Biomass (g)", 2, line=2.25, cex=0.9)
  abline(h=seq(0,40,5), col=scales::alpha("gray",0.4), lty=3)
  
  cols <- c("goldenrod", "purple", "tomato")
    betaIters <- extract(fit_full_withInteraction, "beta")$beta 
    lineageMCdat <- data.frame(PepAsh=rep(NA, 3000),
                               AlpWBYgum=rep(NA, 3000),
                               Wgum=rep(NA, 3000))
    for(i in 1:3000){
      lineageMCdat$PepAsh[i] <- mean(betaIters[i,species$no[1:11],2])
      lineageMCdat$AlpWBYgum[i] <- mean(betaIters[i,species$no[12:19],2])
      lineageMCdat$Wgum[i] <- mean(betaIters[i,species$no[21:23],2])
    }
    (exp(quantile(lineageMCdat$PepAsh, probs = 0.1))-1)*100
    (exp(quantile(lineageMCdat$AlpWBYgum, probs = 0.1))-1)*100
    (exp(quantile(lineageMCdat$Wgum, probs = 0.1))-1)*100
    ltys <- c(1,1,1)
  for(i in 1:3){
    lines(xs,getys(overallBetaEstsLinInt$mean[i], 
                   overallBetaEstsLinSlope1$mean[i]),type="l", col=cols[i], lty=ltys[i])
  }
  
  ######## Lineage responses to N in high P ########
  
  
  getysInt<-function(intercept,slope1,slope2,interact){
    predictedy<-NA
    for(i in 1:1000){
      predictedy[i] <- exp(intercept + slope1*xs[i] + slope2 + interact*xs[i])
    }
    return(predictedy)
  }
  
  plot(0,0, ylim=c(0,40), xlim=c(0,10.5), type="n", las=1, 
       mgp=c(2.25, 0.65, 0),
       xlab="", xaxt="n", yaxt="n")
  abline(h=seq(0,40,5), col=scales::alpha("gray",0.4), lty=3)
  
    betaIters <- extract(fit_full_withInteraction, "beta")$beta 
    lineageMCdat <- data.frame(PepAsh=rep(NA, 3000),
                              AlpWBYgum=rep(NA, 3000),
                               Wgum=rep(NA, 3000))
   for(i in 1:3000){
     lineageMCdat$PepAsh[i] <- mean(betaIters[i,species$no[1:11],2] + betaIters[i,species$no[1:11],4])
     lineageMCdat$AlpWBYgum[i] <- mean(betaIters[i,species$no[12:19],2] + betaIters[i,species$no[12:19],4])
     lineageMCdat$Wgum[i] <- mean(betaIters[i,species$no[21:23],2] + betaIters[i,species$no[21:23],4])
   }
   (exp(quantile(lineageMCdat$PepAsh, probs = 0.1))-1)*100
   (exp(quantile(lineageMCdat$AlpWBYgum, probs = 0.1))-1)*100
   (exp(quantile(lineageMCdat$Wgum, probs = 0.1))-1)*100
   ltys <- c(1,1,2)
  for(i in 1:3){
    lines(xs,getysInt(overallBetaEstsLinInt$mean[i], 
                      overallBetaEstsLinSlope1$mean[i],
                      overallBetaEstsLinSlope2$mean[i],
                      overallBetaEstsLinSlope3$mean[i]),type="l", col=cols[i], lty=ltys[i])
  }
  
   mtext(c("Lineage:","Wgum", "AlpWBYgum", "PepAsh"), side=4, at=c(40,36,32,28),las=1,col=c("black",cols[c(3,1,2)]), font=c(1,3,3,3), cex = 0.6, line=0.5)

  
  
  ######## Species responses to N in low P ########
  
  plot(0,0, ylim=c(0,40), xlim=c(0,10.5), type="n", las=1, 
       mgp=c(2.25, 0.65, 0),
       xaxt="n")
  axis(1, at=c(0, 5, 10), cex=0.7)
  axis(1, at=c(1.2, 2.4), cex=0.7)
  mtext(expression(paste("Nitrogen (g m"^-2, " yr"^-1,")")), 1, line=3, cex=0.8, at = 10.5)
  mtext("Biomass (g)", 2, line=2.25, cex=0.9)
  abline(h=seq(0,40,5), col=scales::alpha("gray",0.4), lty=3)
  cols <- distinctColorPalette(23)
    mcNit <- rep(NA, nospp)
    names(mcNit) <- species$sp[order(species$sp)]
    for(i in 1:nospp){
       mcNit[i] <- sum(betaIters[,i,2] > 0) / length(betaIters[,i,2])
    }
    mcNit <- mcNit[rev(species$no)]
    mcNit <- round(mcNit, 2)
  ltys <- rep(1,23); ltys[10] <- 2
  for(i in 1:23){
    lines(xs,getys(overallBetaEsts[species$no[i],1], overallBetaEsts[nospp+species$no[i],1]),type="l", col=cols[i], lty=ltys[i])
  }
  
  ######## Species responses to N in high P ########
  
  
  plot(0,0, ylim=c(0,40), xlim=c(0,10.5), type="n", las=1, 
       mgp=c(2.25, 0.65, 0),
       xaxt="n", yaxt="n")
  axis(1, at=c(0, 5, 10), cex=0.7)
  axis(1, at=c(1.2, 2.4), cex=0.7)
  abline(h=seq(0,40,5), col=scales::alpha("gray",0.4), lty=3)
  
    mcNit <- rep(NA, nospp)
    names(mcNit) <- species$sp[order(species$sp)]
    for(i in 1:nospp){
      mcNit[i] <- sum(betaIters[,i,2] + betaIters[,i,4] > 0) / length(betaIters[,i,2])
    }
    mcNit <- mcNit[rev(species$no)]
    mcNit <- round(mcNit, 2)
    ltys <- rep(1,23); ltys[c(2,8,9,10,17,21,22)] <- 2
    for(i in 1:23){
    lines(xs,getysInt(overallBetaEsts[species$no[i],1], 
                      overallBetaEsts[nospp+species$no[i],1], 
                      overallBetaEsts[(nospp*2)+species$no[i],1], 
                      overallBetaEsts[(nospp*3)+species$no[i],1]),type="l", col=cols[i], lty=ltys[i])
    }
    mtext(c("Species:",rev(as.character(species$sp))), side=4, at=seq(40,-10,length.out = 24),las=1,col=c("black",rev(cols)), font=c(1,rep(3,23)), cex = 0.5, line=0.5)
    

  dev.off()
