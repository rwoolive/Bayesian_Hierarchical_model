#=======================
## Things to do:
## 1) Understand what the model code is doing
## 1.5) Model validation/sensitivity analysis/model selection
## 2) Incorporate phylo corr matrix
## 3?) Put model and data on stan board and ask about ways to deal with skew
## 3) Try to back-transform model in gen. quantities for interpretability
## 4) Imputation of missing data
## 5) Try to incorporate trait slopes to predict biomass response to N?

#=======================


library(ape)
library(shinystan)
library(rstan)
library(loo)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Bayesian_Hierarchical_model")


eucal <- read.csv("EucsGH.csv", stringsAsFactors=FALSE)
eucal$n.g.m2.yr <- eucal$n.kg.ha.yr/10
eucal$p.g.m2.yr <- eucal$p.kg.ha.yr/10

tempDat <- eucal[,c("sp","dead","subgenus","n.g.m2.yr","p.g.m2.yr","sla.cm2.g","srl.cm.g","ssd.mg.mm3","total.biomass.g")]
o <- tempDat$sp!= "E. morrisbyi" & tempDat$dead!= "y"
tempDat2 <- tempDat[o,]

#######################################
tree29 <- read.tree("Euc.final.dated.tre")
tree29$tip.label <- paste("E.", tree29$tip.label, sep=" ")
tree <- drop.tip(tree29, c("E. morrisbyi","E. nitida","E. vernicosa" ,"E. johnstonii","E. perriniana","E. nitens" ))
phylocor <- vcv.phylo(tree,corr=TRUE)
corNames <- order(colnames(phylocor))
phyCor <- phylocor[corNames, corNames]

##########################################
#lower level biomass
tempBiom <- tempDat2[, c("sp", "n.g.m2.yr", "p.g.m2.yr", "total.biomass.g", "sla.cm2.g", "srl.cm.g", "ssd.mg.mm3")]

tempBiom2 <- na.omit(tempBiom)
naDrop = which(is.na(tempBiom[,4]))

# str(tempBiom2)
# head(tempBiom2)


Biom <- tempBiom2$total.biomass.g # Folded normal?

nitroBiom <- tempBiom2$n.g.m2.yr # nitrogen treatment
phosBiom <- tempBiom2$p.g.m2.yr  # phosphorus treatment


sppBiom <- as.integer(as.factor(tempBiom2 $sp)) # species identity as integers: in alph order
funTraits <- log(tempBiom2[,5:7])

indMatBiom <- cbind(rep(1, length(nitroBiom)), nitroBiom, phosBiom)
indMatBiomNP <- cbind(rep(1, length(nitroBiom)), nitroBiom) 

groupVecBiom <- rep(1, length(unique(sppBiom)))

groupMatBiom <- matrix(rep(1, length(unique(sppBiom))), nrow=length(unique(sppBiom)), ncol=1)

UmatBiom <- cbind(groupVecBiom)


biomData <- list(
            nSpp = length(unique(sppBiom)),
            KBiom = ncol(indMatBiomNP),
            LBiom = ncol(UmatBiom),
            nObsBiom = length(Biom),  
            nFun = ncol(funTraits),
            UmatBiom = UmatBiom[,1],
            sppBiom = sppBiom,
            x = indMatBiomNP,
            funTraits = funTraits,
            ObsBiom = Biom,
            LKJParam = 2)

logBiom = log(Biom)

#fit <- stan(file="HBM_CRPNextStep.stan", data= biomData, iter=1000, chains=3, seed=5, control = list(adapt_delta = .9, max_treedepth = 13), cores = 3) 

fit <- stan(file="HBM_CRP_AddN_Normal.stan", data= biomData, iter=1000, chains=3, seed=5, control = list(adapt_delta = .9, max_treedepth = 13), cores = 3) 

#fit <- stan(fit = fit, data= biomData, iter=1000, chains=6, seed=5, control = list(adapt_delta = .9, max_treedepth = 13), cores = 3) 

#fit <- stan(file.path(getwd(), "Try 3/HBM_CRP_AddN_SN.stan"), data= biomData, iter=10, chains=1, seed=5, control = list(adapt_delta = .9, max_treedepth = 13), cores = 1) 

#initSkewParam = function(i) list(rawSkewParam = runif(length(unique(sppBiom)), -3, -.1), hierSkewParam = runif(1,-2.5, -1.5))
#fit <- stan(fit=fit, data= biomData, iter=500, chains=1, 
  seed=4, control = list(adapt_delta = .99, max_treedepth = 13), cores = 1) 


launch_shinystan(fit)






### Better SN PP check
library(sn)
etas = rstan::extract(fit, "eta")$eta
sigmas = rstan::extract(fit, "sigma")$sigma
skewParams = rstan::extract(fit, "skewParam")$skewParam
str(skewParams)
nIter = dim(etas)[1]
nObs = dim(etas)[2]
nSpp = dim(sigmas)[2]
yNew = do.call(rbind, lapply(1:nIter, function(i) 
  rsn(nObs, etas[i,], sigmas[i,sppBiom], skewParams[i,sppBiom])))
  #lapply(1:nObs, function(j) {
yFoo = function(foo) sapply(1:nIter, function(i) foo(yNew[i,]))
yFoo2 = function(foo) sapply(1:nObs, function(i) foo(yNew[,i]))

yNewMin = yFoo(min)
yNewMax = yFoo(max)
yNewSD = yFoo(sd)
yNewMean = yFoo2(mean)
yNewMed = yFoo2(median)
yNewMean2 = yFoo(mean)

plot(logBiom, yNewMean)
sapply(1:nIter, function(i) points(logBiom, yNew[i,]))
abline(0,1)

hist(yNewMin); abline(v=min(logBiom))
hist(yNewMax); abline(v=max(logBiom))
hist(yNewSD); abline(v=sd(logBiom))
hist(yNewMean); abline(v=mean(logBiom))


library(dplyr)
df = data_frame(logBiom, sppBiom)
df %>% group_by(sppBiom) %>% summarise(skew = skewness(logBiom)) %$% skew %>% hist

