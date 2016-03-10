#=======================
# THIS MODEL USES MULTI-SPECIES DATA ON INDIVIDUAL FUNCTIONAL TRAITS AND BIOMASS ACROSS NUTRIENT TREATMENTS TO ESTIMATE 1) BIOMASS AND FUNCTIONAL TRAIT INTERCEPTS AND RESPONSES TO NUTRIENTS FOR EACH SPECIES AND 2) EFFECTS OF FUNCTIONAL TRAIT INTERCEPTS ON BIOMASS INTERCEPTS AND RESPONSES TO NUTRIENTS
#=======================


library(ape)
library(shinystan)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("PATH TO WORKING DIRECOTRY")


eucal <- read.csv("EucsGH.csv", stringsAsFactors=FALSE)
eucal$n.g.m2.yr <- eucal$n.kg.ha.yr/10
eucal$p.g.m2.yr <- eucal$p.kg.ha.yr/10

tempDat <- eucal[,c("sp","dead","subgenus","n.g.m2.yr","p.g.m2.yr","sla.cm2.g","srl.cm.g","ssd.mg.mm3","total.biomass.g")]
o <- tempDat$sp!= "E. morrisbyi" 
tempDat1 <- tempDat[o,]
p <- tempDat1$dead!= "y"
tempDat2 <- tempDat1[p,]

#######################################
tree29 <- read.tree("Euc.final.dated.tre")
tree29$tip.label <- paste("E.", tree29$tip.label, sep=" ")
tree <- drop.tip(tree29, c("E. morrisbyi","E. nitida","E. vernicosa" ,"E. johnstonii","E. perriniana","E. nitens" ))
tree$tip.label[8] <- "E. tenuiramis"
phylocor <- vcv.phylo(tree,corr=TRUE)
corNames <- order(colnames(phylocor))
phyCor <- phylocor[corNames, corNames]

##########################################
#lower level biomass
tempBiom <- tempDat2[, c("sp", "n.g.m2.yr", "p.g.m2.yr", "total.biomass.g", "sla.cm2.g", "srl.cm.g", "ssd.mg.mm3")]

tempBiom2 <- na.omit(tempBiom) # 722 total individuals
naDrop = which(is.na(tempBiom[,4]))

# str(tempBiom2)
# head(tempBiom2)


Biom <- tempBiom2$total.biomass.g # Folded normal?

nitroBiom <- tempBiom2$n.g.m2.yr # nitrogen treatment
phosBiom <- tempBiom2$p.g.m2.yr  # phosphorus treatment


sppBiom <- as.integer(as.factor(tempBiom2 $sp)) # species identity as integers: in alph order
funTraits <- log(tempBiom2[,5:7])

indMatBiom <- cbind(rep(1, length(nitroBiom)), nitroBiom, phosBiom)

groupVecBiom <- rep(1, length(unique(sppBiom)))

groupMatBiom <- matrix(rep(1, length(unique(sppBiom))), nrow=length(unique(sppBiom)), ncol=1)

UmatBiom <- cbind(groupVecBiom)
UmatFun <- rbind(groupVecBiom, groupVecBiom, groupVecBiom)


biomData <- list(
            nSpp = length(unique(sppBiom)), # 23 species
            KBiom = ncol(indMatBiom), # 3 individual level predictors
            LBiom = ncol(UmatBiom), # 1 species level predictor
            nObsBiom = length(Biom), # 722 total individuals
            nFun = ncol(funTraits), # 3 functional traits
            UmatBiom = UmatBiom[,1], # species predictor matrix (biomass)
            #UmatFun = UmatFun,
            sppBiom = sppBiom, # species identities
            x = indMatBiom, # individual predictor matrix
            funTraits = funTraits, # functional trait matrix
            ObsBiom = Biom, # individual biomass
            LKJParam = 2, 
            phyloCor = phyCor,
            d = diag(diag(phyCor)))

logBiom = log(Biom)


fit_phy_withPhos <- stan(file.path(getwd(), "HBM.stan"), data= biomData, iter=1000, chains=3, seed=5, control = list(adapt_delta = .9, max_treedepth = 13), cores = 3)  
