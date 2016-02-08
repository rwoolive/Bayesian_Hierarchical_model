data {
  int<lower=1> nSpp; # no of Eucalyptus spp.
  int<lower=0> KBiom; # no. of ind-lev predictors (and intercept) for biomass
  int<lower=0> LBiom; # no. of spp-lev predictors (and intercept) for biomass
  int<lower=1> nObsBiom; # no of samples for biomass
  
  int nFun; # number of functional traits
  vector[nSpp] UmatBiom; # spp. level design matrix
  int<lower=1, upper=nSpp> sppBiom[nObsBiom]; # vector of spp identities for biomass
  matrix<lower=0>[nObsBiom, KBiom] x; # indMatBiom; # indiv. level design matrix
  #matrix<lower=0>[nObsBiom, KBiom] x; # x-values, including intercept
  row_vector[nFun] funTraits[nObsBiom]; 
  /* All for imputing missing data
  	int nMissingFun;
  	int whereFunMissingK[nMissingFun];
  	int whereFunMissingN[nMissingFun];
  	int nMissingY;
  	int whereY;
  	int whereMissingY;
  */
  
  vector<lower=0>[nObsBiom] ObsBiom; # vec of biomass measurements
  real LKJParam; // Usually a 2, but good to have in data section for flexibility
}

transformed data{
	vector[nObsBiom] lObsBiom;
	lObsBiom <- log(ObsBiom);
}

parameters {
  matrix[KBiom, nSpp] betaZ; # rand normal deviates for mv matt trick
  cholesky_factor_corr[KBiom] omegaL; # prior correlation		
  cholesky_factor_corr[nFun] omegaLPhi; # prior correlation	for Phi
  
  vector<lower=0>[KBiom] tau; 	# prior scale 
  vector<lower=0>[nFun] tauPhi; 	# prior scale 

  vector[KBiom] gamma;  # spp-level coeffs
  vector<lower=0>[nSpp] sigmaRaw;	    # prediction error scale
  real<lower=0> hierSigma;
  matrix[nFun, nSpp] phiZ; // Species-level functional trait means; raw.
  vector[nFun] hierPhi; // hyperprior for species-level functional traits
  row_vector[nFun] delta; // coefficients for for the effect of the Phi's on beta.
  
  // vector[nMissingY] yMissing;
  // vector[nMissingFun] funMissing;
}
transformed parameters {
  matrix[nSpp, KBiom] beta; # ind-level coeffs  
  matrix[nSpp, nFun] phi; # species-level functional traits  
  vector[nSpp] betaNraw;
  vector<lower=0>[nSpp] sigma;	    # prediction error scale
  
  sigma <- sigmaRaw * hierSigma;
  
  beta <- UmatBiom * gamma' + 
      (diag_pre_multiply(tau, omegaL) * betaZ)';
  
  betaNraw <- col(beta, 2);   // This extracts the beta column for nitrogen before the phi's are taken into account.
  
  phi <- UmatBiom * hierPhi' + 
      (diag_pre_multiply(tauPhi, omegaLPhi) * phiZ)';
  #####
  
  beta <- beta + append_col(rep_vector(0, nSpp), 
    phi * delta);
  }
model{
  /* For missing data imputation
  
	  row_vector[nFun] funTraitsFull[nObsBiom]; 
	  funTraitsFull <- funTraits;
	  for(i in 1:nMissingFun){
	    
	  }
  */
  
  tauPhi ~ cauchy(0, 2.5);
  tau ~ cauchy(0, 2.5);
  to_vector(betaZ) ~ normal(0,1); 
  to_vector(phiZ) ~ normal(0,1); 
  omegaL ~ lkj_corr_cholesky(LKJParam);
  omegaLPhi ~ lkj_corr_cholesky(LKJParam);
  gamma ~ normal(0,3);
  delta ~ normal(0,3);
  sigmaRaw ~ cauchy(0, 1);
  hierSigma ~ cauchy(0, 1);
  
  { vector[nObsBiom] eta; // Linear Predictor
    vector[nObsBiom] err; // Error Term

    for(n in 1: nObsBiom){
  	  eta[n] <- beta[sppBiom[n]] * x[n]';
      err[n] <- sigma[sppBiom[n]];
      
    }
    lObsBiom ~ normal(eta, err);
    // yMissing ~ normal(eta[whereMissingY], err[whereMissingY]);
    
  }
}
generated quantities {
  vector[nObsBiom] log_lik;
  vector[nObsBiom] yNew;
  vector[nObsBiom] eta;
  for(n in 1: nObsBiom) {
    eta[n] <- sum(beta[sppBiom[n]] .* x[n]);
  	  log_lik[n] <- normal_log(ObsBiom[n], eta[n], sigma[sppBiom[n]]);
  	  yNew[n] <- normal_rng(eta[n], sigma[sppBiom[n]]);
   }
}


