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
//  row_vector[nFun] funTraits[nObsBiom]; 
  matrix[nObsBiom,nFun] funTraits; 

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
  cholesky_factor_corr[2] omegaLPhi[nFun]; # prior correlation	for Phi
  //cholesky_factor_corr[nFun] omegaLXi; # prior correlation	for Phi
  
  vector<lower=0>[KBiom] tau; 	# prior scale 
  matrix<lower=0>[nFun,2] tauPhi; 	# prior scale 

  vector[KBiom] gamma;  # spp-level coeffs
  vector<lower=0>[nSpp] sigmaRaw;	    # prediction error scale
  matrix<lower=0>[nFun, nSpp] sigmaRawFun;	    # prediction error scale
  vector<lower=0>[nFun] hierSigmaFun;
  real<lower=0> hierSigma;

  matrix[2, nSpp] phiZ[nSpp]; // Species-level functional trait means; raw.
  matrix[nFun,2] hierPhi; // hyperprior for species-level functional traits
  matrix[nFun,KBiom] delta; // coefficients for for the effect of the Phi's on beta.

  // Add The New Stuff...
  // Essentially, add in a set of slopes for each of the funTraits.
  // For now, duplicate the current approach, and assume no correlation between phi and xi. 
  // vector[nMissingY] yMissing;
  // vector[nMissingFun] funMissing;
}
transformed parameters {
  matrix[nSpp, KBiom] beta; # ind-level coeffs  
  matrix[nSpp, 2] phi[nFun]; # species-level functional traits  

  vector<lower=0>[nSpp] sigma;	    # prediction error scale
  row_vector<lower=0>[nSpp] sigmaFun[nFun];	    # prediction error scale

  sigma <- sigmaRaw * hierSigma;
  for(i in 1:nFun)
    sigmaFun[i] <- sigmaRawFun[i] * hierSigmaFun[i];

  for(i in 1:nFun)
    phi[i] <- UmatBiom * hierPhi[i] + 
      (diag_pre_multiply(tauPhi[i], omegaLPhi[i]) * phiZ[i])';
  beta <- UmatBiom * gamma' + // gamma is the intercept
  append_col(append_col(col(phi[1],1), col(phi[2],1)), col(phi[3],1)) * delta + // delta is functional trait coefficients
    (diag_pre_multiply(tau, omegaL) * betaZ)';  // this is the error/noise
  }
model{
  /* For missing data imputation
  
	  row_vector[nFun] funTraitsFull[nObsBiom]; 
	  funTraitsFull <- funTraits;
	  for(i in 1:nMissingFun){
	    
	  }
  */
#  skewParam ~ student_t(7,0, 2.5);
  to_vector(tauPhi) ~ cauchy(0, 2.5);
  tau ~ cauchy(0, 2.5);
  to_vector(betaZ) ~ normal(0,1); 
  for(i in 1:nSpp) 
    to_vector(phiZ[i]) ~ normal(0,1);  
  omegaL ~ lkj_corr_cholesky(LKJParam);
for(i in 1:nSpp)  omegaLPhi[i] ~ lkj_corr_cholesky(LKJParam);
  gamma ~ normal(0,3);
  to_vector(hierPhi) ~ normal(0,3);
  to_vector(delta) ~ normal(0,3);
  sigmaRaw ~ cauchy(0, 1);
  hierSigma ~ cauchy(0, 1);
  sigmaRaw ~ cauchy(0, 1);
  hierSigma ~ cauchy(0, 1);
  to_vector(sigmaRawFun) ~ cauchy(0, 1);
  hierSigmaFun ~ cauchy(0, 1);

  { vector[nObsBiom] eta; // Linear Predictor
    vector[nObsBiom] err; // Error Term
    matrix[nObsBiom, nFun] etaFun;
    matrix[nObsBiom, nFun] errFun;
    for(n in 1: nObsBiom){
  	  eta[n] <- beta[sppBiom[n]] * x[n]';
      err[n] <- sigma[sppBiom[n]];
    }
    for(i in 1:nFun) { for(n in 1: nObsBiom){
  	  etaFun[n,i] <- phi[i,sppBiom[n]] * x[n]';
      errFun[n,i] <- sigmaFun[i,sppBiom[n]];
    }}
    #lObsBiom ~ normal(eta, err);
// We're going to try a skew-normal with this, because log(Biom) is defniitely skewed.
    lObsBiom ~ normal(eta, err, skewParam[sppBiom]);
    to_vector(funTraits) ~ normal(to_vector(etaFun), to_vector(errFun));
    }
    // yMissing ~ normal(eta[whereMissingY], err[whereMissingY]);
    
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