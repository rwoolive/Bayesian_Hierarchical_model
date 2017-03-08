data {
  int<lower=1> nSpp; # no. of Eucalyptus spp. (23)
  int<lower=0> KBiom; # no. of ind-lev predictors for biomass (4)
  int<lower=1> nObsBiom; # no of samples for biomass (722)
  
  int nFun; # no. of functional traits (3)
  vector[nSpp] UmatBiom; # spp. level design matrix (column of all 1's)
  int<lower=1, upper=nSpp> sppBiom[nObsBiom]; # vector of spp identities
  matrix<lower=0>[nObsBiom, KBiom] x; # indiv. level design matrix (722 x 4)
  matrix[nObsBiom,nFun] funTraits; # indiv. fun. traits (722 x 3)
    
  vector<lower=0>[nObsBiom] ObsBiom; # vec of biomass measurements (722)
  real LKJParam; # Usually a 2, but good to have in data section for flexibility
  
  corr_matrix[nSpp] phyloCor; # phylo correlation matrix
  matrix[nSpp, nSpp] d;
  
}
transformed data{
  vector[nObsBiom] lObsBiom;
  lObsBiom = log(ObsBiom); # log transform biomass
  
}
parameters {
  real<lower=0, upper=1> lambda;
  matrix[KBiom, nSpp] betaZ; # rand normal deviates for mv matt trick
  cholesky_factor_corr[KBiom] omegaL; # prior correlation		
  cholesky_factor_corr[KBiom] omegaLPhi[nFun]; # prior correlation	for Phi
  
  vector<lower=0>[KBiom] tau; 	# prior scale 
  matrix<lower=0>[nFun,KBiom] tauPhi; 	# prior scale 
  
  vector[KBiom] gamma;  # spp-level coeffs
  vector<lower=0>[nSpp] sigmaRaw;	    # prediction error scale
  matrix<lower=0>[nFun, nSpp] sigmaRawFun;	    # prediction error scale
  vector<lower=0>[nFun] hierSigmaFun;
  real<lower=0> hierSigma;
  
  matrix[KBiom, nSpp] phiZ[nFun]; // Species-level functional trait means; raw.
  matrix[nFun,KBiom] hierPhi; // hyperprior for species-level functional traits
  matrix[nFun,KBiom] delta; // coefficients for for the effect of the Phi on beta.
}


transformed parameters {
  matrix[nSpp, nSpp] pCor;
  matrix[nSpp, nSpp] lamCor;
  
  matrix[nSpp, KBiom] beta; # ind-level coeffs  
  matrix[nSpp, KBiom] phi[nFun]; # species-level functional traits  
  
  vector<lower=0>[nSpp] sigma;	    # prediction error scale
  row_vector<lower=0>[nSpp] sigmaFun[nFun];	    # prediction error scale
  
  lamCor = lambda * (phyloCor - d) + d;
  pCor = cholesky_decompose(lamCor);
  
  
  sigma = sigmaRaw * hierSigma;
  for(i in 1:nFun)
    sigmaFun[i] = sigmaRawFun[i] * hierSigmaFun[i];
  
  
  for(i in 1:nFun)
    phi[i] = UmatBiom * hierPhi[i] + 
      (diag_pre_multiply(tauPhi[i], omegaLPhi[i]) * phiZ[i])';
  
  
  beta = UmatBiom * gamma' + // gamma is the intercept
  append_col(append_col(col(phi[1],1), col(phi[2],1)), col(phi[3],1)) * delta + // delta is functional trait coefficients
  (diag_pre_multiply(tau, omegaL) * betaZ * pCor)';  // this is the error/noise
}






model{
to_vector(tauPhi) ~ cauchy(0, 2.5); # change to cauchy(0, 10) or cauchy(0, 20)
tau ~ cauchy(0, 2.5); # change to cauchy(0, 10) or cauchy(0, 20)

to_vector(betaZ) ~ normal(0,1); 
for(i in 1:nFun) 
to_vector(phiZ[i]) ~ normal(0,1);  
omegaL ~ lkj_corr_cholesky(LKJParam);
for(i in 1:nFun)  omegaLPhi[i] ~ lkj_corr_cholesky(LKJParam);
gamma ~ normal(0,1);
to_vector(hierPhi) ~ normal(0,3);
to_vector(delta) ~ normal(0,3);
sigmaRaw ~ cauchy(0, 1);
hierSigma ~ cauchy(0, 1);
to_vector(sigmaRawFun) ~ cauchy(0, 1);
hierSigmaFun ~ cauchy(0, 1);

{ vector[nObsBiom] eta; // Linear Predictor
vector[nObsBiom] err; // Error Term
matrix[nObsBiom, nFun] etaFun;
matrix[nObsBiom, nFun] errFun;
for(n in 1: nObsBiom){
eta[n] = beta[sppBiom[n]] * x[n]';
err[n] = sigma[sppBiom[n]];
}
for(i in 1:nFun) { for(n in 1: nObsBiom){
  etaFun[n,i] = phi[i,sppBiom[n]] * x[n]';
  errFun[n,i] = sigmaFun[i,sppBiom[n]];
}}
  lObsBiom ~ normal(eta, err);
  to_vector(funTraits) ~ normal(to_vector(etaFun), to_vector(errFun));
}
  
}
generated quantities {
vector[nObsBiom] log_lik;
vector[nObsBiom] yNew;
vector[nObsBiom] eta;
for(n in 1: nObsBiom) {
eta[n] = sum(beta[sppBiom[n]] .* x[n]);
log_lik[n] = normal_log(ObsBiom[n], eta[n], sigma[sppBiom[n]]);
yNew[n] = normal_rng(eta[n], sigma[sppBiom[n]]);
}
}