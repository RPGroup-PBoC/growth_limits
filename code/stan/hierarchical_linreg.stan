data {
    // Dimensional Parameters
    int<lower=1> J;  // Number of unique datasets
    int<lower=1> N; // Number of total data points
    real<lower=1, upper=J> idx[N]; // ID vector for measurements

    // Prior specifications
    real<lower=0> hyper_m_mu_val;
    real<lower=0> hyper_m_sigma_val;
    real<lower=0> hyper_ms_mu_val;
    real<lower=0> hyper_ms_sigma_val;
    real<lower=0> hyper_b_mu_val;
    real<lower=0> hyper_b_sigma_val;
    real<lower=0> hyper_bs_mu_val;
    real<lower=0> hyper_bs_sigma_val;

    // Observed data
    vector<lower=1>[N] x; 
    vector<lower=1>[N] y;
}

parameters {
    // Level-0 parameters
    real<lower=0> hyper_m;
    real<lower=0> hyper_m_sigma;
    real<lower=0> hyper_b;
    real<lower=0> hyper_b_sigma;

    // Level-1 parameters
    vector<lower=0>[J] m; 
    vector<lower=0>[J] b;

    // Homoscedastic error 
    vector<lower=0>[J] sigma;
}

model {
    vector[N] mu;
    // Level-0 prior
    hyper_m ~ normal(hyper_m_mu_val, hyper_m_sigma_val);
    hyper_m_sigma ~ normal(hyper_ms_mu_val, hyper_ms_sigma_val);
    hyper_b ~ normal(hyper_b_mu_val, hyper_b_sigma_val);
    hyper_b_sigma ~ normal(hyper_bs_mu_val, hyper_bs_sigma_val);

    // Level-1 priors
    m ~ normal(hyper_m, hyper_m_sigma);
    b ~ normal(hyper_b, hyper_b_sigma);
    
    // Homoscedastic error 
    sigma ~ normal(0, 1);

   // Define the likelihood 
   y ~ normal(m[idx] * x + b[idx], sigma[idx]);
}