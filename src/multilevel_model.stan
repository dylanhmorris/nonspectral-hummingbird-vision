/* name: multilevel_model.stan
   author: Dylan Morris <dhmorris@princeton.edu>

   description: Stan code for a multilevel
   model of hummingbird vision in which estimates
   for color discrimination are partially pooled
   for experiments in the same treatment type
   (spectral, nonspectral, null/control)
*/


functions{

  vector full_to_half(vector x){
    return(x/2 + 0.5);
  }
    
  real half_to_full(real x){
    return((x - 0.5) * 2);
  }
  
  real indef_integral_E(real equilibrium_E,
                        real saturation_rate,
                        real t) {
    return(equilibrium_E *
           ((exp(-saturation_rate * t) / saturation_rate) + t));
  }

  real average_E(real equilibrium_E,
                 real saturation_rate,
                 real t_initial,
                 real t_final){
    if(t_final > t_initial){
      return((indef_integral_E(equilibrium_E, saturation_rate, t_final) -
              indef_integral_E(equilibrium_E, saturation_rate, t_initial)) /
             (t_final - t_initial));
    }
    else{
      return(-equilibrium_E * exp(-saturation_rate * t_initial) +
             equilibrium_E);
    }
  }

  real bdist_alpha(real mode,
                   real certainty){
    return(mode * (certainty - 2) + 1);
  }
  
  real bdist_beta(real mode,
                  real certainty){
    return((1 - mode) * (certainty - 2) + 1);
  }

}

data {
  
  int<lower=0> n_experiments;
  int<lower=0> n_trials_total;
  int<lower=0> n_locations;
  int<lower=0> n_treatments;
  int<lower=0> cumulative_hits_at_trial_start[n_trials_total];
  int<lower=0> n_successes[n_trials_total];
  int<lower=0> n_failures[n_trials_total];
  
  int<lower=-1, upper=1> correct_side[n_trials_total]; 
  // 1 = left
  // -1 = right
  
  int<lower=1, upper=n_experiments> experiment_id[n_trials_total];
  int<lower=1, upper=n_locations> location[n_experiments];
  int<lower=1, upper=n_treatments> treatment_id[n_experiments];

  real min_learning_effect;

  // hyperparameters set at runtime; they should be chosen
  // to be weakly informative



  // for equilibrium fraction experienced,
  // real number between 0 and 1 ~ Beta(alpha, beta)
  real hyper_alpha_equilibrium_E;
  real hyper_beta_equilibrium_E;

  // for half-saturation time of approach to equilibrium E
  // positive-constrained Gaussian 
  real<lower=0> hyper_mean_experience_half_saturation_time;
  real<lower=0> hyper_sd_experience_half_saturation_time;

  // for the mode correctness in each treatment
  real<lower=inv_logit(min_learning_effect), upper=1> hyper_mode_treatment_mode;
  real<lower=0> hyper_certainty_treatment_mode;

  
  // for error sds-- key to make sure
  // sd does not exceed 1.5
  // as this corresponds to uniform
  // in probability-space, and
  // greater sds push most mass to the
  // endpoints (0 and 1)
  real hyper_mean_sd_error;
  real hyper_sd_sd_error;

  real<lower=0> hyper_mean_sd_left_bias;
  real<lower=0> hyper_sd_sd_left_bias;
  
  real<lower=0> hyper_mean_treatment_certainty;
  
  // for each location's mean left bias (weakly informative)
  real<lower=0, upper=1> hyper_mode_loc_left_bias;
  real<lower=2> hyper_certainty_loc_left_bias;
    
}

transformed data{
  int<lower=0> n_hits[n_trials_total];
  int<lower=0> n_hits_pred[n_trials_total];
  real<lower=0, upper=1> min_prob;
  
  for(trial_id in 1:n_trials_total){
    n_hits[trial_id] = n_failures[trial_id] + n_successes[trial_id];
    n_hits_pred[trial_id] = n_hits[trial_id];
  }

  min_prob = inv_logit(min_learning_effect);
}


parameters{
  vector[n_trials_total] error_correct;
  vector[n_experiments] error_bias;

  real<lower=0> sd_error;
  real<lower=0> sd_left_bias;
  
  vector<lower=0>[n_treatments] treatment_certainty;

  vector<lower=0, upper=1>[n_locations] loc_p_left_naive;

  real<lower=0, upper=1> equilibrium_E;
  real<lower=0> experience_half_saturation_time;

  vector<lower=min_prob, upper=1>[n_treatments] treatment_mode;
  vector<lower=0, upper=1>[n_experiments] interval_learning;

}

transformed parameters{
  vector<lower=min_learning_effect>[n_experiments] learning_effect;
  vector<lower=min_prob, upper=1>[n_experiments] unbiased_p_correct_experienced;
  real<lower=0> experience_saturation_rate;

  
  vector<lower=0, upper=1>[n_trials_total] p_correct_naive;
  vector<lower=0, upper=1>[n_trials_total] p_correct_experienced;
  vector<lower=0, upper=1>[n_trials_total] E;
  vector<lower=0, upper=1>[n_trials_total] p_correct_mean;

  vector<lower=0, upper=1>[n_trials_total] p_correct_overall;

  vector[n_locations] loc_mode_left_bias;
  vector<lower=0, upper=1>[n_experiments] p_left_naive;
  vector[n_experiments] left_bias;

  loc_mode_left_bias = logit(loc_p_left_naive);
  
  for(exp_id in 1:n_experiments)
    left_bias[exp_id] = loc_mode_left_bias[location[exp_id]] +
      error_bias[exp_id] *
      sd_left_bias;
  p_left_naive = inv_logit(left_bias);

  unbiased_p_correct_experienced = full_to_half(interval_learning);
  learning_effect = logit(unbiased_p_correct_experienced);
 
  experience_saturation_rate = log(2) / experience_half_saturation_time;

  for(trial_id in 1:n_trials_total){
    int exp_id = experiment_id[trial_id];

    // assuming E(0) ~= 0 for the moment
    // and dE/dt = becoming_experienced_rate * (1 - E) - turnover_rate * E
      
    p_correct_experienced[trial_id] =
      inv_logit(left_bias[exp_id]
                * correct_side[trial_id]
                + learning_effect[exp_id]);
    
    p_correct_naive[trial_id] =
      inv_logit(left_bias[exp_id] * correct_side[trial_id]);

    E[trial_id] = average_E(equilibrium_E,
                            experience_saturation_rate,
                            cumulative_hits_at_trial_start[trial_id],
                            cumulative_hits_at_trial_start[trial_id] +
                            n_hits[trial_id]);
    
    p_correct_mean[trial_id] =
      E[trial_id] * p_correct_experienced[trial_id] +
      (1 - E[trial_id]) * p_correct_naive[trial_id];
    
    p_correct_overall[trial_id] = 1 / (1 + (exp(-error_correct[trial_id] * sd_error) * (1 - p_correct_mean[trial_id]) / p_correct_mean[trial_id]));
  }
}

model {

  // binomial count data
  n_successes ~ binomial(n_hits,
                         p_correct_overall);
  
  // non-centered. sd will be multiplied
  error_correct ~ normal(0, 1);
  error_bias ~ normal(0, 1);

  treatment_certainty ~ exponential(1 / hyper_mean_treatment_certainty);
  
  // prior on the treatment mode correctnesses
  treatment_mode ~
  beta(bdist_alpha(hyper_mode_treatment_mode,
                  hyper_certainty_treatment_mode),
       bdist_beta(hyper_mode_treatment_mode,
                  hyper_certainty_treatment_mode));
       
  
  // weakly informative sd priors
  // sds should be kept below 1.5
  // in logit space, as this
  // corresponds to uniform in
  // probability space
  
  sd_error ~ normal(hyper_mean_sd_error,
                    hyper_sd_sd_error);
  sd_left_bias ~ normal(hyper_mean_sd_left_bias,
                        hyper_sd_sd_left_bias);

  // weakly informative priors for equilibrium E
  // and half saturation of E
  equilibrium_E ~ beta(hyper_alpha_equilibrium_E,
                       hyper_beta_equilibrium_E);
  
  experience_half_saturation_time ~
     normal(hyper_mean_experience_half_saturation_time,
            hyper_sd_experience_half_saturation_time);

  for (exp_id in 1:n_experiments){
    int t_id = treatment_id[exp_id];
    real full_mode = half_to_full(treatment_mode[t_id]);
    interval_learning[exp_id] ~
      beta(bdist_alpha(full_mode,
                       2 + treatment_certainty[t_id]),
           bdist_beta(full_mode,
                      2 + treatment_certainty[t_id]));
  }
  
  loc_p_left_naive ~ beta(bdist_alpha(hyper_mode_loc_left_bias,
                                      hyper_certainty_loc_left_bias),
                          bdist_beta(hyper_mode_loc_left_bias,
                                     hyper_certainty_loc_left_bias));

}

generated quantities {}
