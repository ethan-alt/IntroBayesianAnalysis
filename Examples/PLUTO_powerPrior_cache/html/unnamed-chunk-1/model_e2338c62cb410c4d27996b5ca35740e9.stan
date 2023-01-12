data {
  int<lower=0> ntrt;                 // num. treated patients in current data set
  int<lower=0> nctrl;                // num. placebo patients in current data set
  int<lower=0,upper=ntrt> ytrt;      // num. responders in treated group
  int<lower=0,upper=nctrl> yctrl;    // num. responders in control group
  int<lower=0> J;                    // num. of historical data sets
  array[J] int<lower=0> ntrt0;       // num. treated patients in historical data sets
  array[J] int<lower=0> nctrl0;      // num. placebo patients in historical data sets
  array[J] int<lower=0> ytrt0;       // num. responders in treated group for historical data sets
  array[J] int<lower=0> yctrl0;      // num. responders in placebo group for historical data sets
  array[J] real<lower=0,upper=1> a0; // discounting parameter (could also be vector type)
  real<lower=0> theta_shape1;        // first shape hyperparameter for initial beta prior
  real<lower=0> theta_shape2;        // second shape hyperparameter for initial beta prior
}
parameters {
  real<lower=0,upper=1> theta1;  // probability of responding among treated
  real<lower=0,upper=1> theta0;  // probability of responding among placebo
}
model {
  // power priors
  for ( j in 1:J ) {
    target += a0[j] * binomial_lpmf(ytrt0[j] | ntrt0[j], theta1);    // power prior for treated in j-th data set
    target += a0[j] * binomial_lpmf(yctrl0[j] | nctrl0[j], theta0);  // power prior for placebo in j-th data set
  }
  theta1 ~ beta(theta_shape1, theta_shape2);   // initial prior
  theta0 ~ beta(theta_shape1, theta_shape2);   // initial prior
  
  // likelihood
  ytrt ~ binomial(ntrt, theta1);
  yctrl ~ binomial(nctrl, theta0);
}
generated quantities {
  real risk_difference = theta1 - theta0;
  real risk_ratio = theta1 * inv(theta1);  // = theta1 / theta0
  real odds_ratio = (theta1 * inv(1-theta1)) / (theta0 * inv(1-theta0));
}
