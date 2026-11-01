// ============================================================================
// Bayesian false-positive occupancy model, multiple-methods (site-confirmation)
// design of Miller et al. (2011): an automated method that can return false
// positives, plus verified detections that cannot.
//
//   z_i        ~ Bernoulli(psi_i),        logit(psi_i)  = X_i beta
//   y2_ik      ~ Bernoulli(z_i p11_ik + (1 - z_i) p10_i)
//                logit(p11_ik) = Z_ik alpha + u_rec[ik] [+ d_i]
//                logit(p10_i)  = (alpha_1 - fp_gap) + Xfp_i gamma_fp [+ e_i]
//   y1_im      ~ Bernoulli(z_i r11_im),   logit(r11_im) = delta0 + X1s_im delta_s
//
// The latent z is marginalized. fp_gap > 0 keeps p11 > p10 at the reference
// (all covariates at 0), which resolves label switching. Column 1 of Z must be
// the intercept. Kfp = 0 gives the constant-rate model; estimate_fp = 0 gives
// the standard occupancy model (no false positives). estimate_fp_re = 1 adds a
// site random effect e_i ~ Normal(0, sigma_fp) to the logit of p10, and
// estimate_det_re = 1 one d_i ~ Normal(0, sigma_det) to the logit of p11
// (of the automated method; r11 has none). Both can be on at once.
// fp_by_row = 1 takes the false-positive covariates from each automated row
// (Xfp2, and Xfp3 for confirmed occasions) instead of from the site (Xfp),
// so p10 varies within sites. Kver_s > 0 gives the verified method
// covariates, with its rows sorted by site. Runs of detections (fp_occu(
// markov = )) need no code here: the previous detection is a column of Z.
//
// Observation-confirmation design (Chambert et al. 2015), when obs_confirm = 1:
// the observations of some automated occasions are checked, giving a state
// v3 = 0 (nothing), 1 (false observations only), 2 (true only), 3 (both). The
// detection and false-positive predictors are then those of s1 and s0, the
// probabilities of at least one true and at least one false observation on an
// occasion, and p11 = 1 - (1 - s1)(1 - s0), p10 = s0 for unchecked occasions.
//
// Automated occasions enter as rows of identical covariates within a site
// (n2 trials, y2 detections), which is exact for a product of Bernoulli
// trials and much faster when detection covariates are constant within
// sites. Rows must be sorted by site (fp_occu() collapses and sorts them), so
// that each site's rows are contiguous. The unoccupied branch only needs the totals
// of occasions and detections of each site, because the false-positive
// rate is constant within a site. The occupied branch uses one
// binomial_logit_glm_lpmf() call per site when detection has covariates,
// which evaluates the linear predictor, the log-likelihood and its gradient
// without intermediate autodiff variables; otherwise (rows already
// collapsed, or the observation-confirmation design) the rows are
// evaluated as vectors and summed with segment().
//
// The data of the likelihood travel in tuples (built at the end of
// transformed data) through site_loglik() and unit_emissions() to
// emissions_sites(); tuples as arguments of reduce_sum() need CmdStan 2.40
// or later. The emissions (the likelihood of each site's observations given
// that it is occupied and given that it is not) are in
// fp_emissions.stanfunctions, shared with fp_occu_dyn.stan.
//
// The sum over sites runs with reduce_sum(), which splits the sites among
// threads when the model is compiled with stan_threads (threads_per_chain in
// fp_occu()). Each slice of sites only touches its own rows of the automated
// and confirmed occasions. With one thread the result is the same sum.
// ============================================================================
functions {
#include fp_emissions.stanfunctions

  // Marginal log-likelihood of sites a to b, with the latent state z summed
  // out, from the emissions of the shared functions (fp_emissions.stanfunctions).
  vector site_loglik(int a, int b, vector beta, vector alpha, real gamma0,
                     vector gamma_fp, vector e_fp, vector u, vector d_det, real delta0,
                     vector delta_s,
                     tuple(matrix, matrix, matrix, matrix, matrix) design,
                     tuple(array[] int, array[] int, array[] int, array[] int,
                           vector, vector, array[] int, array[] int,
                           array[] int, array[] int, vector, vector, vector) automated,
                     tuple(vector, vector, array[] int, matrix, vector,
                           array[] int, array[] int) verified,
                     tuple(array[] int, array[] int, array[] int, matrix,
                           array[] int, array[] int, matrix) confirmed,
                     array[] int options) {
    int ns = b - a + 1;
    tuple(vector[ns], vector[ns]) em = unit_emissions(a, b, alpha, gamma0, gamma_fp, e_fp, u,
                                                      d_det, delta0, delta_s, design,
                                                      automated, verified, confirmed, options);
    vector[ns] lpsi = design.1[a:b] * beta;
    vector[ns] ll;
    for (i in 1:ns) {
      real occ_branch = log_inv_logit(lpsi[i]) + em.1[i];
      ll[i] = verified.3[a + i - 1] == 1
              ? occ_branch
              : log_sum_exp(occ_branch, log1m_inv_logit(lpsi[i]) + em.2[i]);
    }
    return ll;
  }

  // partial sum over a slice of sites, for reduce_sum()
  real partial_loglik(array[] int slice, int a, int b, vector beta, vector alpha,
                      real gamma0, vector gamma_fp, vector e_fp, vector u, vector d_det,
                      real delta0, vector delta_s,
                      tuple(matrix, matrix, matrix, matrix, matrix) design,
                      tuple(array[] int, array[] int, array[] int, array[] int,
                            vector, vector, array[] int, array[] int,
                            array[] int, array[] int, vector, vector, vector) automated,
                      tuple(vector, vector, array[] int, matrix, vector,
                            array[] int, array[] int) verified,
                      tuple(array[] int, array[] int, array[] int, matrix,
                            array[] int, array[] int, matrix) confirmed,
                      array[] int options) {
    return sum(site_loglik(a, b, beta, alpha, gamma0, gamma_fp, e_fp, u, d_det, delta0,
                           delta_s, design, automated, verified, confirmed, options));
  }
}

data {
  int<lower=1> nsite;
  int<lower=1> Kbeta;
  matrix[nsite, Kbeta] X;                    // occupancy design (with intercept)

  int<lower=0> N2;                           // rows of automated occasions
  array[N2] int<lower=1> n2;                 // occasions in each row
  array[N2] int<lower=0> y2;                 // detections in each row
  array[N2] int<lower=1, upper=nsite> site2;  // sorted by site
  array[nsite] int<lower=1> start2;          // first row of each site in y2
  array[nsite] int<lower=0> len2;            // number of rows of each site
  int<lower=1> Kalpha;
  matrix[N2, Kalpha] Z;                      // detection design (column 1 = intercept)
  int<lower=1> nrec;
  array[N2] int<lower=1, upper=nrec> rec2;
  int<lower=0, upper=1> estimate_re;         // recorder random intercept on p11

  int<lower=0> Kfp;
  int<lower=0, upper=1> fp_by_row;           // false-positive covariates vary by occasion
  matrix[nsite, fp_by_row == 1 ? 0 : Kfp] Xfp;  // false-positive design (no intercept)
  matrix[N2, fp_by_row == 1 ? Kfp : 0] Xfp2;    // ... by automated row
  int<lower=0, upper=1> estimate_fp;

  int<lower=0> N1;                           // verified occasions (long format)
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=nsite> site1;
  int<lower=0> Kver_s;                       // slopes of the verified detection
  matrix[N1, Kver_s] X1s;                    // their design (verified rows sorted by site)
  array[nsite] int<lower=0, upper=1> known_present;

  int<lower=0, upper=1> obs_confirm;         // observation-confirmation design
  int<lower=0> N3;                           // confirmed occasions
  array[N3] int<lower=0, upper=3> v3;        // confirmation state
  array[N3] int<lower=1, upper=nsite> site3;
  matrix[N3, Kalpha] Z3;                     // detection design of confirmed occasions
  array[N3] int<lower=1, upper=nrec> rec3;
  matrix[N3, fp_by_row == 1 ? Kfp : 0] Xfp3; // false-positive design of confirmed occasions

  real<lower=0> prior_sd_int;                // sd of Normal priors on intercepts
  real<lower=0> prior_sd_slope;              // sd of Normal priors on slopes
  real gap_mean;                             // Normal(gap_mean, gap_sd) prior on fp_gap, truncated at 0
  real<lower=0> gap_sd;
  real<lower=0> recorder_sd;                 // half-Normal(0, recorder_sd) on sigma_u

  int<lower=0, upper=1> estimate_fp_re;      // site random effect on logit(p10)
  real<lower=0> fp_site_sd;                  // half-Normal(0, fp_site_sd) on sigma_fp
  int<lower=0, upper=1> estimate_det_re;     // site random effect on logit(p11)
  real<lower=0> det_site_sd;                 // half-Normal(0, det_site_sd) on sigma_det
}

transformed data {
  vector[N2] y2v = to_vector(y2);
  vector[N2] n2v = to_vector(n2);
  // sites with an automated detection cannot be unoccupied without false positives
  array[nsite] int any_det2 = rep_array(0, nsite);
  for (n in 1:N2) {
    if (y2[n] > n2[n]) reject("y2 cannot exceed n2");
    if (n > 1 && site2[n] < site2[n - 1]) reject("automated rows must be sorted by site");
    if (y2[n] > 0) any_det2[site2[n]] = 1;
  }
  if (obs_confirm == 1 && estimate_fp == 0) {
    reject("the observation-confirmation design needs a false-positive process (estimate_fp = 1)");
  }
  // verified occasions and detections per site
  vector[nsite] n1v = rep_vector(0, nsite);
  vector[nsite] s1v = rep_vector(0, nsite);
  for (n in 1:N1) {
    n1v[site1[n]] += 1;
    s1v[site1[n]] += y1[n];
  }
  // with covariates on the verified detection, its rows are used one by one
  vector[N1] y1v = to_vector(y1);
  array[nsite] int len1 = rep_array(0, nsite);
  for (n in 1:N1) {
    if (Kver_s > 0 && n > 1 && site1[n] < site1[n - 1]) reject("verified occasions must be sorted by site");
    len1[site1[n]] += 1;
  }
  array[nsite] int start1;
  start1[1] = 1;
  for (i in 2:nsite) start1[i] = start1[i - 1] + len1[i - 1];
  // first row and number of rows of each site among the confirmed occasions
  array[nsite] int len3 = rep_array(0, nsite);
  for (n in 1:N3) {
    if (n > 1 && site3[n] < site3[n - 1]) reject("confirmed occasions must be sorted by site");
    len3[site3[n]] += 1;
  }
  array[nsite] int start3;
  start3[1] = 1;
  for (i in 2:nsite) start3[i] = start3[i - 1] + len3[i - 1];
  // sites, sliced among threads by reduce_sum()
  array[nsite] int site_index = linspaced_int_array(nsite, 1, nsite);
  // site effects on the false-positive and true-detection rates when there
  // are none (see model)
  vector[nsite] zero_site = rep_vector(0, nsite);
  // totals of automated occasions and detections per site (unoccupied branch)
  vector[nsite] Ysite = rep_vector(0, nsite);
  vector[nsite] Nsite = rep_vector(0, nsite);
  for (n in 1:N2) {
    Ysite[site2[n]] += y2[n];
    Nsite[site2[n]] += n2[n];
  }
  // the occupied branch uses one GLM call per site when detection has
  // covariates (otherwise the rows are already collapsed) and the design is
  // the multiple-methods one; Z_slopes is Z without its intercept column, and lchoose_site
  // the binomial constants that the GLM includes and the likelihood omits
  int use_glm = (Kalpha > 1 && obs_confirm == 0) ? 1 : 0;
  matrix[N2, Kalpha - 1] Z_slopes = rep_matrix(0, N2, Kalpha - 1);
  if (Kalpha > 1) Z_slopes = block(Z, 1, 2, N2, Kalpha - 1);
  vector[nsite] lchoose_site = rep_vector(0, nsite);
  for (n in 1:N2) lchoose_site[site2[n]] += lchoose(n2[n], y2[n]);
  // recorder of each site, when it does not change within sites
  int rec_by_site = 1;
  array[nsite] int rec_site = rep_array(1, nsite);
  for (n in 1:N2) {
    if (n == 1 || site2[n] != site2[n - 1]) rec_site[site2[n]] = rec2[n];
    else if (rec2[n] != rec_site[site2[n]]) rec_by_site = 0;
  }
  // the data of the likelihood, grouped for site_loglik() and reduce_sum();
  // the order of the fields is the order of the unpacking in site_loglik()
  tuple(matrix[nsite, Kbeta], matrix[N2, Kalpha], matrix[N2, Kalpha - 1],
        matrix[nsite, fp_by_row == 1 ? 0 : Kfp], matrix[N2, fp_by_row == 1 ? Kfp : 0])
    design = (X, Z, Z_slopes, Xfp, Xfp2);
  tuple(array[N2] int, array[N2] int, array[N2] int, array[N2] int, vector[N2],
        vector[N2], array[nsite] int, array[nsite] int, array[nsite] int,
        array[nsite] int, vector[nsite], vector[nsite], vector[nsite])
    automated = (y2, n2, site2, rec2, y2v, n2v, start2, len2, any_det2, rec_site,
                 Ysite, Nsite, lchoose_site);
  tuple(vector[nsite], vector[nsite], array[nsite] int, matrix[N1, Kver_s],
        vector[N1], array[nsite] int, array[nsite] int)
    verified = (n1v, s1v, known_present, X1s, y1v, start1, len1);
  tuple(array[N3] int, array[N3] int, array[N3] int, matrix[N3, Kalpha],
        array[nsite] int, array[nsite] int, matrix[N3, fp_by_row == 1 ? Kfp : 0])
    confirmed = (v3, site3, rec3, Z3, start3, len3, Xfp3);
  array[8] int options = {estimate_re, estimate_fp, Kfp, obs_confirm, use_glm,
                          rec_by_site, Kver_s, fp_by_row};
}

parameters {
  vector[Kbeta] beta;
  vector[Kalpha] alpha;
  vector[nrec] u_raw;
  real<lower=0> sigma_u;
  real<lower=0> fp_gap;
  vector[Kfp] gamma_fp;
  real delta0;
  vector[Kver_s] delta_s;
  vector[estimate_fp_re == 1 ? nsite : 0] e_raw;
  array[estimate_fp_re] real<lower=0> sigma_fp;
  vector[estimate_det_re == 1 ? nsite : 0] d_raw;
  array[estimate_det_re] real<lower=0> sigma_det;
}

transformed parameters {
  vector[nrec] u = sigma_u * u_raw;
  real gamma0 = alpha[1] - fp_gap;
}

model {
  beta[1]        ~ normal(0, prior_sd_int);
  if (Kbeta > 1)  beta[2:Kbeta]   ~ normal(0, prior_sd_slope);
  alpha[1]       ~ normal(0, prior_sd_int);
  if (Kalpha > 1) alpha[2:Kalpha] ~ normal(0, prior_sd_slope);
  u_raw          ~ std_normal();
  sigma_u        ~ normal(0, recorder_sd);
  fp_gap         ~ normal(gap_mean, gap_sd);
  gamma_fp       ~ normal(0, prior_sd_slope);
  delta0         ~ normal(0, prior_sd_int);
  delta_s        ~ normal(0, prior_sd_slope);
  e_raw          ~ std_normal();
  sigma_fp       ~ normal(0, fp_site_sd);
  d_raw          ~ std_normal();
  sigma_det      ~ normal(0, det_site_sd);

  // The four calls differ only in the site effects on the false-positive
  // and true-detection rates: parameters when estimated, otherwise the data
  // vector zero_site, which reduce_sum() neither copies nor differentiates
  // for each slice.
  if (estimate_fp_re == 1 && estimate_det_re == 1) {
    vector[nsite] e_fp = sigma_fp[1] * e_raw;
    vector[nsite] d_det = sigma_det[1] * d_raw;
    target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                         e_fp, u, d_det, delta0, delta_s, design, automated, verified,
                         confirmed, options);
  } else if (estimate_fp_re == 1) {
    vector[nsite] e_fp = sigma_fp[1] * e_raw;
    target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                         e_fp, u, zero_site, delta0, delta_s, design, automated, verified,
                         confirmed, options);
  } else if (estimate_det_re == 1) {
    vector[nsite] d_det = sigma_det[1] * d_raw;
    target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                         zero_site, u, d_det, delta0, delta_s, design, automated, verified,
                         confirmed, options);
  } else {
    target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                         zero_site, u, zero_site, delta0, delta_s, design, automated,
                         verified, confirmed, options);
  }
}

generated quantities {
  vector[nsite] log_lik;
  real p11_ref = obs_confirm == 1
                 ? 1 - (1 - inv_logit(alpha[1])) * (1 - inv_logit(gamma0))
                 : inv_logit(alpha[1]);
  real p10_ref = estimate_fp == 1 ? inv_logit(gamma0) : 0.0;
  real p10_bar = 0;
  real r11 = inv_logit(delta0);
  {
    vector[nsite] e_fp = rep_vector(0, nsite);
    if (estimate_fp_re == 1) e_fp = sigma_fp[1] * e_raw;
    vector[nsite] d_det = rep_vector(0, nsite);
    if (estimate_det_re == 1) d_det = sigma_det[1] * d_raw;
    vector[nsite] lp10 = rep_vector(gamma0, nsite);
    if (Kfp > 0 && fp_by_row == 0) lp10 += Xfp * gamma_fp;
    lp10 += e_fp;
    if (estimate_fp == 1 && fp_by_row == 0) p10_bar = mean(inv_logit(lp10));
    // with occasion-level covariates: mean rate over the automated occasions
    if (estimate_fp == 1 && fp_by_row == 1 && N2 > 0) {
      p10_bar = dot_product(n2v, inv_logit(lp10[site2] + Xfp2 * gamma_fp)) / sum(n2v);
    }
    log_lik = site_loglik(1, nsite, beta, alpha, gamma0, gamma_fp, e_fp, u, d_det, delta0,
                          delta_s, design, automated, verified, confirmed, options);
  }
}
