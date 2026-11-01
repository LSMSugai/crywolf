// ============================================================================
// Dynamic (multi-season) false-positive occupancy model.
//
//   z_i1       ~ Bernoulli(psi1_i),               logit(psi1_i)  = X_i beta
//   z_it       ~ Bernoulli(z_i,t-1 (1 - eps_it) + (1 - z_i,t-1) gam_it), t > 1
//                logit(gam_it) = Xcol_it beta_col, logit(eps_it) = Xext_it beta_ext
//
// Within each season the observation model is that of fp_occu.stan: each
// site-season (unit) is a site of that program, and its emissions (the
// log-likelihood of its observations given z = 1 and given z = 0) come from
// the shared functions in fp_emissions.stanfunctions. Units are numbered
// (i - 1) * nseason + t, so the units of a site are contiguous and in time
// order; a unit without data has emissions 0 in both states and enters only
// through the transitions. known_present marks units with a verified
// detection or a confirmed true observation, which cannot be unoccupied.
//
// The states are summed out by the forward algorithm of a two-state hidden
// Markov model, written by hand because the transition matrix varies by
// site and season (hmm_marginal() takes one per sequence). Site random
// effects on the false-positive rate and on true detection belong to sites,
// not units. The sum over sites runs with reduce_sum(), each slice of sites
// with all their seasons. With one season the log density equals that of
// fp_occu.stan, apart from the priors of beta_col and beta_ext.
//
// Generated quantities add, per season: the expected occupancy averaged over
// sites (psi_mean), the realized occupancy, the share of sites occupied with
// the states drawn from their joint posterior by forward filtering, backward
// sampling (psi_fs), the turnover, the probability that an occupied site was
// unoccupied in the season before (turnover, seasons 2 to nseason), and, with
// intercept-only colonization and extinction, the equilibrium occupancy
// gamma / (gamma + eps) (psi_eq).
// ============================================================================
functions {
#include fp_emissions.stanfunctions

  // Marginal log-likelihood of sites a to b over their seasons
  vector dyn_loglik(int a, int b, int nseason, vector beta, vector beta_col,
                    vector beta_ext, matrix Xcol, matrix Xext, vector alpha,
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
    int ns = b - a + 1;
    int u0 = (a - 1) * nseason + 1;
    int u1 = b * nseason;
    int nu = u1 - u0 + 1;
    tuple(vector[nu], vector[nu]) em = unit_emissions(u0, u1, alpha, gamma0, gamma_fp, e_fp,
                                                      u, d_det, delta0, delta_s, design,
                                                      automated, verified, confirmed, options);
    vector[ns] lpsi = design.1[a:b] * beta;
    vector[nu] lcol = Xcol[u0:u1] * beta_col;
    vector[nu] lext = Xext[u0:u1] * beta_ext;
    vector[ns] ll;
    for (i in 1:ns) {
      int k = (i - 1) * nseason + 1;
      // log P(z, observations so far) for z = 1 and z = 0
      real l1 = log_inv_logit(lpsi[i]) + em.1[k];
      real l0 = log1m_inv_logit(lpsi[i])
                + (verified.3[u0 + k - 1] == 1 ? negative_infinity() : em.2[k]);
      for (t in 2:nseason) {
        k += 1;
        real n1 = log_sum_exp(l1 + log1m_inv_logit(lext[k]), l0 + log_inv_logit(lcol[k]))
                  + em.1[k];
        real n0 = log_sum_exp(l1 + log_inv_logit(lext[k]), l0 + log1m_inv_logit(lcol[k]))
                  + (verified.3[u0 + k - 1] == 1 ? negative_infinity() : em.2[k]);
        l1 = n1;
        l0 = n0;
      }
      ll[i] = log_sum_exp(l1, l0);
    }
    return ll;
  }

  // partial sum over a slice of sites, for reduce_sum()
  real partial_loglik_dyn(array[] int slice, int a, int b, int nseason, vector beta,
                          vector beta_col, vector beta_ext, matrix Xcol, matrix Xext,
                          vector alpha, real gamma0, vector gamma_fp, vector e_fp,
                          vector u, vector d_det, real delta0, vector delta_s,
                     tuple(matrix, matrix, matrix, matrix, matrix) design,
                     tuple(array[] int, array[] int, array[] int, array[] int,
                           vector, vector, array[] int, array[] int,
                           array[] int, array[] int, vector, vector, vector) automated,
                     tuple(vector, vector, array[] int, matrix, vector,
                           array[] int, array[] int) verified,
                     tuple(array[] int, array[] int, array[] int, matrix,
                           array[] int, array[] int, matrix) confirmed,
                     array[] int options) {
    return sum(dyn_loglik(a, b, nseason, beta, beta_col, beta_ext, Xcol, Xext, alpha,
                          gamma0, gamma_fp, e_fp, u, d_det, delta0, delta_s, design,
                          automated, verified, confirmed, options));
  }
}

data {
  int<lower=1> nsite;
  int<lower=1> nseason;
  int<lower=1> nunit;                        // site-seasons: nsite * nseason
  int<lower=1> Kbeta;
  matrix[nsite, Kbeta] X;                    // initial occupancy design (with intercept)
  int<lower=1> Kcol;
  matrix[nunit, Kcol] Xcol;                  // colonization design (with intercept)
  int<lower=1> Kext;
  matrix[nunit, Kext] Xext;                  // extinction design (with intercept)

  int<lower=0> N2;                           // rows of automated occasions
  array[N2] int<lower=1> n2;                 // occasions in each row
  array[N2] int<lower=0> y2;                 // detections in each row
  array[N2] int<lower=1, upper=nunit> site2;  // sorted by site
  array[nunit] int<lower=1> start2;          // first row of each site in y2
  array[nunit] int<lower=0> len2;            // number of rows of each site
  int<lower=1> Kalpha;
  matrix[N2, Kalpha] Z;                      // detection design (column 1 = intercept)
  int<lower=1> nrec;
  array[N2] int<lower=1, upper=nrec> rec2;
  int<lower=0, upper=1> estimate_re;         // recorder random intercept on p11

  int<lower=0> Kfp;
  int<lower=0, upper=1> fp_by_row;           // false-positive covariates vary by occasion
  matrix[nunit, fp_by_row == 1 ? 0 : Kfp] Xfp;  // false-positive design (no intercept)
  matrix[N2, fp_by_row == 1 ? Kfp : 0] Xfp2;    // ... by automated row
  int<lower=0, upper=1> estimate_fp;

  int<lower=0> N1;                           // verified occasions (long format)
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=nunit> site1;
  int<lower=0> Kver_s;                       // slopes of the verified detection
  matrix[N1, Kver_s] X1s;                    // their design (verified rows sorted by site)
  array[nunit] int<lower=0, upper=1> known_present;

  int<lower=0, upper=1> obs_confirm;         // observation-confirmation design
  int<lower=0> N3;                           // confirmed occasions
  array[N3] int<lower=0, upper=3> v3;        // confirmation state
  array[N3] int<lower=1, upper=nunit> site3;
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
  if (nunit != nsite * nseason) reject("nunit must be nsite * nseason");
  vector[N2] y2v = to_vector(y2);
  vector[N2] n2v = to_vector(n2);
  // sites with an automated detection cannot be unoccupied without false positives
  array[nunit] int any_det2 = rep_array(0, nunit);
  for (n in 1:N2) {
    if (y2[n] > n2[n]) reject("y2 cannot exceed n2");
    if (n > 1 && site2[n] < site2[n - 1]) reject("automated rows must be sorted by site");
    if (y2[n] > 0) any_det2[site2[n]] = 1;
  }
  if (obs_confirm == 1 && estimate_fp == 0) {
    reject("the observation-confirmation design needs a false-positive process (estimate_fp = 1)");
  }
  // verified occasions and detections per site
  vector[nunit] n1v = rep_vector(0, nunit);
  vector[nunit] s1v = rep_vector(0, nunit);
  for (n in 1:N1) {
    n1v[site1[n]] += 1;
    s1v[site1[n]] += y1[n];
  }
  // with covariates on the verified detection, its rows are used one by one
  vector[N1] y1v = to_vector(y1);
  array[nunit] int len1 = rep_array(0, nunit);
  for (n in 1:N1) {
    if (Kver_s > 0 && n > 1 && site1[n] < site1[n - 1]) reject("verified occasions must be sorted by site");
    len1[site1[n]] += 1;
  }
  array[nunit] int start1;
  start1[1] = 1;
  for (i in 2:nunit) start1[i] = start1[i - 1] + len1[i - 1];
  // first row and number of rows of each site among the confirmed occasions
  array[nunit] int len3 = rep_array(0, nunit);
  for (n in 1:N3) {
    if (n > 1 && site3[n] < site3[n - 1]) reject("confirmed occasions must be sorted by site");
    len3[site3[n]] += 1;
  }
  array[nunit] int start3;
  start3[1] = 1;
  for (i in 2:nunit) start3[i] = start3[i - 1] + len3[i - 1];
  // sites, sliced among threads by reduce_sum()
  array[nsite] int site_index = linspaced_int_array(nsite, 1, nsite);
  // site of each unit
  array[nunit] int unit_site;
  for (n in 1:nunit) unit_site[n] = (n - 1) %/% nseason + 1;
  // site effects on the false-positive and true-detection rates when there
  // are none (see model)
  vector[nunit] zero_site = rep_vector(0, nunit);
  // totals of automated occasions and detections per site (unoccupied branch)
  vector[nunit] Ysite = rep_vector(0, nunit);
  vector[nunit] Nsite = rep_vector(0, nunit);
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
  vector[nunit] lchoose_site = rep_vector(0, nunit);
  for (n in 1:N2) lchoose_site[site2[n]] += lchoose(n2[n], y2[n]);
  // recorder of each site, when it does not change within sites
  int rec_by_site = 1;
  array[nunit] int rec_site = rep_array(1, nunit);
  for (n in 1:N2) {
    if (n == 1 || site2[n] != site2[n - 1]) rec_site[site2[n]] = rec2[n];
    else if (rec2[n] != rec_site[site2[n]]) rec_by_site = 0;
  }
  // the data of the likelihood, grouped for site_loglik() and reduce_sum();
  // the order of the fields is the order of the unpacking in site_loglik()
  // design.1 is not used by the emissions
  tuple(matrix[nsite, Kbeta], matrix[N2, Kalpha], matrix[N2, Kalpha - 1],
        matrix[nunit, fp_by_row == 1 ? 0 : Kfp], matrix[N2, fp_by_row == 1 ? Kfp : 0])
    design = (X, Z, Z_slopes, Xfp, Xfp2);
  tuple(array[N2] int, array[N2] int, array[N2] int, array[N2] int, vector[N2],
        vector[N2], array[nunit] int, array[nunit] int, array[nunit] int,
        array[nunit] int, vector[nunit], vector[nunit], vector[nunit])
    automated = (y2, n2, site2, rec2, y2v, n2v, start2, len2, any_det2, rec_site,
                 Ysite, Nsite, lchoose_site);
  tuple(vector[nunit], vector[nunit], array[nunit] int, matrix[N1, Kver_s],
        vector[N1], array[nunit] int, array[nunit] int)
    verified = (n1v, s1v, known_present, X1s, y1v, start1, len1);
  tuple(array[N3] int, array[N3] int, array[N3] int, matrix[N3, Kalpha],
        array[nunit] int, array[nunit] int, matrix[N3, fp_by_row == 1 ? Kfp : 0])
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
  vector[Kcol] beta_col;
  vector[Kext] beta_ext;
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
  beta_col[1]    ~ normal(0, prior_sd_int);
  if (Kcol > 1)   beta_col[2:Kcol] ~ normal(0, prior_sd_slope);
  beta_ext[1]    ~ normal(0, prior_sd_int);
  if (Kext > 1)   beta_ext[2:Kext] ~ normal(0, prior_sd_slope);

  // as in fp_occu.stan: site effects (expanded to units) when estimated,
  // otherwise the data vector zero_site
  if (estimate_fp_re == 1 && estimate_det_re == 1) {
    vector[nunit] e_fp = (sigma_fp[1] * e_raw)[unit_site];
    vector[nunit] d_det = (sigma_det[1] * d_raw)[unit_site];
    target += reduce_sum(partial_loglik_dyn, site_index, 1, nseason, beta, beta_col,
                         beta_ext, Xcol, Xext, alpha, gamma0, gamma_fp, e_fp, u, d_det,
                         delta0, delta_s, design, automated, verified, confirmed, options);
  } else if (estimate_fp_re == 1) {
    vector[nunit] e_fp = (sigma_fp[1] * e_raw)[unit_site];
    target += reduce_sum(partial_loglik_dyn, site_index, 1, nseason, beta, beta_col,
                         beta_ext, Xcol, Xext, alpha, gamma0, gamma_fp, e_fp, u, zero_site,
                         delta0, delta_s, design, automated, verified, confirmed, options);
  } else if (estimate_det_re == 1) {
    vector[nunit] d_det = (sigma_det[1] * d_raw)[unit_site];
    target += reduce_sum(partial_loglik_dyn, site_index, 1, nseason, beta, beta_col,
                         beta_ext, Xcol, Xext, alpha, gamma0, gamma_fp, zero_site, u, d_det,
                         delta0, delta_s, design, automated, verified, confirmed, options);
  } else {
    target += reduce_sum(partial_loglik_dyn, site_index, 1, nseason, beta, beta_col,
                         beta_ext, Xcol, Xext, alpha, gamma0, gamma_fp, zero_site, u,
                         zero_site, delta0, delta_s, design, automated, verified,
                         confirmed, options);
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
  vector[nseason] psi_mean;
  vector[nseason] psi_fs;
  vector[nseason - 1] turnover;
  array[Kcol == 1 && Kext == 1 && nseason > 1 ? 1 : 0] real psi_eq;
  if (size(psi_eq) == 1) {
    psi_eq[1] = inv_logit(beta_col[1]) / (inv_logit(beta_col[1]) + inv_logit(beta_ext[1]));
  }
  {
    vector[nunit] e_fp = rep_vector(0, nunit);
    if (estimate_fp_re == 1) e_fp = (sigma_fp[1] * e_raw)[unit_site];
    vector[nunit] d_det = rep_vector(0, nunit);
    if (estimate_det_re == 1) d_det = (sigma_det[1] * d_raw)[unit_site];
    // as in fp_occu.stan, with site-seasons in place of sites
    vector[nunit] lp10 = rep_vector(gamma0, nunit);
    if (Kfp > 0 && fp_by_row == 0) lp10 += Xfp * gamma_fp;
    lp10 += e_fp;
    if (estimate_fp == 1 && fp_by_row == 0) p10_bar = mean(inv_logit(lp10));
    if (estimate_fp == 1 && fp_by_row == 1 && N2 > 0) {
      p10_bar = dot_product(n2v, inv_logit(lp10[site2] + Xfp2 * gamma_fp)) / sum(n2v);
    }
    // forward over the seasons of each site, as in dyn_loglik(), keeping
    // log P(z_t, observations up to t) of every unit for the backward pass
    tuple(vector[nunit], vector[nunit]) em = unit_emissions(1, nunit, alpha, gamma0, gamma_fp,
                                                            e_fp, u, d_det, delta0, delta_s,
                                                            design, automated, verified,
                                                            confirmed, options);
    vector[nsite] lpsi = X * beta;
    vector[nunit] lcol = Xcol * beta_col;
    vector[nunit] lext = Xext * beta_ext;
    vector[nunit] f1;
    vector[nunit] f0;
    vector[nunit] pz;                        // occupancy from the covariates alone
    array[nunit] int z;
    for (i in 1:nsite) {
      int k = (i - 1) * nseason + 1;
      f1[k] = log_inv_logit(lpsi[i]) + em.1[k];
      f0[k] = log1m_inv_logit(lpsi[i]) + (known_present[k] == 1 ? negative_infinity() : em.2[k]);
      pz[k] = inv_logit(lpsi[i]);
      for (t in 2:nseason) {
        k += 1;
        f1[k] = log_sum_exp(f1[k - 1] + log1m_inv_logit(lext[k]),
                            f0[k - 1] + log_inv_logit(lcol[k])) + em.1[k];
        f0[k] = log_sum_exp(f1[k - 1] + log_inv_logit(lext[k]),
                            f0[k - 1] + log1m_inv_logit(lcol[k]))
                + (known_present[k] == 1 ? negative_infinity() : em.2[k]);
        pz[k] = pz[k - 1] * (1 - inv_logit(lext[k])) + (1 - pz[k - 1]) * inv_logit(lcol[k]);
      }
      log_lik[i] = log_sum_exp(f1[k], f0[k]);
      // backward sampling: the last season, then each season given the next
      z[k] = bernoulli_rng(inv_logit(f1[k] - f0[k]));
      for (t in 1:(nseason - 1)) {
        k -= 1;
        real w = z[k + 1] == 1
                 ? (f1[k] + log1m_inv_logit(lext[k + 1])) - (f0[k] + log_inv_logit(lcol[k + 1]))
                 : (f1[k] + log_inv_logit(lext[k + 1])) - (f0[k] + log1m_inv_logit(lcol[k + 1]));
        z[k] = bernoulli_rng(inv_logit(w));
      }
    }
    // per season, over sites
    for (t in 1:nseason) {
      array[nsite] int k = linspaced_int_array(nsite, t, t + (nsite - 1) * nseason);
      psi_mean[t] = mean(pz[k]);
      psi_fs[t] = mean(to_vector(z[k]));
      if (t > 1) {
        array[nsite] int k0 = linspaced_int_array(nsite, t - 1, t - 1 + (nsite - 1) * nseason);
        turnover[t - 1] = dot_product(1 - pz[k0], inv_logit(lcol[k])) / sum(pz[k]);
      }
    }
  }
}
