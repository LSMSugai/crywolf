#' Output names and interpretation
#'
#' A reference for the parameter names and result fields returned by crywolf.
#' Use `help("crywolf-output")` to open this page. Coefficients keep the same
#' names in `as_draws_df(fit)`, `summary(fit)$table`, and `coef(fit)`.
#'
#' @section Stability:
#' From crywolf 1.9.0, the names documented on this page are frozen: the
#' readable parameter and derived-quantity names, the columns of the
#' summaries and of the other results described below, and the fields of
#' the fitted object listed under "Readable draws and raw Stan names". They
#' change only in a major version, announced in the previous minor version.
#' New options may add names (as `random_det` added `sd_det_site`) without
#' changing existing ones. Raw Stan names and the fields of `design` and
#' `stan_data` are implementation details and are not covered. From 2.0.0
#' the guarantee extends to arguments, classes and the results of the
#' other stable functions: see [crywolf-api].
#'
#' @section Coefficient names:
#' Suffixes come from the columns of R's [stats::model.matrix()]: `(Intercept)`
#' for an intercept, the covariate name for a numeric slope, factor contrasts
#' such as `habitatforest`, and interaction columns such as `temp:rh`.
#' They refer to the covariate scale supplied to the fit; crywolf does not
#' automatically standardize formula covariates.
#'
#' | Name | Meaning and availability |
#' | :--- | :--- |
#' | `occ_(Intercept)`, `occ_<term>` | Coefficients of `occ`, on the logit scale of occupancy probability. Always present. |
#' | `det_(Intercept)`, `det_<term>` | Coefficients of `det`, on the logit scale of automated detection at occupied sites in the Miller design; in the Chambert design, of the true-detection component before combining true and false detections. Always present. |
#' | `det_y_prev` | Coefficient of the previous automated detection when `markov` is supplied. A log odds change within the same site, not a transition in occupancy. |
#' | `fp_(Intercept)`, `fp_<term>` | Coefficients of `fp`, on the logit scale of the false-positive probability at unoccupied sites. Present when `fp` is not `NULL`. |
#' | `verified_(Intercept)`, `verified_<term>` | Coefficients of `ver`, on the logit scale of verified-method detection at occupied sites. Present only when verified occasions enter the fit. |
#' | `fp_gap` | Positive difference `det_(Intercept) - fp_(Intercept)` on the logit scale. Present when `fp` is not `NULL`. It orders the reference intercepts, not all covariate-specific probabilities. |
#' | `sd_recorder` | Standard deviation of the recorder random intercept on the detection logit scale, when the recorder effect is enabled. |
#' | `sd_fp_site` | Standard deviation of the site random intercept on the false-positive logit scale, when `random_fp = TRUE`. |
#' | `sd_det_site` | Standard deviation of the site random intercept on the true-detection logit scale, when `random_det = TRUE`. |
#'
#' `coef(fit, "occ")`, `coef(fit, "det")`, and `coef(fit, "fp")`
#' select the corresponding prefixes. `coef(fit, "fp")` excludes `fp_gap`.
#' `coef(fit)` returns posterior means of all readable draw variables,
#' including derived probabilities and standard deviations.
#'
#' @section Derived probabilities in posterior draws:
#' Reference probabilities set design columns other than the intercept to
#' zero and random effects to zero. This means the reference factor level
#' under treatment contrasts, and `y_prev = 0` under Markov detection.
#' Reference values need not describe an observed site.
#'
#' | Name | Meaning and availability |
#' | :--- | :--- |
#' | `p11_ref` | Probability of any automated detection at an occupied reference site. In Miller, `plogis(det_(Intercept))`. In Chambert, `1 - (1 - plogis(det_(Intercept))) * (1 - plogis(fp_(Intercept)))`. Always present. |
#' | `p10_ref` | `plogis(fp_(Intercept))`: false-positive probability at the reference covariates, with zero site effect. Present when `fp` is not `NULL`. |
#' | `p10_bar` | Average false-positive probability over the fitted data, including fitted site effects. With site-level FP covariates, an unweighted mean over sites; with occasion-level FP covariates, a mean over unchecked automated occasions weighted by their trial counts. Confirmed occasions are excluded from this occasion average. Present when `fp` is not `NULL`. |
#' | `p_verified` | `plogis(verified_(Intercept))`: verified-method detection probability at reference covariates, not an average over verified occasions. Present only with verified occasions. |
#'
#' When `fp = NULL`, false-positive probabilities are zero, but `fp_*`,
#' `p10_ref`, and `p10_bar` are omitted from the readable draws.
#' The Chambert observation-confirmation design requires a false-positive
#' process; `fp = NULL` is rejected for that design.
#'
#' @section Posterior summaries and fit diagnostics:
#' `summary(fit)` returns `table`, `prob`, `diagnostics`, and `call`.
#' Its `table` has one row per readable variable:
#'
#' | Column | Meaning |
#' | :--- | :--- |
#' | `variable` | Readable parameter or derived-quantity name. |
#' | `mean`, `sd` | Posterior mean and standard deviation. |
#' | `lower`, `upper` | Equal-tailed credible interval, with width `prob`. |
#' | `rhat` | R-hat for this variable. |
#' | `ess_bulk`, `ess_tail` | Bulk and tail effective sample sizes. |
#'
#' `fit$diagnostics` contains `divergences` (total divergent transitions),
#' `max_treedepth` (the number of transitions reaching the maximum depth,
#' not the maximum depth observed), `max_rhat` (largest R-hat),
#' `min_ess_bulk`, and `min_ess_tail` (smallest effective sample sizes).
#' R-hat and the effective sample sizes are summarized across the readable
#' variables; `divergences` and `max_treedepth` count sampler transitions.
#'
#' @section Predictions and site estimates:
#' `predict(fit)` appends `estimate`, `lower`, and `upper` to `newdata`:
#' the posterior mean probability and equal-tailed interval for `process`.
#' With `summary = FALSE`, it returns a matrix with posterior draws in rows
#' and prediction rows in columns. These predictions use fixed effects;
#' site or recorder random effects are not added or integrated out.
#'
#' `fp_effects()` returns `covariate`, `process`, `value` (the covariate on
#' the fitting scale), `x` (its display value after any `transform`),
#' `estimate`, `lower`, and `upper`, with interval width in attribute `prob`.
#'
#' `fp_sites()` retains the input site columns and adds:
#'
#' | Column | Meaning |
#' | :--- | :--- |
#' | `psi`, `psi_lower`, `psi_upper` | Posterior mean and interval for occupancy probability given site covariates. |
#' | `occupied`, `occupied_lower`, `occupied_upper` | Posterior mean and interval of the probability of occupancy conditional on the site's observations and each parameter draw. The interval describes parameter uncertainty in that probability, not an interval for a binary state. |
#' | `detections` | Number of detections in unchecked automated occasions. Confirmed occasions are excluded. |
#' | `verified` | Whether occupancy is certain from a verified detection or a confirmed true observation (Chambert). Despite the name, it includes both sources of evidence. |
#' | `p10` | Posterior mean FP probability including the site's fitted random effect, if any. With occasion-level FP covariates, the trial-weighted average over the site's unchecked automated occasions. Present only with false positives. |
#'
#' For occasion-level FP covariates, the current implementation returns
#' `p10 = 0` for a site without unchecked automated occasions; this is an
#' empty-data convention, not evidence for a zero false-positive rate.
#'
#' @section Checks and threshold sweeps:
#' `fp_ppc()` returns `y` (observed detection counts by retained site),
#' `yrep` (replicated counts, draws by sites), `fit_obs`, `fit_rep`
#' (discrepancies per draw), `p_value` (the proportion with
#' `fit_rep > fit_obs`), `site` (retained identifiers), `method`, and `fit_stat`.
#'
#' `fp_prior_predictive()` has columns `naive_occupancy` (share of all sites
#' with an automated detection), `automated_rate`, `verified_rate`,
#' `mean_psi` (mean occupancy probability), and `p10_bar`. Attribute `observed`
#' holds the first three summaries for the observed data. Each row is one
#' simulated dataset from the priors.
#'
#' `fp_threshold_sweep()` has columns `threshold`, `detection_rate` (share
#' of occasions detected), `parameter` (readable name), `mean`, `lower`,
#' `upper`, `divergences`, and `max_rhat`. Attribute `prob` holds the interval
#' width. Named threshold selection adds `threshold_selection` (one row per
#' method) and `validation_metrics` (one row per candidate score):
#' `threshold`, `tp`, `fp`, `fn`, `precision`, `recall`, and `f1`;
#' `threshold_selection` also has `method`. Here `fp` is a false-positive
#' count, not a model coefficient. See [fp_threshold_sweep()] for definitions.
#'
#' @section Calibration and simulation:
#' `fp_calibrate()` and `fp_calibrate_fit()` return `cell`, `rep`, `seed`,
#' `model`, `parameter`, `truth`, `mean`, `lower`, `upper`, `covered`,
#' `divergences`, `max_rhat`, `ok`, and `seconds`; `fp_calibrate()` adds the
#' columns of its `design`.
#' `ok` records whether fitting returned a result, not convergence.
#' `summary()` groups by cell, model and parameter and returns the design
#' columns plus `n` (successful fits), `failed`, `coverage`, `mc_se`
#' (Monte Carlo standard error of coverage), `bias`, `rmse`, and
#' `prop_divergent` (fraction of successful fits with divergences).
#' `fp_verification_curve()` returns the calibration object invisibly;
#' its printed curve selects `model`, `frac_verified`, `coverage`,
#' `mc_se`, and `bias` for `occ_x`.
#'
#' `fp_sbc()` returns `rep`, `parameter`, `truth`, `rank` (number of
#' retained draws below the true value), `divergences`, and `ok`, with
#' attribute `ndraws`. Its `summary()` returns `parameter`, `n`, `mean_rank`
#' (mean rank divided by `ndraws`), `p_uniform` (chi-squared rank-histogram
#' test p-value), and `divergent` (fraction with divergences); attribute
#' `failed` counts failed replicates. In SBC only, parameter `log_lik` is
#' the total log likelihood, rather than a site's contribution.
#'
#' `simulate_fp()` returns `data` and `truth`. The latter records the
#' generating arguments (`psi`, `b1`, `p11`, `p10`, `r11`, `sigma_p10`,
#' `p10_conc`, `conf_effect`, `conf_cor`, `a_bio`, `a_aci`, `sd_recorder`,
#' `frac_confirmed`), `b0` (occupancy intercept), `regime`,
#' `conf_cor_realized`, `n_verified_sites`, `z` (binary state by site),
#' `p10_site`, and `p10_marginal` (mean of `p10_site`).
#' `simulate(fit)` returns a list of datasets, each with `data` and `truth`;
#' its truth fields are `pars` (named generating coefficients), `z`,
#' `p10_site`, and `p10_bar`.
#'
#' @section Multispecies model (in development):
#' [fp_occu_ms()] gives the community means under the names above
#' (`occ_*`, `det_*`, `fp_*`, `fp_gap`, which is the median gap among
#' species, and `p11_ref`, `p10_ref` at the community means), the
#' standard deviations among species as `sd_occ_<term>`,
#' `sd_det_(Intercept)` and `sd_log_fp_gap` (of the logarithm of the gap),
#' and the values of each species with its name in brackets:
#' `occ_<term>[species]`, `det_(Intercept)[species]`, `fp_gap[species]`,
#' `p11_ref[species]` and `p10_ref[species]`. It has no `p10_bar`.
#' `log_lik[i]` is the log-likelihood of unit `i`, site
#' `(i - 1) %/% nspecies + 1` and species `(i - 1) %% nspecies + 1`.
#'
#' @section Dynamic model (experimental):
#' [fp_occu_dyn()] uses the names above, except that the coefficients of
#' initial occupancy are `psi1_(Intercept)` and `psi1_<term>` (in place of
#' `occ_*`), and adds the coefficients of colonization, `col_(Intercept)`
#' and `col_<term>`, and of extinction, `ext_(Intercept)` and
#' `ext_<term>`, on the logit scale of the probability per transition
#' between consecutive seasons. `p10_bar` averages over site-seasons in
#' place of sites. `log_lik[i]` is the log-likelihood of site `i` over all
#' its seasons. The raw Stan names are `beta[k]`, `beta_col[k]` and
#' `beta_ext[k]`.
#'
#' Derived quantities, one per season `t` (numbered in the order of the
#' seasons; `summary()` gives the season of each in its `season` column):
#' `psi_mean[t]`, the occupancy probability expected from the covariates
#' (initial occupancy carried forward by colonization and extinction),
#' averaged over sites; `psi_fs[t]`, the realized occupancy, the share of
#' sites occupied, with the states drawn from their joint posterior given
#' the data by forward filtering, backward sampling; `turnover[t]`, for
#' `t >= 2`, the probability that a site occupied in season `t` was
#' unoccupied in season `t - 1`, \eqn{\sum_i (1 - \psi_{i,t-1})
#' \gamma_{it} / \sum_i \psi_{it}}{sum_i (1 - psi_i,t-1) gamma_it / sum_i
#' psi_it}, with the expected occupancy \eqn{\psi}{psi}; and, when
#' colonization and extinction have only an intercept, `psi_eq`, the
#' equilibrium occupancy \eqn{\gamma / (\gamma + \epsilon)}{gamma / (gamma
#' + eps)}. Being experimental, these names may change in a minor version
#' (see [crywolf-api]).
#'
#' @section Readable draws and raw Stan names:
#' The fitted object contains `draws`, `fit` (the backend fit), `data`,
#' `formulas`, `markov`, `random_recorder`, `random_fp`, `random_det`, `priors`,
#' `sampler`, `diagnostics`, `call`, `design`, and `stan_data`.
#' `sampler` records `chains`, `parallel_chains`, `threads_per_chain`,
#' `iter_warmup`, `iter_sampling`, and `adapt_delta`.
#'
#' `as_draws_df(fit)` returns the readable variables above, with posterior
#' metadata `.chain`, `.iteration`, and `.draw`. The backend's
#' `fit$fit$draws()` also exposes raw Stan names. Current mappings are:
#'
#' | Stan name | Readable name |
#' | :--- | :--- |
#' | `beta[k]` | `occ_<term>` |
#' | `alpha[k]` | `det_<term>` |
#' | `gamma0`, `gamma_fp[k]` | `fp_(Intercept)`, `fp_<term>` |
#' | `delta0`, `delta_s[k]` | `verified_(Intercept)`, `verified_<term>` |
#' | `r11` | `p_verified` |
#' | `sigma_u`, `sigma_fp[1]`, `sigma_det[1]` | `sd_recorder`, `sd_fp_site`, `sd_det_site` |
#'
#' `fp_gap`, `p11_ref`, `p10_ref`, and `p10_bar` keep their names.
#' `log_lik[i]` is the state-marginalized log likelihood of site `i`, in
#' input site-table order; [loo()] uses sites as its pointwise units.
#' It is available from `fit$fit$draws("log_lik")`, not the readable draws.
#' With a fitted site random effect, state marginalization does not also
#' integrate that effect out; [loo()] does, by quadrature (its argument
#' `integrate`). Raw random effects (`u`, `u_raw`, `e_raw`, `d_raw`)
#' and sampler variable `lp__` are backend outputs, not additional public
#' coefficient names. `design` and `stan_data` describe internal encodings;
#' their field names and raw Stan names are implementation details.
#'
#' @seealso [fp_occu()], [fp_occu-methods], [fp_sites()], [fp_effects()],
#'   [fp_ppc()], [fp_threshold_sweep()], [fp_calibrate()], [fp_sbc()]
#' @name crywolf-output
#' @aliases output-names
#' @keywords documentation
NULL
