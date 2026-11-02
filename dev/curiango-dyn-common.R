# Curiango data for the dynamic model (the nights of each site as seasons)
# and for the model of independent site-nights, as described in the header
# of dev/dyn-kfold/run_kfold.R. Sourced by that script and by
# dev/precompute-dynamic.R.
zs <- function(v) as.numeric(scale(v))
m <- curiango$minutes
nt <- curiango$nights
daily <- stats::aggregate(cbind(temp, rh) ~ night, nt, mean)
daily <- data.frame(night = daily$night, temp_d = zs(daily$temp), rh_d = zs(daily$rh))
i <- match(m$site_night, nt$site_night)
occasions <- data.frame(file = paste(m$site_night, format(m$time, "%Y%m%d_%H%M")),
                        site = m$site, site_night = m$site_night, night = nt$night[i],
                        time = m$time, bio = zs(m$bio), aci = zs(m$aci),
                        temp = zs(nt$temp[i]), rh = zs(nt$rh[i]), recorder = m$site,
                        y = as.integer(m$conf_custom >= 0.40))
verified <- data.frame(site = curiango$verified$site,
                       site_night = curiango$verified$site_night,
                       night = nt$night[match(curiango$verified$site_night, nt$site_night)],
                       y = curiango$verified$y)
stopifnot(!anyNA(occasions$night), !anyNA(verified$night))
site_nights <- merge(nt[c("site_night", "site", "night")], daily, by = "night")

forms <- list(det = ~ bio + aci, fp = ~ temp + rh, occ = ~ temp_d + rh_d,
              psi1 = ~ 1, col = ~ temp_d + rh_d, ext = ~ temp_d + rh_d)
