# Simulations of the multispecies model (`fp_occu_ms()`)

`recovery.R` is the short recovery check of PR 2: one simulated community
(12 species, 80 sites, 8 automated occasions, 3 verified occasions at 30%
of the sites; occupancy `~ x` and detection `~ bio` by species, constant
false-positive rate with a gap by species), fitted once with 4 chains of
1000 + 1000 iterations. The simulator in the script is provisional;
`simulate_fp_ms()` comes in PR 4, and the simulations and SBC on AWS in
PR 5. Results in `recovery.rds`.

## Run of 02/11/2026

On the Mac, 2 threads per chain, 7.1 minutes. 43 of the 46 95% intervals
(community means, standard deviations among species, occupancy and
detection intercepts and `p10_ref` of each species) hold the true value;
the three misses are two occupancy intercepts and one `p10_ref`. Largest
R-hat 1.007, no divergences.
