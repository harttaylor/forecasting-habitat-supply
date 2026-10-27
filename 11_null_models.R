# =============================================================================
# 11_null_models.R
# First bird models: QPAD detection offsets + survey design only, no habitat variables. 
# also get a baseline density for each species (males/ha), and choose the model structure 
# (Poisson vs negative binomial, random effects) before adding habitat covariates
#
# Offsets: BAM QPAD version 3, "constant" models (one singing rate and one
# detection distance per species). Upgrade to time-of-day / day-of-year /
# forest-cover models later
#
# QPAD install:  remotes::install_github("psolymos/QPAD")
# Input: data/derived/bird_data_fit_candidates.csv created in 10_fit_region_birds_and_map.R
# =============================================================================

library(dplyr)
library(glmmTMB)
library(QPAD)

proj_dir <- "C:/Users/hartt/Documents/Chapter 5/Analysis"
out_dir  <- file.path(proj_dir, "data/derived")

one_visit_per_year <- TRUE    # TRUE = one random survey per location per year (ARU sites can have 50+)

# ---- 1. prep data ---------------------------------------------------------------------------
birds <- read.csv(file.path(out_dir, "bird_data_fit_candidates.csv"))
fit   <- filter(birds, keep)

# Maximum survey duration (min)
dur_txt  <- sub("min.*", "", fit$survey_duration_method)
fit$dur  <- suppressWarnings(as.numeric(sub(".*-", "", dur_txt)))

# Maximum survey distance (m): unlimited for ARUs and "INF" methods, e.g. "0m-50m" -> 50
dist_txt     <- sub(".*-", "", gsub("m", "", fit$survey_distance_method))
fit$maxdist  <- suppressWarnings(as.numeric(dist_txt))
fit$maxdist[grepl("INF", fit$survey_distance_method)] <- Inf

fit$method <- ifelse(grepl("ARU", fit$survey_distance_method), "ARU", "human")
fit <- filter(fit, !is.na(dur), dur > 0, !is.na(maxdist))

if (one_visit_per_year) {
  set.seed(1)
  fit <- fit %>% group_by(location, year) %>% slice_sample(n = 1) %>% ungroup() %>% as.data.frame()
}
cat("Surveys used:", nrow(fit), "  locations:", n_distinct(fit$location), "\n")
print(table(fit$method))

# ---- 2. QPAD offsets ---------------------------------------------------------------------------
# offset = log(p * A * q): p = chance a bird sings during the survey (from duration),
# A * q = effective area sampled in ha (from detection distance tau and max distance).
load_BAM_QPAD(version = 3)
r <- fit$maxdist / 100                               # QPAD works in 100 m units -> areas in ha

cf    <- coefBAMspecies("BTNW", 0, 0)
phi   <- exp(cf$sra);  tau <- exp(cf$edr)
p     <- 1 - exp(-fit$dur * phi)
Aq    <- ifelse(is.infinite(r), pi * tau^2, pi * tau^2 * (1 - exp(-r^2 / tau^2)))
fit$off_BTNW <- log(p * Aq)

cf    <- coefBAMspecies("BRCR", 0, 0)
phi   <- exp(cf$sra);  tau <- exp(cf$edr)
p     <- 1 - exp(-fit$dur * phi)
Aq    <- ifelse(is.infinite(r), pi * tau^2, pi * tau^2 * (1 - exp(-r^2 / tau^2)))
fit$off_BRCR <- log(p * Aq)

cf    <- coefBAMspecies("CAWA", 0, 0)
phi   <- exp(cf$sra);  tau <- exp(cf$edr)
p     <- 1 - exp(-fit$dur * phi)
Aq    <- ifelse(is.infinite(r), pi * tau^2, pi * tau^2 * (1 - exp(-r^2 / tau^2)))
fit$off_CAWA <- log(p * Aq)


cat("\nOffset ranges (log ha):\n")
print(round(sapply(fit[, c("off_BTNW", "off_BRCR", "off_CAWA")], range), 2))

fit$year_s <- (fit$year - 2015) / 10     # year centred on 2015, per decade

# ---- 3. Null models: BTNW ------------------------------------------------------------------------
b0 <- glmmTMB(BTNW ~ 1 + offset(off_BTNW), family = poisson, data = fit)
b1 <- glmmTMB(BTNW ~ method + offset(off_BTNW) + (1 | location), family = poisson, data = fit)
b2 <- glmmTMB(BTNW ~ method + NSRNAME + offset(off_BTNW) + (1 | location), family = poisson, data = fit)
b3 <- glmmTMB(BTNW ~ method + NSRNAME + offset(off_BTNW) + (1 | location), family = nbinom2, data = fit)
b4 <- glmmTMB(BTNW ~ method + NSRNAME + year_s + offset(off_BTNW) + (1 | location), family = poisson, data = fit)

# ---- 4. Null models: BRCR ------------------------------------------------------------------------
r0 <- glmmTMB(BRCR ~ 1 + offset(off_BRCR), family = poisson, data = fit)
r1 <- glmmTMB(BRCR ~ method + offset(off_BRCR) + (1 | location), family = poisson, data = fit)
r2 <- glmmTMB(BRCR ~ method + NSRNAME + offset(off_BRCR) + (1 | location), family = poisson, data = fit)
r3 <- glmmTMB(BRCR ~ method + NSRNAME + offset(off_BRCR) + (1 | location), family = nbinom2, data = fit)
r4 <- glmmTMB(BRCR ~ method + NSRNAME + year_s + offset(off_BRCR) + (1 | location), family = poisson, data = fit)

# ---- 5. Null models: CAWA ------------------------------------------------------------------------
c0 <- glmmTMB(CAWA ~ 1 + offset(off_CAWA), family = poisson, data = fit)
c1 <- glmmTMB(CAWA ~ method + offset(off_CAWA) + (1 | location), family = poisson, data = fit)
c2 <- glmmTMB(CAWA ~ method + NSRNAME + offset(off_CAWA) + (1 | location), family = poisson, data = fit)
c3 <- glmmTMB(CAWA ~ method + NSRNAME + offset(off_CAWA) + (1 | location), family = nbinom2, data = fit)
c4 <- glmmTMB(CAWA ~ method + NSRNAME + year_s + offset(off_CAWA) + (1 | location), family = poisson, data = fit)

# ---- 6. Check it out! --------------------------------------------------------------------------------------
model_names <- c("0 intercept only", "1 + method + location RE", "2 + subregion",
                 "3 as 2, negative binomial", "4 as 2 + year trend")
aic_b <- c(AIC(b0), AIC(b1), AIC(b2), AIC(b3), AIC(b4))
aic_r <- c(AIC(r0), AIC(r1), AIC(r2), AIC(r3), AIC(r4))
aic_c <- c(AIC(c0), AIC(c1), AIC(c2), AIC(c3), AIC(c4))
aic <- data.frame(model = model_names,
                  BTNW = round(aic_b - min(aic_b, na.rm = TRUE), 1),
                  BRCR = round(aic_r - min(aic_r, na.rm = TRUE), 1),
                  CAWA = round(aic_c - min(aic_c, na.rm = TRUE), 1),
                  # TRUE = model converged properly (NA AIC probably means it did not)
                  ok_BTNW = c(b0$sdr$pdHess, b1$sdr$pdHess, b2$sdr$pdHess, b3$sdr$pdHess, b4$sdr$pdHess),
                  ok_BRCR = c(r0$sdr$pdHess, r1$sdr$pdHess, r2$sdr$pdHess, r3$sdr$pdHess, r4$sdr$pdHess),
                  ok_CAWA = c(c0$sdr$pdHess, c1$sdr$pdHess, c2$sdr$pdHess, c3$sdr$pdHess, c4$sdr$pdHess))
cat("\nDelta AIC (0 = best):\n"); print(aic)

# Average density across the survey points (males/ha) from the intercept-only model
dens <- data.frame(species = c("BTNW", "BRCR", "CAWA"),
                   males_per_ha = round(exp(c(fixef(b0)$cond[1], fixef(r0)$cond[1], fixef(c0)$cond[1])), 3),
                   naive_det_rate = round(c(mean(fit$BTNW > 0), mean(fit$BRCR > 0), mean(fit$CAWA > 0)), 3))
cat("\nBaseline density (no habitat):\n"); print(dens)

cat("\nYear trend (per decade, log scale) from model 4:\n")
print(round(rbind(BTNW = summary(b4)$coefficients$cond["year_s", ],
                  BRCR = summary(r4)$coefficients$cond["year_s", ],
                  CAWA = summary(c4)$coefficients$cond["year_s", ]), 3))

cat("\nSubregion effects (model 2; log scale, relative to the first subregion):\n")
print(round(cbind(BTNW = fixef(b2)$cond, BRCR = fixef(r2)$cond, CAWA = fixef(c2)$cond), 2))

write.csv(fit, file.path(out_dir, "fit_data_with_offsets.csv"), row.names = FALSE)
saveRDS(list(b0 = b0, b1 = b1, b2 = b2, b3 = b3, b4 = b4, r0 = r0, r1 = r1, r2 = r2, r3 = r3, r4 = r4,
             c0 = c0, c1 = c1, c2 = c2, c3 = c3, c4 = c4), file.path(out_dir, "null_models.rds"))
