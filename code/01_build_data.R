# =============================================================================
# 01_build_data.R — builds data/analysis_sample.csv from the raw inputs
#   data/raw/stage1_sample.csv   : 739 country-years (37 economies, 2003-2023), final DTI version
#   data/raw/panel_full.csv      : full 38-economy panel 2000-2023 (DTI sub-components, enrolment)
#   data/raw/wdi_IT.MLT.MAIN.P2_fixed_telephone.csv  : WDI fixed-telephone subscriptions /100
#   data/raw/wdi_IT.CEL.SETS.P2_mobile_subscriptions.csv : WDI mobile subscriptions /100
# Instruments:
#   z_tel = (mean fixed-telephone density 1990-1995) x (leave-one-out mean DTI of other countries in year t)
#   z_mob = (mean mobile density 1995-1999)          x (same leave-one-out global DTI path)
#   ict_imports_l1 (already in stage1_sample)
# =============================================================================
source("code/functions.R")
s1    <- read.csv("data/raw/stage1_sample.csv", stringsAsFactors = FALSE)
panel <- read.csv("data/raw/panel_full.csv", stringsAsFactors = FALSE)
ft    <- read.csv("data/raw/wdi_IT.MLT.MAIN.P2_fixed_telephone.csv")
ms    <- read.csv("data/raw/wdi_IT.CEL.SETS.P2_mobile_subscriptions.csv")

pre <- ft %>% filter(year %in% 1990:1995) %>% group_by(iso3c) %>% summarise(ft_pre = mean(value, na.rm = TRUE)) %>%
  left_join(ms %>% filter(year %in% 1995:1999) %>% group_by(iso3c) %>% summarise(mob_pre = mean(value, na.rm = TRUE)), by = "iso3c")

# The full panel stores an earlier scaling of the DTI; map it onto the final scale (r = 0.99995)
jj  <- s1 %>% select(iso3c, year, dti_s1 = dti_score) %>% left_join(panel %>% select(iso3c, year, dti_df = dti_score), by = c("iso3c", "year"))
map <- coef(lm(dti_s1 ~ dti_df, jj))
glob <- panel %>% mutate(dti_fin = map[1] + map[2] * dti_score) %>%
  group_by(year) %>% mutate(n = sum(!is.na(dti_fin)), s = sum(dti_fin, na.rm = TRUE)) %>% ungroup() %>%
  mutate(glob_loo = (s - ifelse(is.na(dti_fin), 0, dti_fin)) / (n - !is.na(dti_fin))) %>% select(iso3c, year, glob_loo)

# Alternative DTI without digital financial inclusion (robustness): PCA of fixed_bb, mobile_bb, egdi
sub3 <- c("fixed_bb", "mobile_bb", "egdi"); cc <- complete.cases(panel[, sub3])
p3 <- prcomp(panel[cc, sub3], center = TRUE, scale. = TRUE)
panel$dti3 <- NA; panel$dti3[cc] <- p3$x[, 1] * sign(p3$rotation[1, 1])

an <- s1 %>% left_join(pre, by = "iso3c") %>% left_join(glob, by = c("iso3c", "year")) %>%
  left_join(panel %>% select(iso3c, year, dti3), by = c("iso3c", "year")) %>%
  mutate(z_tel = ft_pre * glob_loo, z_mob = mob_pre * glob_loo,
         dti3 = (dti3 - mean(dti3)) / sd(dti3) * sd(dti_score) + mean(dti_score))   # rescaled to the DTI's mean/SD
stopifnot(nrow(an) == 739, !anyNA(an$z_tel))
write.csv(an, "data/analysis_sample.csv", row.names = FALSE)

# knots: 5th / 50th / 95th percentiles of tertiary enrolment in the full panel
KN <- unname(quantile(panel$tert_enrol, c(.05, .5, .95), na.rm = TRUE))
saveRDS(list(KN = KN, enrol_full = panel$tert_enrol), "data/knots.rds")
cat("analysis_sample.csv written:", nrow(an), "rows; knots =", round(KN, 2), "\n")
