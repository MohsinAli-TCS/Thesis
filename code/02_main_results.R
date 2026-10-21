# =============================================================================
# 02_main_results.R — Tables 1-3, Figures 1-2 (main design)
# Main design: CF spline, instruments z_tel + ict_imports_l1, country + year FE
# =============================================================================
source("code/functions.R")
an <- read.csv("data/analysis_sample.csv", stringsAsFactors = FALSE)
KN <- readRDS("data/knots.rds")$KN
MAIN <- list(form = "ns", kn_int = KN[2], kn_bd = KN[c(1, 3)], cf = TRUE, inst = c("z_tel", "ict_imports_l1"),
             country_fe = TRUE, pak_edu = PAK_EDU)
saveRDS(MAIN, "output/main_spec.rds")

# ---- Table 1: estimation samples --------------------------------------------
fm <- fit_cf(an, MAIN, return_models = TRUE)
T1 <- do.call(rbind, lapply(OUTS, function(y) { d <- fm[[y]]$data[obs(fm[[y]]$model), ]
  data.frame(Outcome = OUT_LAB[y], N = nrow(d), Countries = n_distinct(d$iso3c), Years = paste(range(d$year), collapse = "-"),
             Outcome_mean = f2(mean(d[[y]])), Outcome_SD = f2(sd(d[[y]])), DTI_SD = f2(sd(d$dti_score)),
             LowEnrol_obs_countries = paste0(sum(d$tert_enrol < 20), " (", n_distinct(d$iso3c[d$tert_enrol < 20]), ")")) }))
write.csv(T1, "output/tables/table1_samples.csv", row.names = FALSE)
SIGMA <- sd(fm$poverty_inv$data$dti_score[obs(fm$poverty_inv$model)])   # 1.4661: SD of DTI in the poverty sample

# ---- Table 2: first stage ------------------------------------------------------
fsr <- function(lab, z, fe = "year_f + iso3c") { f <- first_stage(an, z, fe); ct <- coeftable(f$model)
  g <- function(v) if (v %in% rownames(ct)) sprintf("%.4f (%.4f)", ct[v, 1], ct[v, 2]) else ""
  data.frame(Design = lab, Telephones_x_globalDTI = g("z_tel"), Mobile_x_globalDTI = g("z_mob"), ICT_imports_l1 = g("ict_imports_l1"),
             Reg_quality_l1 = g("rq_est_l1"), Govt_total_expenditure_l1 = g("govcapex_l1"), Cluster_F = f2(f$F), N = nobs(f$model)) }
T2 <- rbind(fsr("Telephones + ICT imports (main)", c("z_tel", "ict_imports_l1")), fsr("Telephones only", "z_tel"), fsr("Mobile only", "z_mob"),
            fsr("ICT imports only", "ict_imports_l1"), fsr("Original three instruments", c("rq_est_l1", "ict_imports_l1", "govcapex_l1")),
            fsr("Original three, year FE only", c("rq_est_l1", "ict_imports_l1", "govcapex_l1"), "year_f"))
write.csv(T2, "output/tables/table2_first_stage.csv", row.names = FALSE)

# ---- Main-design bootstrap (all outcomes) -> Table 3 -------------------------
bm <- boot_cf(an, MAIN)
G <- bm$grid[, , "poverty_inv"]; ix <- function(e) which(EDU_GRID == e); pt <- fm$poverty_inv$grid
rows <- t(sapply(c(11, 30, 60, 90), function(e) c(pt[ix(e)], ci_q(G[, ix(e)]))))
dd <- G[, ix(11)] - G[, ix(90)]; rows <- rbind(rows, c(pt[ix(11)] - pt[ix(90)], ci_q(dd)))
T3 <- data.frame(Enrolment = c("11 (Pakistan)", "30", "60", "90", "Difference 11 vs 90"),
                 ME_per_unit = f2(rows[, 1]), CI_per_unit = paste0("[", f2(rows[, 2]), ", ", f2(rows[, 3]), "]"),
                 ME_per_SD = f2(SIGMA * rows[, 1]), CI_per_SD = paste0("[", f2(SIGMA * rows[, 2]), ", ", f2(SIGMA * rows[, 3]), "]"))
write.csv(T3, "output/tables/table3_poverty_me_by_enrolment.csv", row.names = FALSE)
vhat_p <- sapply(OUTS, function(y) summary(fm[[y]]$model, cluster = ~iso3c)$coeftable["v_hat", 4])
saveRDS(list(fm_pt = lapply(fm, function(x) x[c("me", "grid", "n", "vhat")]), boot = bm, SIGMA = SIGMA, vhat_p = vhat_p), "output/main_results.rds")

# ---- Figure 1: first-stage binned scatters (FWL) ---------------------------------
an$year_f <- factor(an$year); m1 <- first_stage(an, c("z_tel", "ict_imports_l1"))$model; ct <- coeftable(m1)
res <- function(v, others) { m <- feols(as.formula(paste(v, "~", paste(c(others, CTRL1), collapse = " + "), "| year_f + iso3c")), data = an, notes = FALSE)
  r <- rep(NA, nrow(an)); r[obs(m)] <- resid(m); r }
mk <- function(z, other, lab) { d <- data.frame(x = res(z, other), y = res("dti_score", other)) %>% filter(!is.na(x), !is.na(y))
  list(raw = d %>% mutate(panel = lab), bins = d %>% mutate(bin = ntile(x, 20)) %>% group_by(bin) %>% summarise(x = mean(x), y = mean(y)) %>% mutate(panel = lab)) }
a <- mk("z_tel", "ict_imports_l1", "Pre-1995 telephone lines x global DTI"); b <- mk("ict_imports_l1", "z_tel", "ICT goods imports (t-1)")
lab <- data.frame(panel = c(a$raw$panel[1], b$raw$panel[1]), txt = c(sprintf("slope = %.4f (SE %.4f)", ct["z_tel", 1], ct["z_tel", 2]),
                                                                     sprintf("slope = %.3f (SE %.3f)", ct["ict_imports_l1", 1], ct["ict_imports_l1", 2])))
p1 <- ggplot(rbind(a$raw, b$raw), aes(x, y)) + geom_point(colour = "grey75", size = .6, alpha = .6) +
  geom_smooth(method = "lm", se = FALSE, colour = NAVY, linewidth = .7, formula = y ~ x) + geom_point(data = rbind(a$bins, b$bins), colour = NAVY, size = 1.8) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~panel, scales = "free_x") + labs(x = "Instrument (residualised)", y = "DTI (residualised)",
  title = sprintf("Both instruments predict the DTI within countries (joint first-stage F = %.1f)", first_stage(an, c("z_tel", "ict_imports_l1"))$F)) + theme_paper
ggsave("output/figures/fig1_first_stage.png", p1, width = 7, height = 3.4, dpi = 300, bg = "white")

# ---- Figure 2: poverty ME curve --------------------------------------------------
c4 <- data.frame(edu = EDU_GRID, me = pt, lo = apply(G, 2, quantile, .025, na.rm = TRUE), hi = apply(G, 2, quantile, .975, na.rm = TRUE))
rug <- data.frame(edu = fm$poverty_inv$data$tert_enrol[obs(fm$poverty_inv$model)])
p2 <- ggplot(c4, aes(edu)) + geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = .3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#C9D3E0") + geom_line(aes(y = me), colour = NAVY, linewidth = .9) +
  geom_vline(xintercept = PAK_EDU, colour = CRIMSON, linetype = "dotted", linewidth = .5) +
  annotate("text", x = 13, y = max(c4$hi), label = "Pakistan (11%)", colour = CRIMSON, size = 2.8, hjust = 0, vjust = 1) +
  geom_rug(data = rug, aes(x = edu), sides = "b", alpha = .3, length = unit(3, "pt"), inherit.aes = FALSE) +
  labs(x = "Gross tertiary enrolment (%)", y = "Marginal effect of DTI (pp)", title = "The poverty effect of digital transformation is largest at low tertiary enrolment") +
  theme_paper + theme(panel.grid.major.y = element_line(colour = "grey92"))
ggsave("output/figures/fig2_poverty_me_curve.png", p2, width = 5.5, height = 3.6, dpi = 300, bg = "white")
cat("02 done. ME(11) poverty =", f2(pt[ix(11)]), " CI", f2(ci_q(G[, ix(11)])), "\n")
