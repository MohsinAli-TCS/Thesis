# =============================================================================
# 04_robustness.R — Table 5, Figure 3, leave-one-country-out jackknife,
# and the DTI x income horse race (Section 5.2). Poverty outcome, main design.
# =============================================================================
source("code/functions.R")
an <- read.csv("data/analysis_sample.csv", stringsAsFactors = FALSE)
kk <- readRDS("data/knots.rds"); MAIN <- readRDS("output/main_spec.rds"); mr <- readRDS("output/main_results.rds")
qk <- function(p) unname(quantile(kk$enrol_full, p, na.rm = TRUE)); spec <- function(...) modifyList(MAIN, list(...))
an_d3 <- an %>% mutate(dti_score = dti3)
specs <- list(
  M0  = list("Main design: telephones + ICT imports, country + year FE", MAIN, an),
  M1  = list("Knots: boundary 1/99, interior 50", spec(kn_int = qk(.5), kn_bd = qk(c(.01, .99))), an),
  M2  = list("Knots: boundary 10/90, interior 50", spec(kn_int = qk(.5), kn_bd = qk(c(.10, .90))), an),
  M3  = list("Knots: 5/95 + interior 33/67 (3 df)", spec(kn_int = qk(c(1/3, 2/3))), an),
  M4  = list("Knots: 5/95 + interior 25/50/75 (4 df)", spec(kn_int = qk(c(.25, .5, .75))), an),
  M5  = list("Linear DTI x enrolment", spec(form = "linear"), an),
  M6  = list("Quadratic DTI x enrolment", spec(form = "quad"), an),
  M7  = list("No control function (FE-OLS)", spec(cf = FALSE), an),
  M8  = list("Emerging economies only", MAIN, an %>% filter(group == "Emerging")),
  M9  = list("Drop Pakistan", MAIN, an %>% filter(iso3c != "PAK")),
  M10 = list("Pre-COVID (2003-2019)", MAIN, an %>% filter(year <= 2019)),
  M11 = list("DTI rebuilt without DFI (3-component PCA)", MAIN, an_d3))
rob <- list()
for (id in names(specs)) { s <- specs[[id]]
  pt <- fit_cf(s[[3]], s[[2]], outcomes = "poverty_inv")$poverty_inv
  me <- if (id == "M0") mr$boot$me[, "poverty_inv"] else boot_cf(s[[3]], s[[2]], outcomes = "poverty_inv")$me[, 1]
  q <- ci_q(me); rob[[id]] <- data.frame(id = id, spec = s[[1]], N = pt$n, ME = pt$me, lo = q[1], hi = q[2], share_le0 = mean(me <= 0, na.rm = TRUE))
  cat(id, "done\n") }
rob <- do.call(rbind, rob)
jk <- do.call(rbind, lapply(sort(unique(an$iso3c)), function(cc) data.frame(dropped = cc, me = fit_cf(an[an$iso3c != cc, ], MAIN, outcomes = "poverty_inv")$poverty_inv$me)))
SIG <- mr$SIGMA
T5 <- rob %>% transmute(Specification = spec, N, ME_per_unit = f2(ME), CI = paste0("[", f2(lo), ", ", f2(hi), "]"),
                        ME_per_SD = f2(SIG * ME), CI_per_SD = paste0("[", f2(SIG * lo), ", ", f2(SIG * hi), "]"), Share_le0 = formatC(share_le0, format = "f", digits = 3))
T5 <- rbind(T5, data.frame(Specification = "Leave-one-country-out jackknife (range of 37 estimates)", N = "", ME_per_unit = paste0(f2(min(jk$me)), " to ", f2(max(jk$me))), CI = "",
                           ME_per_SD = paste0(f2(SIG * min(jk$me)), " to ", f2(SIG * max(jk$me))), CI_per_SD = "", Share_le0 = ""))
write.csv(T5, "output/tables/table5_robustness_poverty.csv", row.names = FALSE); write.csv(jk, "output/tables/jackknife_loco_poverty.csv", row.names = FALSE)

# ---- Horse race: DTI x lagged log income (300 country-cluster reps, seed 31) ----
hr_pt <- fit_cf_income(an, MAIN); set.seed(31); ids <- unique(an$iso3c); spl <- split(seq_len(nrow(an)), an$iso3c)
hr_b <- t(replicate(300, { dr <- sample(ids, length(ids), TRUE)
  bd <- do.call(rbind, lapply(seq_along(dr), function(i) { x <- an[spl[[dr[i]]], ]; x$iso3c <- paste0(dr[i], "_", i); x }))
  tryCatch(fit_cf_income(bd, MAIN), error = function(e) rep(NA, 2)) }))
hr <- data.frame(stat = names(hr_pt), estimate = hr_pt, lo = apply(hr_b, 2, quantile, .025, na.rm = TRUE), hi = apply(hr_b, 2, quantile, .975, na.rm = TRUE),
                 share_le0 = colMeans(hr_b <= 0, na.rm = TRUE))
write.csv(hr, "output/tables/horse_race_income.csv", row.names = FALSE)

# ---- Figure 3 ----------------------------------------------------------------
fr <- rbind(rob %>% transmute(spec, ME, lo, hi), data.frame(spec = "Leave-one-country-out (range of 37 estimates)", ME = NA, lo = min(jk$me), hi = max(jk$me)))
fr$clears <- ifelse(fr$lo > 0, "95% CI excludes zero", "95% CI spans zero"); fr$spec <- factor(fr$spec, levels = rev(fr$spec))
p3 <- ggplot(fr, aes(y = spec)) + geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = .3) +
  geom_vline(xintercept = rob$ME[1], colour = NAVY, linetype = "dotted", linewidth = .4) +
  geom_segment(aes(x = lo, xend = hi, yend = spec, colour = clears), linewidth = .6) + geom_point(aes(x = ME, colour = clears), size = 1.8, na.rm = TRUE) +
  scale_colour_manual(values = c("95% CI excludes zero" = NAVY, "95% CI spans zero" = GREY)) +
  scale_x_continuous(name = "Marginal effect at 11% enrolment, per DTI unit (pp)", sec.axis = sec_axis(~ . * SIG, name = "Per standard deviation of the DTI (pp)")) +
  labs(y = NULL, title = "The poverty effect holds in every robustness check except the small emerging-only sample") + theme_paper + theme(panel.grid.major.y = element_blank())
ggsave("output/figures/fig3_robustness_poverty.png", p3, width = 7, height = 4, dpi = 300, bg = "white")
cat("04 done\n")
