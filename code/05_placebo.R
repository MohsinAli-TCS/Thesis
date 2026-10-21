# =============================================================================
# 05_placebo.R — Appendix Figure A1: thesis (year-FE-only) CF specification with the
# three original instruments replaced by independent N(0,1) noise (500 draws, seed 4242)
# =============================================================================
source("code/functions.R")
an <- read.csv("data/analysis_sample.csv", stringsAsFactors = FALSE)
MAIN <- readRDS("output/main_spec.rds"); dr <- readRDS("output/designs_results.rds")
THREE <- c("rq_est_l1", "ict_imports_l1", "govcapex_l1")
THESIS <- modifyList(MAIN, list(inst = THREE, country_fe = FALSE))
set.seed(4242)
pl <- t(replicate(500, { P <- an; for (v in THREE) P[[v]] <- rnorm(nrow(P)); sapply(fit_cf(P, THESIS), `[[`, "me") }))
write.csv(pl, "output/tables/placebo_noise_instruments.csv", row.names = FALSE)
pld <- do.call(rbind, lapply(OUTS, function(y) data.frame(outcome = factor(OUT_LAB[y], levels = OUT_LAB), me = pl[, y]))) %>%
  group_by(outcome) %>% filter(me >= quantile(me, .01, na.rm = TRUE), me <= quantile(me, .99, na.rm = TRUE)) %>% ungroup()
vl <- do.call(rbind, lapply(OUTS, function(y) data.frame(outcome = factor(OUT_LAB[y], levels = OUT_LAB),
  x = c(dr$T0$pt[[y]]$me, dr$T1$pt[[y]]$me), what = c("Thesis CF estimate (real instruments)", "Naive estimate (no V-hat)"))))
pA <- ggplot(pld, aes(me)) + geom_histogram(bins = 40, fill = "grey80", colour = "white", linewidth = .2) +
  geom_vline(data = vl, aes(xintercept = x, colour = what, linetype = what), linewidth = .7) + facet_wrap(~outcome, scales = "free") +
  scale_colour_manual(values = c("Thesis CF estimate (real instruments)" = NAVY, "Naive estimate (no V-hat)" = GREY)) +
  scale_linetype_manual(values = c("Thesis CF estimate (real instruments)" = "solid", "Naive estimate (no V-hat)" = "dashed")) +
  labs(x = "ME of DTI at 11% enrolment with pure-noise instruments (500 draws; central 98% shown)", y = "Placebo draws",
       title = "Real-instrument CF estimates fall inside the range produced by pure-noise instruments") + theme_paper
ggsave("output/figures/figA1_placebo_instruments.png", pA, width = 7.2, height = 4.6, dpi = 300, bg = "white")
cat("05 done\n")
