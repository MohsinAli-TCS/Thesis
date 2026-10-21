# =============================================================================
# 03_instrument_designs.R — Table 4 and Figure 4 (alternative instrument sets),
# Hansen J tests, and the thesis (year-FE-only) specification used in Fig. 4 / Fig. A1
# =============================================================================
source("code/functions.R")
an <- read.csv("data/analysis_sample.csv", stringsAsFactors = FALSE)
KN <- readRDS("data/knots.rds")$KN; MAIN <- readRDS("output/main_spec.rds"); mr <- readRDS("output/main_results.rds")
spec <- function(...) modifyList(MAIN, list(...))
THREE <- c("rq_est_l1", "ict_imports_l1", "govcapex_l1")      # govcapex = general government TOTAL expenditure (% GDP)
designs <- list(
  E4 = list(lab = "Telephones + ICT imports (main)", inst = c("z_tel", "ict_imports_l1")),
  E2 = list(lab = "Telephones only",                 inst = "z_tel"),
  E3 = list(lab = "Mobile only",                     inst = "z_mob"),
  E5 = list(lab = "ICT imports only",                inst = "ict_imports_l1"),
  E1 = list(lab = "Original three instruments",      inst = THREE))
res <- list()
for (id in names(designs)) { s <- spec(inst = designs[[id]]$inst)
  b <- if (id == "E4") mr$boot else boot_cf(an, s)
  res[[id]] <- list(pt = fit_cf(an, s), me = b$me, F = first_stage(an, designs[[id]]$inst)$F)
  cat(id, "done\n") }
# Thesis specification (weak IVs, year FE only) and naive (no CF), year FE only
THESIS <- spec(inst = THREE, country_fe = FALSE); NAIVE <- spec(cf = FALSE, country_fe = FALSE)
res$T0 <- list(pt = fit_cf(an, THESIS), me = boot_cf(an, THESIS)$me, F = first_stage(an, THREE, "year_f")$F)
res$T1 <- list(pt = fit_cf(an, NAIVE),  me = boot_cf(an, NAIVE)$me,  F = NA)
saveRDS(res, "output/designs_results.rds")

hj <- function(id, y) if (length(designs[[id]]$inst) > 1) hansenJ(an, y, designs[[id]]$inst) else NA
T4 <- do.call(rbind, lapply(names(designs), function(id) do.call(rbind, lapply(OUTS, function(y) { q <- ci_q(res[[id]]$me[, y])
  data.frame(Design = designs[[id]]$lab, F = f2(res[[id]]$F), Outcome = OUT_LAB[y], ME = f2(res[[id]]$pt[[y]]$me),
             CI = paste0("[", f2(q[1]), ", ", f2(q[2]), "]"), Hansen_J_p = ifelse(is.na(hj(id, y)), "-", formatC(hj(id, y), format = "f", digits = 3))) }))))
write.csv(T4, "output/tables/table4_instrument_designs_all_outcomes.csv", row.names = FALSE)

# ---- Figure 4 ----------------------------------------------------------------
labs4 <- c(T0 = "Thesis CF: 3 IVs, year FE", T1 = "No IV (naive), year FE", E1 = "Thesis IVs + country FE",
           E2 = "Pre-1995 telephones x global DTI + country FE", E3 = "Pre-1999 mobile x global DTI + country FE",
           E4 = "Telephones + ICT imports + country FE", E5 = "ICT imports only + country FE")
fd <- do.call(rbind, lapply(names(labs4), function(id) do.call(rbind, lapply(OUTS, function(y) { q <- ci_q(res[[id]]$me[, y])
  data.frame(id = id, row = paste0(labs4[id], ifelse(is.na(res[[id]]$F), "", sprintf("  [F = %.1f]", res[[id]]$F))), outcome = OUT_LAB[y],
             ME = res[[id]]$pt[[y]]$me, lo = q[1], hi = q[2]) }))))
fd <- fd %>% mutate(clears = ifelse(lo > 0 | hi < 0, "95% CI excludes zero", "95% CI spans zero"), outcome = factor(outcome, levels = OUT_LAB),
                    block = factor(ifelse(id %in% c("T0", "T1"), "Thesis specification", "Country + year FE designs"), levels = c("Thesis specification", "Country + year FE designs")))
fd$row <- factor(fd$row, levels = rev(unique(fd$row)))
p4 <- ggplot(fd, aes(y = row)) + geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = .3) +
  geom_segment(aes(x = lo, xend = hi, yend = row, colour = clears), linewidth = .5) + geom_point(aes(x = ME, colour = clears), size = 1.7) +
  facet_grid(block ~ outcome, scales = "free", space = "free_y") + scale_colour_manual(values = c("95% CI excludes zero" = NAVY, "95% CI spans zero" = GREY)) +
  scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
  labs(x = "Marginal effect of DTI at 11% tertiary enrolment (500 country-cluster bootstrap reps)", y = NULL,
       title = "With strong instruments the poverty effect is stable (about 10-12); growth and equality change sign with the instrument") +
  theme_paper + theme(strip.text.y = element_text(angle = 0, size = 7, hjust = 0), panel.grid.major.y = element_blank())
ggsave("output/figures/fig4_instrument_designs.png", p4, width = 7.4, height = 3.8, dpi = 300, bg = "white")
cat("03 done\n")
