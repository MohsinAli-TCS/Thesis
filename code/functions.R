# =============================================================================
# functions.R — shared estimators for
# "Digital Transformation and Poverty Reduction in Low-Human-Capital Economies:
#  Panel Evidence with Pre-Existing Network Instruments"
# Control-function (CF) natural-cubic-spline marginal effects with a
# country-cluster bootstrap that re-estimates both stages.
# =============================================================================
suppressPackageStartupMessages({
  library(fixest); library(dplyr); library(splines); library(ggplot2); library(jsonlite)
})

OUTS     <- c("poverty_inv", "gdp_growth", "emp_ratio", "gini_inv")
OUT_LAB  <- c(poverty_inv = "Poverty (inv.)", gdp_growth = "GDP growth", emp_ratio = "Employment", gini_inv = "Equality (inv.)")
CTRL     <- "log_gdp_pc_l1 + life_exp + tert_edu_exp + findev + inst_composite + trade + fdi + urban"   # stage-2 controls
CTRL1    <- "inst_composite + trade + fdi + urban"                                                      # stage-1 controls
EDU_GRID <- 5:100
PAK_EDU  <- 11            # Pakistan tertiary enrolment, 2021-2023 mean = 11.1%
SEED     <- 20260526      # bootstrap seed used for every 500-rep bootstrap unless stated
R_BOOT   <- 500

# Spline / interaction basis -------------------------------------------------
make_basis <- function(x, spec) {
  if (spec$form == "ns") {
    B <- ns(x, knots = spec$kn_int, Boundary.knots = spec$kn_bd)
    colnames(B) <- paste0("b", seq_len(ncol(B))); B
  } else if (spec$form == "linear") cbind(b1 = x) else cbind(b1 = x, b2 = x^2)
}

# Two-step CF estimator: returns ME at spec$pak_edu and on EDU_GRID for each outcome
#   spec: form ("ns"/"linear"/"quad"), kn_int, kn_bd, cf (TRUE/FALSE), inst (char),
#         country_fe (TRUE/FALSE), pak_edu
fit_cf <- function(dat, spec, outcomes = OUTS, return_models = FALSE) {
  dat$year_f <- factor(dat$year)
  fe <- if (isTRUE(spec$country_fe)) "| year_f + iso3c" else "| year_f"
  if (spec$cf) {
    f1 <- as.formula(paste("dti_score ~", paste(spec$inst, collapse = " + "), "+", CTRL1, fe))
    m1 <- feols(f1, data = dat, vcov = "iid", notes = FALSE, warn = FALSE)
    dat$v_hat <- NA_real_; dat$v_hat[obs(m1)] <- resid(m1)
  }
  dat <- dat[!is.na(dat$tert_enrol), ]
  B <- make_basis(dat$tert_enrol, spec); k <- ncol(B)
  for (j in 1:k) { dat[[paste0("b", j)]] <- B[, j]; dat[[paste0("b", j, "_dti")]] <- B[, j] * dat$dti_score }
  terms <- paste(c("dti_score", paste0("b", 1:k), paste0("b", 1:k, "_dti"), CTRL, if (spec$cf) "v_hat"), collapse = " + ")
  Bp <- make_basis(spec$pak_edu, spec); Bg <- make_basis(EDU_GRID, spec)
  out <- list()
  for (y in outcomes) {
    m <- tryCatch(feols(as.formula(paste(y, "~", terms, fe)), data = dat, vcov = "iid", notes = FALSE, warn = FALSE),
                  error = function(e) NULL)
    if (is.null(m)) { out[[y]] <- list(me = NA, grid = rep(NA, length(EDU_GRID)), n = NA, vhat = NA); next }
    cf <- coef(m); g <- cf[paste0("b", 1:k, "_dti")]; g[is.na(g)] <- 0
    out[[y]] <- list(me = unname(cf["dti_score"] + sum(Bp * g)), grid = as.vector(cf["dti_score"] + Bg %*% g),
                     n = nobs(m), vhat = if (spec$cf) unname(cf["v_hat"]) else NA,
                     model = if (return_models) m else NULL, data = if (return_models) dat else NULL)
  }
  out
}

# Country-cluster bootstrap: resample countries with replacement, relabel duplicates,
# re-estimate BOTH stages; knots held fixed --------------------------------------
boot_cf <- function(dat, spec, R = R_BOOT, seed = SEED, outcomes = OUTS) {
  set.seed(seed); ids <- unique(dat$iso3c); G <- length(ids)
  split_idx <- split(seq_len(nrow(dat)), dat$iso3c)
  me   <- matrix(NA, R, length(outcomes), dimnames = list(NULL, outcomes))
  grid <- array(NA, c(R, length(EDU_GRID), length(outcomes)), dimnames = list(NULL, EDU_GRID, outcomes))
  for (r in 1:R) {
    draw <- sample(ids, G, replace = TRUE)
    bd <- do.call(rbind, lapply(seq_along(draw), function(i) { x <- dat[split_idx[[draw[i]]], ]; x$iso3c <- paste0(draw[i], "_", i); x }))
    f <- fit_cf(bd, spec, outcomes)
    for (y in outcomes) { me[r, y] <- f[[y]]$me; grid[r, , y] <- f[[y]]$grid }
  }
  list(me = me, grid = grid)
}

# Cluster-robust Wald F on excluded instruments in the first stage ------------
first_stage <- function(dat, inst, fe = "year_f + iso3c") {
  dat$year_f <- factor(dat$year)
  m <- feols(as.formula(paste("dti_score ~", paste(inst, collapse = " + "), "+", CTRL1, "|", fe)), data = dat, cluster = ~iso3c, notes = FALSE)
  w <- wald(m, keep = paste0("^(", paste(inst, collapse = "|"), ")$"), print = FALSE)
  list(model = m, F = unname(w$stat), p = unname(w$p))
}

# Hansen J (cluster-robust two-step GMM), linear DTI, country + year FE --------
hansenJ <- function(dat, y, inst) {
  dd <- dat[!is.na(dat$tert_enrol), ]; dd$year_f <- factor(dd$year)
  ok <- complete.cases(dd[, c(y, "dti_score", inst, all.vars(as.formula(paste("~", CTRL))))]); dd <- dd[ok, ]
  f_part <- function(v) { m <- feols(as.formula(paste(v, "~", CTRL, "| year_f + iso3c")), data = dd, notes = FALSE, warn = FALSE)
                          r <- rep(0, nrow(dd)); r[obs(m)] <- resid(m); r }
  yt <- f_part(y); xt <- f_part("dti_score"); Z <- sapply(inst, f_part); g <- dd$iso3c
  Sfun <- function(u) Reduce(`+`, lapply(split(seq_along(g), g), function(i) { s <- t(Z[i, , drop = FALSE]) %*% u[i]; s %*% t(s) }))
  b1 <- solve(t(xt) %*% Z %*% solve(crossprod(Z)) %*% t(Z) %*% xt, t(xt) %*% Z %*% solve(crossprod(Z)) %*% t(Z) %*% yt)
  W  <- solve(Sfun(yt - xt %*% b1))
  b2 <- solve(t(xt) %*% Z %*% W %*% t(Z) %*% xt, t(xt) %*% Z %*% W %*% t(Z) %*% yt); u2 <- yt - xt %*% b2
  gb <- t(Z) %*% u2; J <- drop(t(gb) %*% solve(Sfun(u2)) %*% gb)
  pchisq(J, length(inst) - 1, lower.tail = FALSE)
}

# Horse race: add DTI x (demeaned) lagged log income to the main design --------
fit_cf_income <- function(dat, spec, outcome = "poverty_inv", xint = "log_gdp_pc_l1") {
  dat$year_f <- factor(dat$year)
  f1 <- as.formula(paste("dti_score ~", paste(spec$inst, collapse = " + "), "+", CTRL1, "| year_f + iso3c"))
  m1 <- feols(f1, data = dat, vcov = "iid", notes = FALSE, warn = FALSE); dat$v_hat <- NA_real_; dat$v_hat[obs(m1)] <- resid(m1)
  dat <- dat[!is.na(dat$tert_enrol), ]; B <- make_basis(dat$tert_enrol, spec)
  for (j in 1:2) { dat[[paste0("b", j)]] <- B[, j]; dat[[paste0("b", j, "_dti")]] <- B[, j] * dat$dti_score }
  dat$x_dti <- (dat[[xint]] - mean(dat[[xint]], na.rm = TRUE)) * dat$dti_score
  m <- feols(as.formula(paste(outcome, "~ dti_score + b1 + b2 + b1_dti + b2_dti + x_dti +", CTRL, "+ v_hat | year_f + iso3c")),
             data = dat, vcov = "iid", notes = FALSE, warn = FALSE)
  cf <- coef(m); Bp <- make_basis(c(11, 90), spec); me <- cf["dti_score"] + Bp %*% cf[c("b1_dti", "b2_dti")]
  c(diff_11_90 = unname(me[1] - me[2]), x_dti = unname(cf["x_dti"]))
}

# Helpers ----------------------------------------------------------------------
ci_q <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
f2   <- function(x) formatC(x, format = "f", digits = 2)
NAVY <- "#1F3A5F"; GREY <- "#9A9A9A"; CRIMSON <- "#A51C30"
theme_paper <- theme_minimal(base_size = 9) + theme(
  panel.grid.minor = element_blank(), axis.text = element_text(size = 7), strip.text = element_text(size = 9, hjust = 0),
  legend.text = element_text(size = 8), legend.title = element_blank(), legend.position = "bottom",
  plot.title.position = "plot", plot.title = element_text(hjust = 0, size = 9),
  axis.ticks = element_line(colour = "grey40"), axis.ticks.length = unit(2, "pt"))
