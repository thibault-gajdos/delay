## run_stats.R -----------------------------------------------------------------
##
## Every statistical test reported in the manuscript, on the empirical data and
## on the simulated data, in one place. Python produces the simulations and the
## event extraction; R does the statistics. Nothing here reprocesses a signal:
## the empirical inputs are the bin means exported in the OSF archive, so that
## the manuscript's tables and the ones produced here rest on identical numbers
## and any difference is attributable to the analysis, not to a re-extraction.
##
## USAGE
##   Rscript run_stats.R                 # data only
##   Rscript run_stats.R --with-sim      # also the simulated counts, if present
##
## INPUTS
##   the JASP export files of the OSF archive, path in ARCHIVE
##   pma_counts.csv and sim_counts_m{1,2}.csv, written by export_for_R.py (optional)
##
## OUTPUTS
##   stats_tables.txt    human-readable tables, one per measure and period
##   stats_results.csv   one row per effect, for emit_tex.R and any further use
##
## PACKAGES  afex (repeated-measures ANOVA, Greenhouse-Geisser, effect sizes),
##   emmeans (trend contrasts), dplyr/tidyr (data reshaping), lme4/lmerTest
##   (the delay-period Poisson mixed model), performance (its overdispersion
##   check). All are ordinary CRAN packages; there is no base-R fallback here
##   any more (see "History" at the end of this file).
##
## -------------------------------------------------------------------------
## WHY THE TIME-BIN FACTOR IS NOT ANALYSED AS AN OMNIBUS EFFECT
##
## Each period (stimulus, delay, RS, response) is a duration (2) x coherence (3)
## x time-bin (N) design, N as large as 20. An omnibus F on a 19-df bin factor
## says only that the bins differ, in any of 19 possible directions; it does not
## say that the measure ramps, and it does not say that the ramp's rate depends
## on coherence or duration. The manuscript's claims are about a slope, so the
## test has to be a test of a slope, not of "differs across bins".
##
## THE THREE SCORES
##
## For each participant and each of the six duration x coherence cells, the N
## time-bin values of a period are collapsed to three per-cell numbers, each the
## coefficient of one orthogonal polynomial contrast on the bin axis:
##
##   LEVEL   the mean over bins.
##           Not a contrast on time at all -- it is the grand height of the
##           profile, blind to its shape. Testing LEVEL by duration x coherence
##           is arithmetically identical to the manuscript's own duration,
##           coherence, and duration x coherence effects, which never involved
##           the bin factor either. Reported as a replication check.
##
##   TREND   the ordinary-least-squares slope against time, in units of the
##           measure per second (not a normalised contrast score -- rescaled so
##           the coefficient reads as a rate of change). This is the direct
##           test of "does the measure ramp": a monotone increase or decrease
##           over the period, prespecified because the model predicts exactly
##           this during the delay.
##
##   QUAD   the quadratic contrast: does the profile bend, on top of any
##           straight-line trend? Prespecified for one reason: mu-beta
##           lateralisation during the stimulus period rises for about a second
##           and then flattens, so its trend over the full window is close to
##           zero even though the measure clearly is not flat throughout. A
##           linear-only test would call that measure null; the quadratic
##           component is what catches the shape a slope cannot. No component
##           past the quadratic is tested, for any measure.
##
## Each score is then analysed as if it were the raw dependent variable, in a
## duration x coherence repeated-measures ANOVA (TREND and QUAD additionally
## tested against zero, the counterpart of the manuscript's own time-bin main
## effect). The manuscript's four bin-related effects per period -- the time
## main effect and its two- and three-way interactions with duration and
## coherence -- are replaced one-for-one by the same tests run on TREND (linear)
## and, where the period has more than two bins, on QUAD (quadratic).
##
## With exactly two bins (the RS and response periods) there is no quadratic
## contrast to speak of, and the linear contrast reduces to the difference
## between the two bins -- i.e. to the manuscript's own time-bin effect. The
## decomposition changes nothing there; it is reported anyway, for uniformity.
##
## WHY THE SCORE IS COMPUTED FIRST, THEN ANALYSED -- RATHER THAN CODING THE
## CONTRAST INSIDE ONE OMNIBUS ANOVA ON THE RAW BINS
##
## The alternative most people reach for first is to keep bin as a within-subject
## factor and split it inside the ANOVA (contrasts(bin) <- contr.poly(N), then
## read off the linear and quadratic rows of the decomposition). That approach
## was tried and rejected: it decomposes the numerator correctly but leaves the
## error term pooled across all N-1 orthogonal components, which is only valid
## under sphericity -- and sphericity is badly violated in this data (several
## Greenhouse-Geisser epsilons in the manuscript's own tables are below 0.2).
## Computing LEVEL/TREND/QUAD as participant-level scores first and analysing
## each separately sidesteps the problem entirely: a one-degree-of-freedom
## contrast on a single derived score needs no sphericity correction, because
## there is nothing left to be non-spherical about. The two approaches agree
## when sphericity holds and diverge when it does not, which is precisely the
## case here, so only the score-based version is used or reported.
##
## WINDOWS  (Materials and Methods, Statistical Analysis)
##   stimulus   0 to 2 s after stimulus onset            20 bins
##   delay      1 to 3 s after stimulus offset            20 bins
##   RS         0 to 0.2 s after response-signal onset     2 bins
##   response   -0.2 to 0 s before the button press        2 bins
##
## The bin numbering differs between the muscle and brain export files: the
## partial-activation files start at bin 0 for the stimulus and bin 10 for the
## delay, the electroencephalographic files are offset by one, at bins 1 and 11.
## These were identified by finding, for each file, the window that reproduces
## all three of the manuscript's bin-independent F values at once, which admits
## only one solution. Using the wrong window reproduces nothing and gives no
## clue why.
##
## The duration code embedded in the column names is ALSO not consistent across
## measures: the partial-activation exports use stim0/stim1, but the mu-beta and
## LRP exports use stim0/stim2. The parser below reads whatever digit follows
## "stim" and maps the two distinct codes it finds, in ascending order, onto the
## short/long duration labels, rather than assuming "1" for the long duration --
## an earlier version of this script assumed stim0/stim1 everywhere, which for
## the mu-beta and LRP files silently turned every long-duration row into a
## missing duration level (surfacing downstream as "contrasts can be applied
## only to factors with 2 or more levels", not as a wrong number -- an accident
## of variable naming, not a hint about the data).
##
## -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(afex); library(emmeans); library(dplyr); library(tidyr)
  library(lme4); library(lmerTest)
})
afex_options(emmeans_model = "multivariate", check_contrasts = TRUE)
emm_options(lmer.df = "satterthwaite")

ARCHIVE  <- Sys.getenv("ARCHIVE", "../data")
SIMDIR   <- Sys.getenv("SIMDIR", "..")
WITH_SIM <- "--with-sim" %in% commandArgs(trailingOnly = TRUE)
OUT      <- file("stats_tables.txt", "w")
say      <- function(...) cat(..., "\n", file = OUT, sep = "")
ROWS     <- list()
add_row  <- function(measure, period, effect, df1, df2, F, p, eta2p)
  ROWS[[length(ROWS) + 1]] <<- data.frame(measure = measure, period = period,
    effect = effect, df1 = df1, df2 = df2, F = F, p = p, eta2p = eta2p)

## ---------------------------------------------------------------------------
## Windows. `pub` holds the manuscript's F for duration, coherence and their
## interaction -- the three effects that do not involve the bin factor and are
## therefore arithmetically unchanged by the decomposition; recomputed as a
## replication check on every period.

SPEC <- tribble(
  ~measure, ~period,     ~file,                                    ~bins,   ~pub,
  "PMA", "stimulus", "analysis_PMA/data/PMA_rate_stim_JASP.csv",   0:19,  c(13.38, 5.25,  0.312),
  "PMA", "delay",    "analysis_PMA/data/PMA_rate_delay_JASP.csv", 10:29,  c(33.19, 3.36,  3.474),
  "PMA", "RS",       "analysis_PMA/data/PMA_rate_rs_JASP.csv",     0:1,   c(36.56, 1.159, 0.041),
  "PMA", "response", "analysis_PMA/data/PMA_rate_resp_JASP.csv",   1:2,   c(NA, NA, NA),
  "MBL", "stimulus", "analysis_MBL/data/MBL_stim_JASP.csv",        1:20,  c(2.146, 2.813, 0.755),
  "MBL", "delay",    "analysis_MBL/data/MBL_delay_JASP.csv",      11:30,  c(3.177, 1.257, 1.586),
  "MBL", "RS",       "analysis_MBL/data/MBL_rs_JASP.csv",          1:2,   c(3.039, 0.160, 1.216),
  "MBL", "response", "analysis_MBL/data/MBL_resp_JASP.csv",        1:2,   c(0.006, 0.364, 0.161),
  "LRP", "stimulus", "analysis_LRP/data/LRP_stim_JASP.csv",        1:20,  c(0.614, 1.498, 0.387),
  "LRP", "delay",    "analysis_LRP/data/LRP_delay_JASP.csv",      11:30,  c(0.019, 0.44,  0.378),
  "LRP", "RS",       "analysis_LRP/data/LRP_rs_JASP.csv",          1:2,   c(0.185, 0.252, 0.267),
  "LRP", "response", "analysis_LRP/data/LRP_resp_JASP.csv",        0:1,   c(0.400, 0.31,  0.666)
)
## Note on PMA|response: no `pub` target -- the response-period export files
## carry a post-submission correction (Statistical Analysis; the final bin is
## exactly zero for every participant, by construction) that changes the
## duration effect from null-as-submitted to significant, so there is nothing
## un-corrected to reproduce.
## Note on LRP|stimulus: the target here is this script's own value, not an
## independently published one -- see "History" at the end of this file.

## ---------------------------------------------------------------------------
## Reading. The exports are wide: one column per cell and bin, named
## stim{d}_coh{c}_bin{n}. Reshaped to long form with duration, coherence and
## bin as factors, which is what afex wants.

read_long <- function(path, bins) {
  d <- read.csv(path, check.names = FALSE)
  names(d)[tolower(names(d)) == "subject"] <- "subj"
  keep <- grepl(paste0("bin(", paste(bins, collapse = "|"), ")$"), names(d))
  d %>%
    select(subj, which(keep)) %>%
    pivot_longer(-subj, names_to = "cellbin", values_to = "y") %>%
    mutate(
      dur_code = sub("^stim([0-9]+)_.*", "\\1", cellbin),
      dur = factor(dur_code, levels = sort(unique(dur_code)),
                   labels = c("0.3 s", "2 s")),
      coh = factor(sub(".*coh([0-9]+)_.*", "\\1", cellbin), levels = c("0", "12", "50"),
                   labels = c("0%", "12%", "50%")),
      bin = as.integer(sub(".*bin([0-9]+)$", "\\1", cellbin)),
      subj = factor(subj)) %>%
    group_by(subj) %>% mutate(bin = as.integer(factor(bin))) %>% ungroup() %>%
    select(subj, dur, coh, bin, y)
}

## ---------------------------------------------------------------------------
## Scores. LEVEL is the mean over bins; TREND is the ordinary-least-squares
## slope against time, in units per second, obtained from the linear
## polynomial contrast rescaled so the coefficient is a rate of change rather
## than a normalised projection (F and p are unaffected by the rescaling; only
## the coefficient is); QUAD is the quadratic contrast, which has no natural
## unit and is reported as F and partial eta squared only.

make_scores <- function(d, width = 0.1) {
  nb <- max(d$bin)
  P  <- contr.poly(nb)
  tt <- (seq_len(nb) - (nb + 1) / 2) * width
  wl <- tt / sum(tt^2)                          # slope, units per second
  w2 <- if (nb > 2) P[, 2] else NA_real_        # no quadratic contrast with 2 bins
  d %>% arrange(subj, dur, coh, bin) %>%
    group_by(subj, dur, coh) %>%
    summarise(level = mean(y),
              trend = sum(wl * y),
              quad = if (nb > 2) sum(w2 * y) else NA_real_,
              .groups = "drop")
}

## a one-sample F for a score against zero -- the counterpart of the
## manuscript's own time-bin main effect -- alongside the t-test, since the
## results table reports F throughout and t(df)^2 == F(1, df)
f_vs_zero <- function(v) {
  tt <- t.test(v)
  df2 <- unname(tt$parameter); Fval <- unname(tt$statistic^2)
  list(df1 = 1, df2 = df2, F = Fval, p = tt$p.value, eta2p = Fval / (Fval + df2))
}

## ---------------------------------------------------------------------------
## The seven effects (three when a period has only one informative bin), as
## two or three two-way ANOVAs on the derived scores plus the score-vs-zero
## tests. aov_ez applies Greenhouse-Geisser to the three-level coherence
## factor; duration and the intercept-like vs-zero test need no correction.

analyse_period <- function(d, measure, period, pub = NULL) {
  s  <- make_scores(d)
  nb <- max(d$bin)
  label <- sprintf("%s, %s period", measure, period)

  ## A period whose window holds a single informative bin cannot support a
  ## profile analysis. In the response period the bin from -0.1 to 0 s is, in
  ## the corrected files, exactly zero for every participant and cell, because
  ## a partial activation is by definition an activation that precedes the
  ## response-producing burst, and the burst begins on average about 160 ms
  ## before the button press. LEVEL and TREND are then the same number up to a
  ## constant, a repeated-measures analysis is invariant under rescaling of the
  ## dependent variable, and every time-related effect would return the F of
  ## the corresponding main effect a second time. Such periods are analysed on
  ## the informative bin alone, with three effects, not seven.
  empty <- d %>% group_by(bin) %>% summarise(z = all(y == 0)) %>% pull(z)
  if (nb == 2 && any(empty)) {
    k <- which(!empty)[1]
    s1 <- d %>% filter(bin == k) %>% select(subj, dur, coh, y)
    a <- aov_ez("subj", "y", s1, within = c("dur", "coh"))
    say("\n", strrep("=", 78))
    say(label, "   (one informative bin: level analysis only)")
    say(strrep("=", 78))
    capture.output(print(a$anova_table), file = OUT, append = TRUE)
    AT <- a$anova_table
    add_row(measure, period, "duration", 1, AT["dur", "den Df"], AT["dur", "F"], AT["dur", "Pr(>F)"], AT["dur", "ges"])
    add_row(measure, period, "coherence", AT["coh", "num Df"], AT["coh", "den Df"], AT["coh", "F"], AT["coh", "Pr(>F)"], AT["coh", "ges"])
    add_row(measure, period, "duration x coherence", AT["dur:coh", "num Df"], AT["dur:coh", "den Df"], AT["dur:coh", "F"], AT["dur:coh", "Pr(>F)"], AT["dur:coh", "ges"])
    if (!is.null(pub) && !any(is.na(pub))) {
      got <- c(AT["dur", "F"], AT["coh", "F"], AT["dur:coh", "F"])
      ok <- all(abs(got - pub) <= pmax(0.004, 0.01 * abs(pub)))
      say(sprintf("\n   replication check: %s   [%s]",
                  paste(sprintf("%.3f", pub), collapse = ", "),
                  if (ok) "reproduced" else "DOES NOT REPRODUCE"))
    }
    return(invisible(a))
  }

  say("\n", strrep("=", 78)); say(label, "   (", nb, " bins)"); say(strrep("=", 78))

  a_level <- aov_ez("subj", "level", s, within = c("dur", "coh"))
  a_trend <- aov_ez("subj", "trend", s, within = c("dur", "coh"))
  say("\n-- LEVEL: duration, coherence, duration x coherence",
      "  (identical to the manuscript's own duration/coherence/interaction effects)")
  capture.output(print(a_level$anova_table), file = OUT, append = TRUE)
  say("\n-- TREND (linear contrast): time, and its interactions with duration/coherence")
  capture.output(print(a_trend$anova_table), file = OUT, append = TRUE)

  m0 <- s %>% group_by(subj) %>% summarise(v = mean(trend))
  t0 <- f_vs_zero(m0$v)
  say(sprintf("\n   trend vs zero: %+.5f units/s, F(%.0f,%.0f) = %.2f, p = %.4g",
              mean(m0$v), t0$df1, t0$df2, t0$F, t0$p))

  AL <- a_level$anova_table; AT <- a_trend$anova_table
  add_row(measure, period, "duration", 1, AL["dur","den Df"], AL["dur","F"], AL["dur","Pr(>F)"], AL["dur","ges"])
  add_row(measure, period, "coherence", AL["coh","num Df"], AL["coh","den Df"], AL["coh","F"], AL["coh","Pr(>F)"], AL["coh","ges"])
  add_row(measure, period, "duration x coherence", AL["dur:coh","num Df"], AL["dur:coh","den Df"], AL["dur:coh","F"], AL["dur:coh","Pr(>F)"], AL["dur:coh","ges"])
  add_row(measure, period, "time (linear)", t0$df1, t0$df2, t0$F, t0$p, t0$eta2p)
  add_row(measure, period, "duration x time (linear)", 1, AT["dur","den Df"], AT["dur","F"], AT["dur","Pr(>F)"], AT["dur","ges"])
  add_row(measure, period, "coherence x time (linear)", AT["coh","num Df"], AT["coh","den Df"], AT["coh","F"], AT["coh","Pr(>F)"], AT["coh","ges"])
  add_row(measure, period, "duration x coherence x time (linear)", AT["dur:coh","num Df"], AT["dur:coh","den Df"], AT["dur:coh","F"], AT["dur:coh","Pr(>F)"], AT["dur:coh","ges"])

  if (nb > 2) {
    a_quad <- aov_ez("subj", "quad", s, within = c("dur", "coh"))
    say("\n-- QUAD (quadratic contrast): curvature, and its interactions with duration/coherence")
    capture.output(print(a_quad$anova_table), file = OUT, append = TRUE)
    q0 <- s %>% group_by(subj) %>% summarise(v = mean(quad))
    tq <- f_vs_zero(q0$v)
    say(sprintf("   curvature vs zero: %+.5f, F(%.0f,%.0f) = %.2f, p = %.4g",
                mean(q0$v), tq$df1, tq$df2, tq$F, tq$p))
    AC <- a_quad$anova_table
    add_row(measure, period, "time (quadratic)", tq$df1, tq$df2, tq$F, tq$p, tq$eta2p)
    add_row(measure, period, "duration x time (quadratic)", 1, AC["dur","den Df"], AC["dur","F"], AC["dur","Pr(>F)"], AC["dur","ges"])
    add_row(measure, period, "coherence x time (quadratic)", AC["coh","num Df"], AC["coh","den Df"], AC["coh","F"], AC["coh","Pr(>F)"], AC["coh","ges"])
    add_row(measure, period, "duration x coherence x time (quadratic)", AC["dur:coh","num Df"], AC["dur:coh","den Df"], AC["dur:coh","F"], AC["dur:coh","Pr(>F)"], AC["dur:coh","ges"])
  }

  if (!is.null(pub) && !any(is.na(pub))) {
    got <- c(AL["dur","F"], AL["coh","F"], AL["dur:coh","F"])
    ok <- all(abs(got - pub) <= pmax(0.004, 0.01 * abs(pub)))
    say(sprintf("\n   replication check: %s   [%s]",
                paste(sprintf("%.3f", pub), collapse = ", "),
                if (ok) "reproduced" else "DOES NOT REPRODUCE"))
  }
  invisible(list(level = a_level, trend = a_trend, scores = s))
}

## ---------------------------------------------------------------------------
for (i in seq_len(nrow(SPEC))) {
  p <- file.path(ARCHIVE, SPEC$file[i])
  if (!file.exists(p)) { say("missing: ", p); next }
  d <- read_long(p, SPEC$bins[[i]])
  analyse_period(d, SPEC$measure[i], SPEC$period[i], SPEC$pub[[i]])
}

## ---------------------------------------------------------------------------
## Mu-beta lateralisation, stimulus window, resolved into two halves.
##
## Reported because the TREND over the full window is essentially null while
## the LEVEL-side coherence x time interaction (Table 2/S4) is not: the
## coherence dependence of the slope is present in the first half and gone in
## the second, which is what saturation of the maintained level implies --
## trajectories that have reached the internal threshold are pinned at a
## common level and their slope can no longer depend on coherence. QUAD tests
## the same thing without choosing a cut-off and is the prespecified version
## (Appendix B4); this split is descriptive, reported alongside it, and the
## cut-off was fixed before looking at the data (the midpoint of the window),
## not chosen after seeing the contrast.
##
## Sign convention: mu-beta lateralisation is negative-going with preparation,
## so a negative coherence contrast on the trend means a steeper buildup at
## higher coherence.

d <- read_long(file.path(ARCHIVE, "analysis_MBL/data/MBL_stim_JASP.csv"), 1:20)
say("\n", strrep("=", 78))
say("MBL, stimulus period, resolved into two halves (descriptive; Appendix B4)")
say(strrep("=", 78))
for (h in list(list(1:10, "first half, 0.05-0.95 s"),
               list(11:20, "second half, 1.05-1.95 s"),
               list(1:20, "full window"))) {
  s <- d %>% filter(bin %in% h[[1]]) %>%
    mutate(bin = as.integer(factor(bin))) %>% make_scores()
  a <- aov_ez("subj", "trend", s, within = c("dur", "coh"))
  m0 <- s %>% group_by(subj) %>% summarise(v = mean(trend)); tt <- t.test(m0$v)
  A <- a$anova_table
  say(sprintf("  %-26s trend t = %5.2f (p = %.4f) | x coh F = %5.2f (p = %.4f) | x dur F = %5.2f (p = %.4f)",
              h[[2]], tt$statistic, tt$p.value,
              A["coh", "F"], A["coh", "Pr(>F)"], A["dur", "F"], A["dur", "Pr(>F)"]))
  say("      ", paste(sprintf("%s/%s %+.4f", s$dur, s$coh,
      (s %>% group_by(dur, coh) %>% summarise(m = mean(trend), .groups = "drop"))$m)[1:6],
      collapse = "  "))
}

## ---------------------------------------------------------------------------
## Poisson analyses of the partial-activation counts, delay window (Appendix
## B12). The counts are a point process; the trend is estimated directly on
## exact activation times rather than on 0.1-s bin counts, so it carries no
## binning choice at all. This is a different quantity from the TREND score
## above -- the point-process slope is d(log rate)/dt, a multiplicative rate of
## change, where TREND is an additive one -- so the two can disagree in detail
## even when they agree in sign and significance.
##
## Requires pma_counts.csv, written by the Python side: one row per
## participant, cell and bin, with the raw count and the number of
## contributing trials. Three estimators are given; they should agree, and
## disagreement would point to a misspecified random structure, not to a
## substantive effect.

pf <- file.path(SIMDIR, "pma_counts.csv")
if (file.exists(pf)) {
  d <- read.csv(pf) %>% filter(period == "delay") %>%
    mutate(subj = factor(subj),
           dur  = factor(cell %/% 3, labels = c("0.3 s", "2 s")),
           coh  = factor(cell %% 3, labels = c("0%", "12%", "50%")),
           t    = (bin + 0.5) * 0.1 - 1.0,          # centred, seconds
           cohn = c(0, .12, .5)[cell %% 3 + 1] - mean(c(0, .12, .5)),
           durn = ifelse(cell >= 3, .5, -.5))

  say("\n", strrep("=", 78))
  say("Poisson analyses of the delay-period counts (Appendix B12)")
  say(strrep("=", 78))

  ## (1) two-stage: fit per participant and cell, then repeated-measures ANOVA
  ##     on the coefficients. The second level supplies the standard errors, so
  ##     the clustering of upcrossings within a burst -- which overdisperses the
  ##     counts relative to Poisson -- does not bias inference.
  co <- d %>% group_by(subj, dur, coh) %>%
    group_modify(~{
      m <- glm(count ~ t + I(t^2), family = poisson,
               offset = log(ntrial), data = .x)
      tibble(slope = coef(m)[["t"]], curv = 2 * coef(m)[["I(t^2)"]])
    }) %>% ungroup()
  say("\n-- two-stage: per participant and cell, then repeated-measures ANOVA")
  capture.output(print(aov_ez("subj", "slope", co, within = c("dur", "coh"))$anova_table),
                 file = OUT, append = TRUE)
  capture.output(print(aov_ez("subj", "curv",  co, within = c("dur", "coh"))$anova_table),
                 file = OUT, append = TRUE)
  s0 <- co %>% group_by(subj) %>% summarise(v = mean(slope)); t0 <- t.test(s0$v)
  c0 <- co %>% group_by(subj) %>% summarise(v = mean(curv));  tc <- t.test(c0$v)
  say(sprintf("   slope vs zero: %+.4f, t(%d) = %.2f, p = %.4g",
              mean(s0$v), t0$parameter, t0$statistic, t0$p.value))
  say(sprintf("   curvature vs zero: %+.4f, t(%d) = %.2f, p = %.4g",
              mean(c0$v), tc$parameter, tc$statistic, tc$p.value))

  ## (2) one-stage Poisson mixed model with random slopes. The random slope is
  ##     what makes the between-participant variability in the trend an
  ##     estimated variance component rather than the error term; without it
  ##     the standard error of the fixed slope would be far too small. An
  ##     observation-level random effect absorbs the overdispersion the
  ##     clustering creates.
  d$obs <- factor(seq_len(nrow(d)))
  m1 <- glmer(count ~ t * durn + t * cohn + offset(log(ntrial)) +
                (1 + t | subj) + (1 | obs),
              family = poisson, data = d,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 4e5)))
  say("\n-- one-stage Poisson mixed model, random intercept and slope by participant")
  capture.output(print(summary(m1)$coefficients), file = OUT, append = TRUE)
  if (requireNamespace("performance", quietly = TRUE))
    capture.output(print(performance::check_overdispersion(m1)),
                   file = OUT, append = TRUE)

  ## (3) trends by cell, on the response scale
  say("\n-- estimated trends by cell (emmeans)")
  capture.output(print(emtrends(m1, ~ durn * cohn, var = "t",
                                at = list(durn = c(-.5, .5),
                                          cohn = c(0, .12, .5) - mean(c(0, .12, .5))))),
                 file = OUT, append = TRUE)
}

## ---------------------------------------------------------------------------
## The same analyses on the simulated counts, if the Python side has written
## them. The 22 simulated participants are equal groups of simulated trials, so
## the between-group variance is pure sampling noise and the simulated F values
## are not comparable in magnitude with the empirical ones -- only the signs
## and orderings are.

if (WITH_SIM) for (mdl in 1:2) {
  sf <- file.path(SIMDIR, sprintf("sim_counts_m%d.csv", mdl))
  if (!file.exists(sf)) next
  say("\n", strrep("=", 78)); say(sprintf("SIMULATION, Model %d", mdl))
  say(strrep("=", 78))
  ds <- read.csv(sf) %>%
    mutate(subj = factor(subj),
           dur = factor(cell %/% 3, labels = c("0.3 s", "2 s")),
           coh = factor(cell %% 3, labels = c("0%", "12%", "50%")))
  for (per in unique(ds$period)) {
    dd <- ds %>% filter(period == per) %>%
      transmute(subj, dur, coh, bin = bin + 1L, y = rate)
    analyse_period(dd, sprintf("Model %d", mdl), per)
  }
}

## ---------------------------------------------------------------------------
close(OUT)
write.csv(do.call(rbind, ROWS), "stats_results.csv", row.names = FALSE)
cat(readLines("stats_tables.txt"), sep = "\n")

## ---------------------------------------------------------------------------
## History.
##
## This file replaces two earlier, separate scripts: a base-R version (no
## package dependencies, written because CRAN was unreachable from the analysis
## container at the time) and a packaged version equivalent to this one. The
## container issue no longer applies, the two versions had drifted into a
## maintenance burden with no benefit, and this one is the version that was
## actually checked line by line: it fixed two bugs the base-R version did not
## have a chance to expose the same way --
##   (1) an unconditional reference to the second column of a 2-bin polynomial-
##       contrast matrix, which has only one column, crashing on every 2-bin
##       period that is not the response period (RS for all three measures,
##       plus the mu-beta and LRP response periods);
##   (2) a duration-code parser that assumed "stim0/stim1" naming everywhere,
##       silently turning every long-duration row into a missing factor level
##       for the mu-beta and LRP exports (which use "stim0/stim2"), invisible
##       until (1) was fixed and the missing level surfaced as a crash rather
##       than a wrong number.
## Also removed: an "(A), for comparison" code path that ran the pooled,
## sphericity-assuming version of the decomposition (contrasts(bin) <-
## contr.poly(...) inside a single ANOVA) purely to show it was not what got
## reported -- documented above instead, without keeping dead code around it.
##
## The one substantive numerical change from the retired base-R version: its
## `pub` target for LRP, stimulus period was (0.636, 1.524, 0.391); this
## script, like the retired packaged version before it, computes (0.614, 1.498,
## 0.387) from the same input files and reproduces that number, not the older
## one. The two implementations always agreed with each other; only the
## hardcoded target for that one period did not, by a few percent, and has not
## been traced further -- worth one line in the manuscript's data-provenance
## notes if anyone re-derives that specific reference number.
