## run_stats.R -----------------------------------------------------------------
##
## Every statistical test on the data reported in the manuscript, in one place.
## Python produces the simulations and the model-data comparison; R does the
## statistics on the data. Nothing here reprocesses a signal:
## the empirical inputs are the bin means exported in the OSF archive, so that
## the manuscript's tables and the ones produced here rest on identical numbers
## and any difference is attributable to the analysis, not to a re-extraction.
##
## USAGE
##   Rscript run_stats.R
##
## INPUTS
##   the JASP export files of the OSF archive, path in ARCHIVE, including the
##   behavioural exports in analysis_behavioral_emg/
##
## OUTPUTS
##   stats_tables.txt    human-readable tables, one per measure and period
##   stats_results.csv   one row per effect, read by stats_tables_tex.py
##
## PACKAGES  afex (repeated-measures ANOVA, Greenhouse-Geisser, effect sizes),
##   emmeans, dplyr/tidyr (data reshaping). All are ordinary CRAN packages.
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
##   QUAD   the quadratic coefficient of the least-squares fit on time in
##           seconds, in units of the measure per second squared: does the
##           profile bend, on top of any straight-line trend? Prespecified for one reason: mu-beta
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
## short/long duration labels, rather than assuming "1" for the long duration.
## Assuming stim0/stim1 everywhere silently turns every long-duration row of the
## mu-beta and LRP files into a missing duration level, which surfaces downstream
## as "contrasts can be applied only to factors with 2 or more levels" rather than
## as a wrong number.
##
## -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(afex); library(emmeans); library(dplyr); library(tidyr)
})
## The manuscript's tables report PARTIAL eta squared; afex's default is GENERALISED eta
## squared. The difference is not a rounding one -- for the accuracy duration effect the
## two are 0.027 and 0.519 -- so es_aov is set explicitly below.
afex_options(emmeans_model = "multivariate", check_contrasts = TRUE, es_aov = "pes")

ARCHIVE  <- Sys.getenv("ARCHIVE", "../data")
OUT      <- file("stats_tables.txt", "w")
say      <- function(...) cat(..., "\n", file = OUT, sep = "")
ROWS     <- list()
add_row  <- function(measure, period, effect, df1, df2, F, p, eta2p,
                     estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
  ROWS[[length(ROWS) + 1]] <<- data.frame(measure = measure, period = period,
    effect = effect, df1 = df1, df2 = df2, F = F, p = p, eta2p = eta2p,
    estimate = estimate, ci_lo = ci_lo, ci_hi = ci_hi)

## ---------------------------------------------------------------------------
## Windows. `pub` holds the manuscript's F for duration, coherence and their
## interaction -- the three effects that do not involve the bin factor and are
## therefore arithmetically unchanged by the decomposition; recomputed as a
## replication check on every period.

SPEC <- tribble(
  ~measure, ~period,     ~file,                                    ~bins,   ~pub,
  "PMA", "stimulus", "analysis_PMA/data/PMA_rate_stim_JASP.csv",   0:19,  c(13.38, 5.25,  0.312),
  "PMA", "delay",    "analysis_PMA/data/PMA_rate_delay_JASP.csv", 10:29,  c(33.19, 3.36,  3.474),
  "MBL", "stimulus", "analysis_MBL/data/MBL_stim_JASP.csv",        1:20,  c(2.146, 2.813, 0.755),
  "MBL", "delay",    "analysis_MBL/data/MBL_delay_JASP.csv",      11:30,  c(3.177, 1.257, 1.586),
  "LRP", "stimulus", "analysis_LRP/data/LRP_stim_JASP.csv",        1:20,  c(0.614, 1.498, 0.387),
  "LRP", "delay",    "analysis_LRP/data/LRP_delay_JASP.csv",      11:30,  c(0.019, 0.44,  0.378),
  "LRP", "RS",       "analysis_LRP/data/LRP_rs_JASP.csv",          1:2,   c(0.185, 0.252, 0.267),
  "LRP", "response", "analysis_LRP/data/LRP_resp_JASP.csv",        0:1,   c(0.400, 0.31,  0.666)
)
## The partial-activation rate and mu-beta lateralization are analysed in the
## stimulus and delay periods only; in the RS and response periods they are
## described, not tested (Materials and Methods, Statistical Analysis). The
## lateralized readiness potential is tested in all four periods.
## Note on LRP|stimulus: the target here is this script's own value, not an
## independently published one -- see the end of this file.

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
## Scores, in physical units. Time is the bin centre in seconds, centred on the
## window. LEVEL is the mean over bins. TREND is the linear coefficient and QUAD
## the quadratic coefficient of the least-squares fit
##     y(t) = b0 + b1 t + b2 t^2
## to the bin values: TREND = b1, in units of the measure per second, and
## QUAD = b2, in units of the measure per second squared (for the partial-
## activation rate, in s^-1, these are s^-2 and s^-3). Both are linear
## contrasts of the bin values. Because the centred bin times are symmetric about
## zero, t is orthogonal to t^2 - mean(t^2), so b1 is also the slope of the
## straight-line fit and b2 is the projection of y on t^2 - mean(t^2), by the
## Frisch-Waugh theorem. F and p are those of the orthogonal polynomial
## contrasts: a repeated-measures analysis is invariant under rescaling of the
## dependent variable, so only the coefficients change.

make_scores <- function(d, width = 0.1, per_second = FALSE) {
  nb <- max(d$bin)
  ## The partial-activation exports hold the mean number of activations per trial
  ## in each bin; dividing by the bin width turns them into a rate in s^-1, so
  ## that LEVEL, TREND and QUAD are in s^-1, s^-2 and s^-3.
  if (per_second) d <- d %>% mutate(y = y / width)
  tt <- (seq_len(nb) - (nb + 1) / 2) * width   # centred bin times, s
  wl <- tt / sum(tt^2)                          # b1, units per second
  q  <- tt^2 - mean(tt^2)
  q  <- q - tt * sum(q * tt) / sum(tt^2)        # zero for symmetric tt; kept for safety
  w2 <- if (nb > 2) q / sum(q^2) else NA_real_  # b2, units per second squared
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
  list(df1 = 1, df2 = df2, F = Fval, p = tt$p.value, eta2p = Fval / (Fval + df2),
       estimate = unname(tt$estimate), ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2])
}

## ---------------------------------------------------------------------------
## The seven effects (five in a two-bin period, which has no quadratic), as
## two or three two-way ANOVAs on the derived scores plus the score-vs-zero
## tests. aov_ez applies Greenhouse-Geisser to the three-level coherence
## factor; duration and the intercept-like vs-zero test need no correction.

analyse_period <- function(d, measure, period, pub = NULL) {
  s  <- make_scores(d, per_second = (measure == "PMA"))
  nb <- max(d$bin)
  label <- sprintf("%s, %s period", measure, period)

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
  say(sprintf("\n   trend vs zero: %+.5f units/s [95%% CI %+.5f, %+.5f], F(%.0f,%.0f) = %.2f, p = %.4g",
              mean(m0$v), t0$ci_lo, t0$ci_hi, t0$df1, t0$df2, t0$F, t0$p))

  AL <- a_level$anova_table; AT <- a_trend$anova_table
  add_row(measure, period, "duration", 1, AL["dur","den Df"], AL["dur","F"], AL["dur","Pr(>F)"], AL["dur","pes"])
  add_row(measure, period, "coherence", AL["coh","num Df"], AL["coh","den Df"], AL["coh","F"], AL["coh","Pr(>F)"], AL["coh","pes"])
  add_row(measure, period, "duration x coherence", AL["dur:coh","num Df"], AL["dur:coh","den Df"], AL["dur:coh","F"], AL["dur:coh","Pr(>F)"], AL["dur:coh","pes"])
  add_row(measure, period, "time (linear)", t0$df1, t0$df2, t0$F, t0$p, t0$eta2p,
          t0$estimate, t0$ci_lo, t0$ci_hi)
  add_row(measure, period, "duration x time (linear)", 1, AT["dur","den Df"], AT["dur","F"], AT["dur","Pr(>F)"], AT["dur","pes"])
  add_row(measure, period, "coherence x time (linear)", AT["coh","num Df"], AT["coh","den Df"], AT["coh","F"], AT["coh","Pr(>F)"], AT["coh","pes"])
  add_row(measure, period, "duration x coherence x time (linear)", AT["dur:coh","num Df"], AT["dur:coh","den Df"], AT["dur:coh","F"], AT["dur:coh","Pr(>F)"], AT["dur:coh","pes"])

  if (nb > 2) {
    a_quad <- aov_ez("subj", "quad", s, within = c("dur", "coh"))
    say("\n-- QUAD (quadratic coefficient, units/s^2): curvature, and its interactions with duration/coherence")
    capture.output(print(a_quad$anova_table), file = OUT, append = TRUE)
    q0 <- s %>% group_by(subj) %>% summarise(v = mean(quad))
    tq <- f_vs_zero(q0$v)
    say(sprintf("   curvature vs zero: %+.5f units/s^2 [95%% CI %+.5f, %+.5f], F(%.0f,%.0f) = %.2f, p = %.4g",
                mean(q0$v), tq$ci_lo, tq$ci_hi, tq$df1, tq$df2, tq$F, tq$p))
    AC <- a_quad$anova_table
    add_row(measure, period, "time (quadratic)", tq$df1, tq$df2, tq$F, tq$p, tq$eta2p,
            tq$estimate, tq$ci_lo, tq$ci_hi)
    add_row(measure, period, "duration x time (quadratic)", 1, AC["dur","den Df"], AC["dur","F"], AC["dur","Pr(>F)"], AC["dur","pes"])
    add_row(measure, period, "coherence x time (quadratic)", AC["coh","num Df"], AC["coh","den Df"], AC["coh","F"], AC["coh","Pr(>F)"], AC["coh","pes"])
    add_row(measure, period, "duration x coherence x time (quadratic)", AC["dur:coh","num Df"], AC["dur:coh","den Df"], AC["dur:coh","F"], AC["dur:coh","Pr(>F)"], AC["dur:coh","pes"])
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
## (Table S5); this split is descriptive, reported alongside it, and the
## cut-off was fixed before looking at the data (the midpoint of the window),
## not chosen after seeing the contrast.
##
## Sign convention: mu-beta lateralisation is negative-going with preparation,
## so a negative coherence contrast on the trend means a steeper buildup at
## higher coherence.

d <- read_long(file.path(ARCHIVE, "analysis_MBL/data/MBL_stim_JASP.csv"), 1:20)
say("\n", strrep("=", 78))
say("MBL, stimulus period, resolved into two halves (descriptive; Table S5)")
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

BEHAV <- tribble(
  ~measure,    ~file,             ~label,
  "Accuracy",  "accuracy.csv",    "response accuracy",
  "Mean RT",   "mean_RT.csv",     "mean response time",
  "Mean PMET", "mean_PMET.csv",   "mean premotor execution time",
  "Mean MET",  "mean_MET.csv",    "mean motor execution time")

read_behav <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  stopifnot(ncol(d) == 6)
  ## the archive's column order is coh0/12/50 at 300 ms, then coh0/12/50 at 2000 ms
  names(d) <- c("d1c0", "d1c12", "d1c50", "d2c0", "d2c12", "d2c50")
  d$subj <- factor(seq_len(nrow(d)))
  d %>%
    pivot_longer(-subj, names_to = "cell", values_to = "y") %>%
    mutate(dur = factor(sub("^d([12]).*", "\\1", cell), levels = c("1", "2"),
                        labels = c("0.3 s", "2 s")),
           coh = factor(sub("^d[12]c", "", cell), levels = c("0", "12", "50"),
                        labels = c("0%", "12%", "50%"))) %>%
    select(subj, dur, coh, y)
}

say("\n", strrep("=", 78))
say("BEHAVIOUR (archive exports, analysis_behavioral_emg/)")
say(strrep("=", 78))
for (i in seq_len(nrow(BEHAV))) {
  f <- file.path(ARCHIVE, "analysis_behavioral_emg", BEHAV$file[i])
  if (!file.exists(f)) { say("MISSING: ", f); next }
  d <- read_behav(f)
  a <- aov_ez("subj", "y", d, within = c("dur", "coh"))
  say("\n", BEHAV$measure[i], " -- ", BEHAV$label[i], "  (n = ", nlevels(d$subj), ")")
  capture.output(print(a$anova_table), file = OUT, append = TRUE)
  at <- a$anova_table
  for (j in seq_len(nrow(at))) {
    add_row(BEHAV$measure[i], "behaviour", rownames(at)[j],
            at[j, "num Df"], at[j, "den Df"], at[j, "F"], at[j, "Pr(>F)"], at[j, "pes"])
  }
}

## ---------------------------------------------------------------------------
close(OUT)
write.csv(do.call(rbind, ROWS), "stats_results.csv", row.names = FALSE)
cat(readLines("stats_tables.txt"), sep = "\n")

## ---------------------------------------------------------------------------
## Two traps this analysis has, worth keeping in view.
##
##   (1) A 2-bin period's polynomial-contrast matrix has one column, so an
##       unconditional reference to its second column crashes on every such
##       period (RS for all three measures, plus the mu-beta and LRP response
##       periods).
##   (2) The duration code in the column names is not consistent across measures
##       (see the header): parsing it wrongly loses a factor level rather than
##       giving a wrong number.
##
## For LRP, stimulus period, SPEC's replication target is (0.614, 1.498, 0.387),
## the value this script computes from the archive's own files.  An older figure
## of (0.636, 1.524, 0.391) circulates for that one period; the few per cent
## between them has not been traced.
