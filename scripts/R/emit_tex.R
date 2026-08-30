## emit_tex.R -- turn stats_results.csv into the LaTeX tables of the note.
## Run after run_stats.R.  Published values are inlined so that each table shows
## the manuscript's number beside the new one.
r <- read.csv("stats_results.csv")
PUB <- list(
 "PMA|stimulus"=c(13.38,5.25,18.05,0.312,4.645,1.718,1.058),
 "PMA|delay"   =c(33.19,3.36,31.42,3.474,4.520,0.737,0.959),
 "PMA|RS"      =c(36.56,1.159,21.96,0.041,7.046,0.407,0.587),
 "PMA|response"=c(0.011,0.724,38.31,0.636,5.013,1.172,0.193),
 "MBL|stimulus"=c(2.146,2.813,6.739,0.755,1.823,2.805,0.810),
 "MBL|delay"   =c(3.177,1.257,5.642,1.586,1.078,1.008,1.199),
 "MBL|RS"      =c(3.039,0.160,0.059,1.216,1.509,0.024,1.458),
 "MBL|response"=c(0.006,0.364,11.66,0.161,2.129,0.252,0.719),
 "LRP|stimulus"=c(0.636,1.524,1.876,0.391,1.181,1.177,0.595),
 "LRP|delay"   =c(0.019,0.44,10.791,0.378,0.555,0.457,0.408),
 "LRP|RS"      =c(0.185,0.252,14.076,0.267,1.238,1.058,0.108),
 "LRP|response"=c(0.400,0.31,54.331,0.666,4.129,0.246,0.589))
TEX <- c("duration"="duration","coherence"="coherence",
 "time (linear)"="time, linear","duration x coherence"="duration $\\times$ coherence",
 "duration x time (linear)"="duration $\\times$ time, linear",
 "coherence x time (linear)"="coherence $\\times$ time, linear",
 "duration x coherence x time (linear)"="duration $\\times$ coherence $\\times$ time, linear",
 "time (quadratic)"="time, quadratic",
 "duration x time (quadratic)"="duration $\\times$ time, quadratic",
 "coherence x time (quadratic)"="coherence $\\times$ time, quadratic",
 "duration x coherence x time (quadratic)"="duration $\\times$ coherence $\\times$ time, quadratic")
fp <- function(p) if (p < .001) "$<$.001" else sub("^0","",sprintf("%.3f",p))
NAME <- c(PMA="Partial-activation rate", MBL="Mu--beta lateralisation",
          LRP="Lateralised readiness potential")
con <- file("_tables.tex","w")
for (m in c("PMA","MBL","LRP")) for (per in c("stimulus","delay","RS","response")) {
  z <- r[r$measure==m & r$period==per,]; if (!nrow(z)) next
  pub <- PUB[[paste(m,per,sep="|")]]
  cat(sprintf("\\paragraph{%s, %s period.}\n\\begin{center}\\small\n\\begin{tabular}{lcccccc}\n\\toprule\n & \\multicolumn{1}{c}{published} & \\multicolumn{5}{c}{trend decomposition}\\\\\n\\cmidrule(lr){2-2}\\cmidrule(l){3-7}\neffect & $F$ & df1 & df2 & $F$ & $p$ & $\\eta_p^2$\\\\\n\\midrule\n",
      NAME[m], per), file=con)
  ## published effects are indexed in the manuscript's order; when a period is
  ## analysed on one bin only, three effects are reported and the published
  ## time-related rows have no counterpart
  pubmap <- c("duration"=1, "coherence"=2, "time (linear)"=3,
              "duration x coherence"=4, "duration x time (linear)"=5,
              "coherence x time (linear)"=6,
              "duration x coherence x time (linear)"=7)
  for (i in seq_len(nrow(z))) {
    k <- as.character(z$effect[i])
    j <- pubmap[k]
    pf <- if (!is.na(j)) sprintf("%.3f", pub[j]) else "---"
    cat(sprintf("%s & %s & %.2f & %.2f & %.3f & %s & %.3f\\\\\n",
        TEX[k], pf, z$df1[i], z$df2[i], z$F[i], fp(z$p[i]), z$eta2p[i]), file=con)
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{center}\n\n", file=con)
}
close(con); cat("wrote _tables.tex\n")
