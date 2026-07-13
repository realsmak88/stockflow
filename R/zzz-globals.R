###############################################################
#
# stockflow
#
# Module : zzz-globals.R
#
# Declaration of variables used in non-standard evaluation
# (ggplot2::aes, formulas) in order to avoid R CMD check NOTES
# of the type "no visible binding for global variable".
#
###############################################################

utils::globalVariables(c(
  # ggplot2 (aes) - lbb.R, lbspr.R, mse.R
  "Year", "val", "ind", "SPR", "Yield_rel", "scenario",
  "BB0", "BBmsy", "FM", "SL50", "SL95",
  # mortality.R / allometry.R
  "M", "method", "L", "W", "logL", "logW", "grp", "mat",
  # growth.R
  "age", "length", "Linf", "K", "phiL", "Rn_max",
  # divers
  "Length", "Freq", "CatchNo", "Lmean",
  # plots.R (aes en evaluation non standard)
  "month_label", "score", "protege", "lo", "hi",
  ".periode", ".gsi", ".mois", ".bin", ".mature",
  "mid", "prop", "prop_f", "n", "used", "composante", "valeur",
  "nom", ".L", ".W", ".lo", ".hi", ".cond",
  "quota_t", "segment", "prorata", "Lc", "etape",
  "lnC", "annee", "groupe", "mois", "moy",
  "mois_lab", "gsi", "sex", "jeu",
  "verdict", "val", ".data"
))
