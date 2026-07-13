# =====================================================================
# [ARCHIVE — inst/legacy/] Old procedural script.
# DO NOT PLACE IN R/ ANYMORE: this file executes code on load
# (library(), stop(), options(), set.seed()) and interrupts
# devtools::load_all(). The active functions live in R/growth.R.
# Kept only for historical reference purposes.
# =====================================================================
#
# growth.R
# Length growth analysis with TropFishR
#
# Objective:
# Estimate von Bertalanffy growth parameters
# (Linf, K) from length frequency (LFQ) data
#
# Functions used exclusively from TropFishR:
# - lfqCreate()
# - lfqModify()
# - lfqRestructure()
# - powell_wetherall()
# - ELEFAN()
# - ELEFAN_GA()
# - ELEFAN_SA()
# - catchCurve()
# - M_empirical()
#
# No internal algorithm is reimplemented.
# =====================================================================


# =====================================================================
# 0. Cleanup and package loading
# =====================================================================

if (!requireNamespace("TropFishR", quietly = TRUE))
  stop("The 'TropFishR' package must be installed.")

library(TropFishR)


# ---------------------------------------------------------------------
# General options
# ---------------------------------------------------------------------

options(stringsAsFactors = FALSE)

set.seed(123)


# =====================================================================
# 1. File organization
# =====================================================================

# The input file must contain at minimum:
#
# - individual length or length class
# - sampling date
#
# Example:
#
# length,date
# 12.5,2020-01-15
# 13.0,2020-01-15
# 18.5,2020-02-15
#
# ---------------------------------------------------------------------

input_file <- "length_frequency_data.csv"


# =====================================================================
# 2. Data import
# =====================================================================

if (file.exists(input_file)) {
  
  raw_data <- read.csv(
    input_file,
    header = TRUE
  )
  
} else {
  
  message(
    "File missing: using the TropFishR synLFQ4 example dataset"
  )
  
  data(TropFishR::synLFQ4)
  
}


# =====================================================================
# 3. Creating the LFQ object
# =====================================================================

# Two possibilities:
#
# A) Personal data:
#    created with lfqCreate()
#
# B) TropFishR dataset:
#    direct use of an existing LFQ object
#
# ---------------------------------------------------------------------


if (exists("raw_data")) {
  
  
  lfq <- TropFishR::lfqCreate(
    
    data = raw_data,
    
    Lname = "length",
    
    Dname = "date",
    
    Fname = NULL,
    
    bin_size = 2,
    
    species = "species_name",
    
    stock = "stock_name",
    
    comment =
      "LFQ object created with TropFishR"
    
  )
  
  
} else {
  
  
  lfq <- TropFishR::synLFQ4
  
  
}



# =====================================================================
# 4. Inspecting the LFQ object
# =====================================================================

print(lfq)


str(lfq)


summary(lfq)



# =====================================================================
# 5. Initial graphical check
# =====================================================================

# Visualization of length frequencies
# before any transformation

plot(
  lfq,
  Fname = "catch",
  date.axis = "modern"
)



# =====================================================================
# 6. Modifying the LFQ object
# =====================================================================

# lfqModify() allows adapting the LFQ structure:
#
# - changing length classes
# - changing parameters
# - preparation for ELEFAN
#
# No statistical processing is performed here.
# ---------------------------------------------------------------------


lfq_mod <- TropFishR::lfqModify(
  
  lfq,
  
  bin_size = 2
  
)



# Check

print(lfq_mod)



# =====================================================================
# 7. Visualization after modification
# =====================================================================


plot(
  
  lfq_mod,
  
  Fname = "catch",
  
  date.axis = "modern"
  
)



# =====================================================================
# 8. Preparation for ELEFAN
# =====================================================================

# The ELEFAN method requires a restructuring
# of the length frequencies.
#
# This step is performed only by:
#
#     lfqRestructure()
#
# from TropFishR.
#
# No score calculation is performed manually.
# ---------------------------------------------------------------------


MA <- 7


lfq_res <- TropFishR::lfqRestructure(
  
  lfq_mod,
  
  MA = MA,
  
  addl.sqrt = FALSE
  
)



# Check

print(lfq_res)



# =====================================================================
# 9. Restructuring plot
# =====================================================================


plot(
  
  lfq_res,
  
  Fname = "rcounts",
  
  date.axis = "modern"
  
)



# =====================================================================
# 10. Intermediate save
# =====================================================================


save(
  
  lfq,
  
  lfq_mod,
  
  lfq_res,
  
  file = "LFQ_preparation.RData"
  
)



# =====================================================================
# End Part 1/4
# =====================================================================

# =====================================================================
# growth.R
# PART 2/4
#
# Estimation of growth parameters
# using the official TropFishR methods
#
# Functions used:
# - powell_wetherall()
# - ELEFAN()
# - ELEFAN_GA()
# - ELEFAN_SA()
#
# =====================================================================



# =====================================================================
# 11. Powell-Wetherall
# =====================================================================

# The Powell-Wetherall method provides a preliminary
# estimate of Linf from the LFQ data.
#
# It is often used as a starting point before ELEFAN.
#
# ---------------------------------------------------------------------

# Powell-Wetherall visualization
pw <- TropFishR::powell_wetherall(
  
  lfq_mod,
  
  catch_columns = 1:ncol(lfq_mod$catch),
  
  reg_int = NULL
  
)


print(pw)

# Optional extraction of the estimated Linf

Linf_pw <- pw$Linf_est


print(Linf_pw)

# =====================================================================
# 12. Defining the search intervals
# =====================================================================

# The bounds are adapted to the stock under study.
#
# They do not replace the TropFishR algorithm.
# They only indicate the search space.
# ---------------------------------------------------------------------

Linf_min <- max(lfq_mod$midLengths)


Linf_max <- max(lfq_mod$midLengths) * 3



K_min <- 0.05


K_max <- 2



lower_parameters <- list(
  
  Linf = Linf_min,
  
  K = K_min
  
)



upper_parameters <- list(
  
  Linf = Linf_max,
  
  K = K_max
  
)



# =====================================================================
# 13. Classic ELEFAN
# =====================================================================

# Original ELEFAN method available in TropFishR.
#
# It searches for the best Linf/K combination
# according to the ELEFAN score.
#
# ---------------------------------------------------------------------


set.seed(123)



elefan_result <- TropFishR::ELEFAN(
  
  lfq_res,
  
  Linf_range =
    seq(
      
      from = Linf_min,
      
      to = Linf_max,
      
      length.out = 100
      
    ),
  
  K_range =
    seq(
      
      from = K_min,
      
      to = K_max,
      
      length.out = 100
      
    ),
  
  MA = MA,
  
  plot = TRUE
  
)



print(elefan_result)



# Estimated parameters

elefan_par <- elefan_result$par


print(elefan_par)



# Score obtained

elefan_score <- elefan_result$score


print(elefan_score)



# =====================================================================
# 14. ELEFAN Genetic Algorithm
# =====================================================================

# ELEFAN_GA uses the genetic optimization built
# into TropFishR.
#
# No genetic algorithm is coded here.
#
# ---------------------------------------------------------------------


set.seed(123)



ga_result <- TropFishR::ELEFAN_GA(
  
  lfq_res,
  
  seasonalised = FALSE,
  
  low_par = lower_parameters,
  
  up_par = upper_parameters,
  
  popSize = 50,
  
  maxiter = 100,
  
  run = 20,
  
  plot = TRUE,
  
  plot.score = TRUE
  
)



print(ga_result)



# GA parameters

ga_par <- ga_result$par


print(ga_par)



# GA score

ga_score <- ga_result$score


print(ga_score)



# =====================================================================
# 15. ELEFAN Simulated Annealing
# =====================================================================

# Simulated annealing optimization built into TropFishR.
#
# ---------------------------------------------------------------------


set.seed(123)



sa_result <- TropFishR::ELEFAN_SA(
  
  lfq_res,
  
  seasonalised = FALSE,
  
  low_par = lower_parameters,
  
  up_par = upper_parameters,
  
  SA_time = 60,
  
  MA = MA,
  
  plot = TRUE
  
)



print(sa_result)



# SA parameters

sa_par <- sa_result$par


print(sa_par)



# SA score

sa_score <- sa_result$score


print(sa_score)



# =====================================================================
# 16. Comparison of results
# =====================================================================


growth_comparison <- data.frame(
  
  Methode = c(
    
    "ELEFAN",
    
    "ELEFAN_GA",
    
    "ELEFAN_SA"
    
  ),
  
  Linf = c(
    
    elefan_par$Linf,
    
    ga_par$Linf,
    
    sa_par$Linf
    
  ),
  
  K = c(
    
    elefan_par$K,
    
    ga_par$K,
    
    sa_par$K
    
  ),
  
  Score = c(
    
    elefan_score,
    
    ga_score,
    
    sa_score
    
  )
  
)



print(growth_comparison)



# =====================================================================
# 17. Automatic selection of the best model
# =====================================================================

# The best ELEFAN score is retained.
#
# The choice is based solely on the result
# provided by TropFishR.
#
# ---------------------------------------------------------------------


best_index <- which.max(
  
  growth_comparison$Score
  
)



best_method <- growth_comparison$Methode[best_index]


cat(
  
  "Method selected:",
  
  best_method,
  
  "\n"
  
)



best_parameters <- list(
  
  Linf =
    growth_comparison$Linf[best_index],
  
  K =
    growth_comparison$K[best_index]
  
)



print(best_parameters)



# =====================================================================
# 18. Saving the estimates
# =====================================================================


save(
  
  pw,
  
  elefan_result,
  
  ga_result,
  
  sa_result,
  
  growth_comparison,
  
  best_parameters,
  
  file =
    "growth_estimations.RData"
  
)



# =====================================================================
# End Part 2/4
# =====================================================================

# =====================================================================
# growth.R
# PART 3/4
#
# Application of growth parameters
# Mortality analysis
#
# TropFishR functions used:
#
# - lfqModify()
# - catchCurve()
# - M_empirical()
#
# No internal statistical calculation is reimplemented.
#
# =====================================================================



# =====================================================================
# 19. Application of the selected growth parameters
# =====================================================================


# The Linf and K parameters come directly
# from the best TropFishR solution obtained previously.


Linf_final <- best_parameters$Linf


K_final <- best_parameters$K



growth_parameters <- list(
  
  Linf = Linf_final,
  
  K = K_final
  
)



print(growth_parameters)



# =====================================================================
# 20. Adding the growth curve to the LFQ
# =====================================================================

# lfqModify() allows associating the growth
# parameters with the LFQ object.
#
# The official TropFishR function is used.
# ---------------------------------------------------------------------


lfq_growth <- TropFishR::lfqModify(
  
  lfq_mod,
  
  par = growth_parameters
  
)



print(lfq_growth)



# =====================================================================
# 21. Visualization of the growth curve
# =====================================================================


plot(
  
  lfq_growth,
  
  Fname = "catch",
  
  date.axis = "modern"
  
)



# =====================================================================
# 22. Preparation for catchCurve
# =====================================================================

# The catch curve requires the length classes
# and the corresponding catches.
#
# The data remain in the TropFishR LFQ structure.
#
# ---------------------------------------------------------------------


lfq_cc <- lfq_growth



# Check the required elements

names(lfq_cc)



# =====================================================================
# 23. Catch curve analysis
# =====================================================================

# catchCurve() estimates total mortality Z
# from the descending part of the curve.
#
# The point selection remains the one performed
# by the TropFishR function.
#
# ---------------------------------------------------------------------


catch_curve <- TropFishR::catchCurve(
  
  param = lfq_cc,
  
  Linf = Linf_final,
  
  K = K_final,
  
  plot = TRUE
  
)



print(catch_curve)



# =====================================================================
# 24. Extraction of mortality parameters
# =====================================================================


Z_estimate <- catch_curve$Z


print(Z_estimate)



# =====================================================================
# 25. Estimation of natural mortality M
# =====================================================================

# M_empirical() applies the empirical equations
# available in TropFishR.
#
# Several methods can be tested depending on
# the ecological characteristics of the stock.
#
# ---------------------------------------------------------------------


M_pauly <- TropFishR::M_empirical(
  
  Linf = Linf_final,
  
  K = K_final,
  
  method = "Pauly"
  
)



print(M_pauly)



# ---------------------------------------------------------------------
# Example with an average water temperature
# (to be adapted to the stock under study)
# ---------------------------------------------------------------------


M_pauly_temp <- TropFishR::M_empirical(
  
  Linf = Linf_final,
  
  K = K_final,
  
  method = "Pauly",
  
  temp = 25
  
)



print(M_pauly_temp)



# =====================================================================
# 26. Calculation of the exploitation rate
# =====================================================================

# Fishing mortality is derived from:
#
# F = Z - M
#
# This operation uses only
# the TropFishR outputs obtained.
#


F_estimate <- Z_estimate - M_pauly_temp



E_estimate <- F_estimate / Z_estimate



mortality_results <- data.frame(
  
  Parametre = c(
    
    "Linf",
    
    "K",
    
    "Z",
    
    "M",
    
    "F",
    
    "E"
    
  ),
  
  Value = c(
    
    Linf_final,
    
    K_final,
    
    Z_estimate,
    
    M_pauly_temp,
    
    F_estimate,
    
    E_estimate
    
  )
  
)



print(mortality_results)



# =====================================================================
# 27. Additional graphical diagnostics
# =====================================================================


# Growth curve alone

plot(
  
  lfq_growth,
  
  Fname = "catch",
  
  date.axis = "modern",
  
  draw = TRUE
  
)



# =====================================================================
# 28. Saving mortality results
# =====================================================================


save(
  
  lfq_growth,
  
  catch_curve,
  
  Z_estimate,
  
  M_pauly,
  
  M_pauly_temp,
  
  mortality_results,
  
  file =
    
    "mortality_results.RData"
  
)



# =====================================================================
# End Part 3/4
# =====================================================================

# =====================================================================
# growth.R
# PART 4/4
#
# Final diagnostics
# Export of results
# Complete save
#
# =====================================================================



# =====================================================================
# 29. Biological check of the growth parameters
# =====================================================================

# The values obtained must be interpreted
# according to the species under study.
#
# This section does not modify the TropFishR results.
# It only provides control alerts.
# ---------------------------------------------------------------------


if (Linf_final <= max(lfq_mod$midLengths)) {
  
  warning(
    
    "Warning: estimated Linf is lower than or close to the maximum observed length."
    
  )
  
}



if (K_final <= 0) {
  
  warning(
    
    "Warning: negative or zero K coefficient."
    
  )
  
}



cat("\n")

cat("------------------------------------\n")

cat("Selected growth parameters\n")

cat("------------------------------------\n")

cat("Linf =", Linf_final, "\n")

cat("K    =", K_final, "\n")

cat("------------------------------------\n")



# =====================================================================
# 30. Checking the ELEFAN score
# =====================================================================


cat("\n")

cat("------------------------------------\n")

cat("Comparison of ELEFAN scores\n")

cat("------------------------------------\n")


print(
  
  growth_comparison
  
)



# =====================================================================
# 31. Export of growth parameters
# =====================================================================


write.csv(
  
  growth_comparison,
  
  file = "ELEFAN_growth_comparison.csv",
  
  row.names = FALSE
  
)



write.csv(
  
  as.data.frame(best_parameters),
  
  file = "best_growth_parameters.csv",
  
  row.names = FALSE
  
)



# =====================================================================
# 32. Export of mortality results
# =====================================================================


write.csv(
  
  mortality_results,
  
  file = "mortality_parameters.csv",
  
  row.names = FALSE
  
)



# =====================================================================
# 33. Creation of an automatic summary
# =====================================================================


summary_text <- paste(
  
  "TropFishR Analysis - Summary",
  
  "",
  
  paste(
    
    "Method selected:",
    
    best_method
    
  ),
  
  paste(
    
    "Linf:",
    
    round(Linf_final, 3)
    
  ),
  
  paste(
    
    "K:",
    
    round(K_final, 3)
    
  ),
  
  paste(
    
    "Z:",
    
    round(Z_estimate, 3)
    
  ),
  
  paste(
    
    "M:",
    
    round(M_pauly_temp, 3)
    
  ),
  
  paste(
    
    "F:",
    
    round(F_estimate, 3)
    
  ),
  
  paste(
    
    "E:",
    
    round(E_estimate, 3)
    
  ),
  
  sep = "\n"
  
)



cat(summary_text)



writeLines(
  
  summary_text,
  
  con = "TropFishR_growth_summary.txt"
  
)



# =====================================================================
# 34. Complete save of the analysis
# =====================================================================


save(
  
  lfq,
  
  lfq_mod,
  
  lfq_res,
  
  pw,
  
  elefan_result,
  
  ga_result,
  
  sa_result,
  
  growth_comparison,
  
  best_parameters,
  
  lfq_growth,
  
  catch_curve,
  
  mortality_results,
  
  file = "TropFishR_complete_growth_analysis.RData"
  
)



# =====================================================================
# 35. Session information
# =====================================================================


sink(
  
  "TropFishR_sessionInfo.txt"
  
)

sessionInfo()

sink()



# =====================================================================
# 36. Final message
# =====================================================================


cat("\n")

cat("====================================================\n")

cat("TropFishR analysis completed successfully\n")

cat("====================================================\n")

cat("\n")

cat("Generated files:\n")

cat("- ELEFAN_growth_comparison.csv\n")

cat("- best_growth_parameters.csv\n")

cat("- mortality_parameters.csv\n")

cat("- TropFishR_growth_summary.txt\n")

cat("- TropFishR_complete_growth_analysis.RData\n")

cat("\n")



# =====================================================================
# END OF SCRIPT growth.R
# =====================================================================

# =====================================================================
# growth.R
# PART 5/5
#
# Multi-method comparison of natural mortality M
#
# Packages used:
#
# - TropFishR
# - FSA
# - fishmethods
#
# Objective:
# compare several empirical estimates of M
# and define a robust synthetic value.
#
# No mortality model is reimplemented.
#
# =====================================================================



# =====================================================================
# 37. Loading complementary packages
# =====================================================================


if (!require("FSA")) {
  
  install.packages("FSA")
  
}


if (!require("fishmethods")) {
  
  install.packages("fishmethods")
  
}



library(FSA)

library(fishmethods)



# =====================================================================
# 38. Biological parameters used
# =====================================================================


# The parameters come directly
# from the previous ELEFAN analysis.


Linf_M <- Linf_final

K_M <- K_final



temperature_M <- 25



# =====================================================================
# 39. Estimation of M with TropFishR
# =====================================================================


M_TropFishR <- TropFishR::M_empirical(
  
  Linf = Linf_M,
  
  K = K_M,
  
  method = "Pauly",
  
  temp = temperature_M
  
)



print(M_TropFishR)



# =====================================================================
# 40. Estimation of M with FSA
# =====================================================================


# Empirical methods available in FSA
#
# Check the arguments according to the version:
# ?Mmeta
# ?Mmethods
#
# ---------------------------------------------------------------------


M_FSA_meta <- tryCatch(
  
  {
    
    FSA::Mmeta(
      
      Linf = Linf_M,
      
      K = K_M
      
    )
    
  },
  
  error = function(e) NA
  
)



M_FSA_methods <- tryCatch(
  
  {
    
    FSA::Mmethods(
      
      Linf = Linf_M,
      
      K = K_M
      
    )
    
  },
  
  error = function(e) NA
  
)



print(M_FSA_meta)

print(M_FSA_methods)



# =====================================================================
# 41. Estimation of M with fishmethods
# =====================================================================


M_fishmethods <- tryCatch(
  
  {
    
    fishmethods::M.empirical(
      
      Linf = Linf_M,
      
      K = K_M
      
    )
    
  },
  
  error = function(e) NA
  
)



print(M_fishmethods)



# =====================================================================
# 42. Building the comparison table
# =====================================================================


M_results <- data.frame(
  
  Methode = c(
    
    "TropFishR_Pauly",
    
    "FSA_Mmeta",
    
    "FSA_Mmethods",
    
    "fishmethods"
    
  ),
  
  M = c(
    
    as.numeric(M_TropFishR),
    
    as.numeric(M_FSA_meta),
    
    as.numeric(M_FSA_methods),
    
    as.numeric(M_fishmethods)
    
  )
  
)



# Removal of impossible results

M_results <- M_results[
  
  is.finite(M_results$M),
  
]



print(M_results)



# =====================================================================
# 43. Graphical analysis of the estimates
# =====================================================================


barplot(
  
  M_results$M,
  
  names.arg = M_results$Methode,
  
  las = 2,
  
  ylab = "Natural mortality M",
  
  main = "Comparison of M estimates"
  
)



# =====================================================================
# 44. Descriptive statistics
# =====================================================================


M_summary <- data.frame(
  
  Minimum =
    min(M_results$M),
  
  Median =
    median(M_results$M),
  
  Moyenne =
    mean(M_results$M),
  
  Maximum =
    max(M_results$M)
  
)



print(M_summary)



# =====================================================================
# 45. Choice of final M
# =====================================================================


# The median is used as a robust estimator
# when several methods are available.
#
# It limits the influence of extreme values.
# ---------------------------------------------------------------------


M_final <- median(
  
  M_results$M
  
)



cat(
  
  "Final M selected =",
  
  M_final,
  
  "\n"
  
)



# =====================================================================
# 46. New estimation of F and E
# =====================================================================


F_final <- Z_estimate - M_final



E_final <- F_final / Z_estimate



exploitation_final <- data.frame(
  
  Parametre = c(
    
    "Z",
    
    "M_final",
    
    "F",
    
    "E"
    
  ),
  
  Value = c(
    
    Z_estimate,
    
    M_final,
    
    F_final,
    
    E_final
    
  )
  
)



print(exploitation_final)



# =====================================================================
# 47. Sensitivity analysis around M
# =====================================================================


M_sensitivity <- data.frame(
  
  M = seq(
    
    min(M_results$M),
    
    max(M_results$M),
    
    length.out = 20
    
  )
  
)



M_sensitivity$F <-
  
  Z_estimate -
  
  M_sensitivity$M



M_sensitivity$E <-
  
  M_sensitivity$F /
  
  Z_estimate



plot(
  
  M_sensitivity$M,
  
  M_sensitivity$E,
  
  type = "b",
  
  xlab = "M",
  
  ylab = "Exploitation rate E",
  
  main = "Sensitivity of E to M"
  
)



# =====================================================================
# 48. Export
# =====================================================================


write.csv(
  
  M_results,
  
  "M_comparison_methods.csv",
  
  row.names = FALSE
  
)



write.csv(
  
  exploitation_final,
  
  "exploitation_final.csv",
  
  row.names = FALSE
  
)



write.csv(
  
  M_sensitivity,
  
  "M_sensitivity.csv",
  
  row.names = FALSE
  
)



# =====================================================================
# 49. Final save
# =====================================================================


save(
  
  M_results,
  
  M_summary,
  
  M_final,
  
  exploitation_final,
  
  M_sensitivity,
  
  file =
    
    "M_multi_methods_results.RData"
  
)



# =====================================================================
# END PART 5/5
# END growth.R
# =====================================================================