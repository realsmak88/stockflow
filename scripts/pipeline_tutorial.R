###############################################################
# stockflow — Complete pipeline tutorial (.R script version)
#
# Run from the ROOT of the package, line by line in RStudio.
# Narrative version (HTML): vignettes/tutoriel_stockflow.Rmd
#
# Pipeline: import -> validation -> LFQ -> growth -> mortality
#            -> allometry -> LBB -> LBSPR -> management measures
###############################################################

## 0. Setup -------------------------------------------------------
# install.packages("devtools")
devtools::load_all(".")

demo_lengths <- system.file("extdata", "penaeus_lengths_demo.csv",
                            package = "stockflow")
demo_bio     <- system.file("extdata", "penaeus_bio_demo.csv",
                            package = "stockflow")
if (demo_lengths == "") demo_lengths <- "inst/extdata/penaeus_lengths_demo.csv"
if (demo_bio == "")     demo_bio     <- "inst/extdata/penaeus_bio_demo.csv"

# Scientific packages (install once):
# install.packages(c("TropFishR", "FSA", "LBSPR", "fishmethods"))


## 1. Import ---------------------------------------------------------
lengths <- read_lengths(demo_lengths)
head(lengths)


## 2. Quality validation --------------------------------------------------
val <- validate_lengths(lengths, species_config = list(),
                        min_length = 0, max_length = 60)
print(val)


## 3. Preparing the LFQ object (TropFishR) ------------------------------
dd <- lengths
names(dd)[names(dd) == "length"] <- "Length"
names(dd)[names(dd) == "year"]   <- "Year"
names(dd)[names(dd) == "month"]  <- "Month"

lfq <- prepare_tropfish(dd, bin_size = 1,
                        species = "Penaeus notialis",
                        stock_name = "SEN-artisanal",
                        length_unit = "cm")

lfq <- TropFishR::lfqRestructure(lfq, MA = 5, addl.sqrt = FALSE)
graphics::plot(lfq, Fname = "rcounts", date.axis = "modern")


## 4. Growth ----------------------------------------------------------
# 4a. Powell-Wetherall
pw <- run_powell(lfq)
summary(pw)

# 4b. ELEFAN GA + SA + automatic selection (1-3 min)
set.seed(42)
growth <- run_growth_analysis(
  lfq,
  methods    = c("elefan_ga", "elefan_sa"),
  Linf_range = c(40, 55), K_range = c(0.4, 2.0),
  criterion  = "Rn_max",
  control = list(
    low_par = list(Linf = 40, K = 0.4, t_anchor = 0),
    up_par  = list(Linf = 55, K = 2.0, t_anchor = 1),
    popSize = 40, maxiter = 30, seasonalised = FALSE
  ),
  verbose = FALSE
)
growth$comparison
best <- growth$best
Linf <- best$Linf ; K <- best$K

# 4c. Growth curve
crb <- vbgf_curve(Linf = Linf, K = K, t_anchor = 0.5, ages = seq(0, 6, 0.1))
plot(crb$age, crb$length, type = "l", lwd = 2, col = "#12507b",
     xlab = "Age", ylab = "Length")


## 5. Mortality -----------------------------------------------------------
# 5a. M — all methods (TropFishR + FSA::metaM + fishmethods + Jensen)
M_all <- estimate_M_all(Linf = Linf, K = K, temp = 27)
M_all$table
M_all$consensus

# 5b. Z — catch curve (choose reg_int from the plot)
# z <- estimate_catchcurve(lfq, Linf = Linf, K = K,
#                          reg_int = c(15, 30), calc_ogive = TRUE, plot = TRUE)

# 5c. Full Z + M + F/E workflow
# mort <- run_mortality(lfq, growth_model = growth$models[["elefan_ga"]],
#                       temp = 27, reg_int = c(15, 30))
# summary(mort)


## 6. Length-weight allometry ---------------------------------------------
bio <- utils::read.csv(demo_bio)
lw  <- fit_length_weight(bio, "length", "weight")
print(lw)


## 7. Maturity (L50 / L95) for LBSPR -------------------------------------
bio$mature <- as.integer(bio$maturity >= 3)
gm  <- stats::glm(mature ~ length, data = bio, family = binomial())
cf  <- stats::coef(gm)
L50 <- as.numeric(-cf[1] / cf[2])
L95 <- as.numeric((log(0.95/0.05) - cf[1]) / cf[2])
c(L50 = L50, L95 = L95)


## 8. LBB -----------------------------------------------------------------
dlbb <- prepare_lbb_data(lengths, "length", "year", stock = "Penaeus_notialis")

# Froese's official backend script (LBB_ggplot.R, ggplot2) — requires R2jags + JAGS
# fit_lbb <- run_lbb_froese(dlbb, stock = "Penaeus_notialis",
#                           mm = TRUE, MK_prior = 1.5, Linf_prior = 45, Lm50 = L50)
# fit_lbb$summary ; fit_lbb$output_dir ; plot(fit_lbb)

lbb_reference_points(Linf = Linf, MK = 1.5, FM = 1)


## 9. LBSPR ---------------------------------------------------------------
pr <- lbspr_pars(Linf = Linf, L50 = L50, L95 = L95,
                 M = M_all$consensus[["geomean"]], K = K,
                 species = "Penaeus notialis")
ln <- lbspr_lengths(lengths, pr, length_col = "length",
                    year_col = "year", bin_width = 1)
fit_sp <- run_lbspr(pr, ln, spr_target = 0.40, spr_limit = 0.20)
fit_sp$summary
plot(fit_sp)


## 10. Management measures (equilibrium MSE) -------------------------------
scenarios <- list(
  mse_measure("Status quo",     "statuquo"),
  mse_measure("-30% effort",   "effort",           effort_reduction = 0.30),
  mse_measure("3-month closure",  "repos_biologique", closure_months  = 3),
  mse_measure("MPA 20%",       "amp",              mpa_fraction     = 0.20),
  mse_measure("Min size 25", "taille_min",       Lc = 25),
  mse_optimal_size(Linf = Linf, MK = 1.5)
)
mse <- run_mse_equilibrium(pr, scenarios, base_FM = 1.5)
mse$comparison
plot(mse)

## 11. Automatic report (HTML / PDF / DOCX) ---------------------------
resultats <- collect_results(
  growth    = growth,
  allometry = lw,
  maturity  = list(L50 = L50, L95 = L95),
  lbspr     = fit_sp,
  mse       = mse,
  meta = list(species = "Penaeus notialis",
              stock   = "Example EEZ",
              period  = "2025-2026",
              n_ind   = nrow(lengths))
)
resultats

fishstock_report(resultats, file = "rapport_penaeus", format = "html",
                 title  = "White shrimp stock assessment",
                 author = "Kamarel Ba")

###############################################################
# End of tutorial — happy modeling!
###############################################################

