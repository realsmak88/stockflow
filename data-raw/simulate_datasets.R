## =====================================================================
## simulate_datasets.R
## ---------------------------------------------------------------------
## Generates the FICTITIOUS datasets shipped with the stockflow package.
##
## IMPORTANT: these data are entirely SIMULATED. They reproduce the
## STRUCTURE (columns, types, factor levels, plausible ranges) of
## fisheries monitoring data, but contain NO real observations. Landing
## sites and fishing zones are fictitious. They exist only to illustrate
## and test the package functions.
##
## Reproducible: fixed set.seed(). To regenerate:
##   source("data-raw/simulate_datasets.R")
## =====================================================================

set.seed(1234)

## ---- FICTITIOUS factor levels (no real place) ----
SITES   <- c("Landing site A", "Landing site B", "Landing site C")
ZONES   <- paste("Fishing zone", sprintf("%02d", 1:10))
SEASONS <- c("Warm", "Cold", "Warm-to-cold", "Cold-to-warm")

## ---- helpers ----
# draw dates over an annual cycle
.rdate <- function(n, from = "2025-01-01", to = "2025-12-31") {
  d0 <- as.Date(from); d1 <- as.Date(to)
  as.Date(sample(as.integer(d0):as.integer(d1), n, replace = TRUE),
          origin = "1970-01-01")
}
# season from month (simple fictitious mapping)
.season <- function(dates) {
  m <- as.integer(format(dates, "%m"))
  out <- character(length(m))
  out[m %in% c(12, 1, 2, 3)] <- "Cold"
  out[m %in% c(6, 7, 8, 9)]  <- "Warm"
  out[m %in% c(4, 5)]        <- "Cold-to-warm"
  out[m %in% c(10, 11)]      <- "Warm-to-cold"
  factor(out, levels = SEASONS)
}
# length -> weight via allometric law W = a L^b (+ lognormal noise)
.lw <- function(L, a, b, cv = 0.12) {
  round(a * L^b * stats::rlnorm(length(L), 0, cv), 1)
}
# maturity ogive: P(mature) = plogis((L - L50)/d)
.mature <- function(L, L50, d = 1.5) stats::plogis((L - L50) / d)

## =====================================================================
## 1. CYMBIUM (gastropod) - bio + freq
## =====================================================================
# 4 species of the genus Cymbium; LCQ lengths ~ 3-40 cm
cym_species <- c("Cymbium cymbium", "Cymbium glans",
                 "Cymbium marmoratum", "Cymbium pepo")

.sim_cymbium_bio <- function(n = 3130) {
  esp <- factor(sample(cym_species, n, replace = TRUE,
                       prob = c(.42, .22, .12, .24)),
                levels = cym_species)
  # Linf by species
  Linf <- c("Cymbium cymbium" = 42, "Cymbium glans" = 38,
            "Cymbium marmoratum" = 34, "Cymbium pepo" = 30)[as.character(esp)]
  LCQ  <- round(pmin(pmax(Linf * (1 - exp(-stats::rgamma(n, 3, 3))), 4.1),
                     Linf), 1)
  dates <- .rdate(n)
  data.frame(
    landing_site = factor(sample(SITES, n, replace = TRUE), levels = SITES),
    Date         = dates,
    fishing_zone = factor(sample(ZONES, n, replace = TRUE), levels = ZONES),
    gear         = factor(sample(c("Bottom set net", "Trammel net"), n,
                          replace = TRUE),
                          levels = c("Bottom set net", "Trammel net")),
    season       = .season(dates),
    species      = esp,
    LCQ          = LCQ,
    weight       = as.integer(.lw(LCQ, a = 0.76, b = 2.63)),
    sex          = factor(sample(c("F", "M"), n, replace = TRUE),
                         levels = c("F", "M")),
    maturity     = factor(
      ifelse(stats::runif(n) < .mature(LCQ, L50 = 9.8),
             sample(c("M", "P"), n, replace = TRUE), "IM"),
      levels = c("IM", "M", "P")),
    stringsAsFactors = FALSE
  )
}

.sim_cymbium_freq <- function(n = 8243) {
  esp <- factor(sample(cym_species, n, replace = TRUE,
                       prob = c(.42, .22, .12, .24)),
                levels = cym_species)
  Linf <- c("Cymbium cymbium" = 42, "Cymbium glans" = 38,
            "Cymbium marmoratum" = 34, "Cymbium pepo" = 30)[as.character(esp)]
  LCQ  <- round(pmin(pmax(Linf * (1 - exp(-stats::rgamma(n, 3, 3))), 2.6),
                     Linf), 1)
  dates <- .rdate(n)
  data.frame(
    landing_site = factor(sample(SITES, n, replace = TRUE), levels = SITES),
    Date         = dates,
    fishing_zone = factor(sample(ZONES, n, replace = TRUE), levels = ZONES),
    gear         = factor(sample(c("Bottom set net", "Trammel net"), n,
                          replace = TRUE),
                          levels = c("Bottom set net", "Trammel net")),
    season       = .season(dates),
    species      = esp,
    LCQ          = LCQ,
    weight       = as.integer(.lw(LCQ, a = 0.76, b = 2.63)),
    sex          = factor(sample(c("F", "M"), n, replace = TRUE),
                         levels = c("F", "M")),
    stringsAsFactors = FALSE
  )
}

## =====================================================================
## 2. OCTOPUS (cephalopod) - bio + freq
## =====================================================================
.sim_octopus_bio <- function(n = 1682) {
  LM <- round(pmin(pmax(25 * (1 - exp(-stats::rgamma(n, 3, 3))), 4), 25))
  dates <- .rdate(n)
  poids <- .lw(LM, a = 0.55, b = 2.9, cv = 0.15)
  data.frame(
    landing_site = factor(sample(SITES[1:2], n, replace = TRUE),
                          levels = SITES[1:2]),
    Date         = dates,
    fishing_zone = factor(sample(ZONES, n, replace = TRUE), levels = ZONES),
    gear         = factor(rep("Pot", n), levels = "Pot"),
    season       = .season(dates),
    species      = factor(rep("Octopus vulgaris", n),
                          levels = "Octopus vulgaris"),
    LM           = LM,
    weight       = poids,
    sex          = factor(sample(c("F", "M"), n, replace = TRUE),
                         levels = c("F", "M")),
    maturity     = factor(sample(as.character(1:4), n, replace = TRUE),
                          levels = as.character(1:4)),
    total_weight_kg = ifelse(stats::runif(n) < .07,
                             round(stats::runif(n, 10, 65), 1), NA_real_),
    sample_weight_kg = ifelse(stats::runif(n) < .03,
                             round(stats::runif(n, 3, 5), 2), NA_real_),
    stringsAsFactors = FALSE
  )
}

.sim_octopus_freq <- function(n = 1689) {
  LM <- round(pmin(pmax(25 * (1 - exp(-stats::rgamma(n, 3, 3))), 4), 25))
  dates <- .rdate(n)
  data.frame(
    landing_site = factor(sample(SITES[1:2], n, replace = TRUE),
                          levels = SITES[1:2]),
    Date         = dates,
    fishing_zone = factor(sample(ZONES, n, replace = TRUE), levels = ZONES),
    gear         = factor(rep("Pot", n), levels = "Pot"),
    season       = .season(dates),
    species      = factor(rep("Octopus vulgaris", n),
                          levels = "Octopus vulgaris"),
    LM           = LM,
    weight       = .lw(LM, a = 0.55, b = 2.9, cv = 0.15),
    sex          = factor(sample(c("F", "M"), n, replace = TRUE),
                         levels = c("F", "M")),
    stringsAsFactors = FALSE
  )
}

## =====================================================================
## 3. PENAEUS (shrimp) - bio (with gonad_weight) + freq
## =====================================================================
.sim_penaeus_bio <- function(n = 6286) {
  LCT <- round(pmin(pmax(44 * (1 - exp(-stats::rgamma(n, 3.5, 3))), 8), 44))
  dates <- .rdate(n)
  sexe  <- factor(sample(c("F", "M"), n, replace = TRUE), levels = c("F", "M"))
  # gonad RGS: females only, seasonal (peak in late cold season)
  moisd <- as.integer(format(dates, "%m"))
  gon_base <- 0.03 + 0.02 * cos((moisd - 2) / 12 * 2 * pi)
  gonad_weight <- ifelse(sexe == "F" & stats::runif(n) < .8,
                     round(pmax(gon_base * (LCT / 20) *
                               stats::rlnorm(n, 0, 0.6), 0.01), 2),
                     NA_real_)
  data.frame(
    landing_site = factor(sample(SITES, n, replace = TRUE), levels = SITES),
    Date         = dates,
    fishing_zone = factor(sample(ZONES, n, replace = TRUE), levels = ZONES),
    gear         = factor(sample(c("Bottom set net", "Fixed net", "Barrier trap"),
                          n, replace = TRUE),
                          levels = c("Bottom set net", "Fixed net", "Barrier trap")),
    season       = .season(dates),
    species      = factor(rep("Penaeus notialis", n),
                          levels = "Penaeus notialis"),
    LCT          = LCT,
    weight       = .lw(LCT, a = 0.01, b = 2.9, cv = 0.12),
    sex          = sexe,
    maturity     = factor(sample(as.character(1:4), n, replace = TRUE),
                          levels = as.character(1:4)),
    gonad_weight = gonad_weight,
    total_weight_kg = ifelse(stats::runif(n) < .95,
                             round(stats::runif(n, 3, 12), 1), NA_real_),
    sample_weight_kg = ifelse(stats::runif(n) < .88,
                             round(stats::runif(n, 0.7, 6.1), 1), NA_real_),
    stringsAsFactors = FALSE
  )
}

.sim_penaeus_freq <- function(n = 31043) {
  LCT <- round(pmin(pmax(44 * (1 - exp(-stats::rgamma(n, 3.5, 3))), 6), 44))
  LCT[sample.int(n, 3)] <- NA          # a few realistic missing values
  dates <- .rdate(n)
  data.frame(
    landing_site = factor(sample(SITES, n, replace = TRUE), levels = SITES),
    Date         = dates,
    fishing_zone = factor(sample(ZONES, n, replace = TRUE), levels = ZONES),
    gear         = factor(sample(c("Bottom set net", "Fixed net", "Barrier trap"),
                          n, replace = TRUE),
                          levels = c("Bottom set net", "Fixed net", "Barrier trap")),
    season       = .season(dates),
    species      = factor(rep("Penaeus notialis", n),
                          levels = "Penaeus notialis"),
    LCT          = LCT,
    sex          = factor(sample(c("F", "M"), n, replace = TRUE),
                         levels = c("F", "M")),
    stringsAsFactors = FALSE
  )
}

## =====================================================================
## 4. ANNUAL SERIES - total catches + CPUE (fictitious)
## =====================================================================
.sim_captures <- function() {
  an <- 1970:2024
  data.frame(
    year    = as.integer(an),
    Cymbium = round(pmax(50, 8000 + 4000 * sin((an - 1970) / 8) +
                         stats::rnorm(length(an), 0, 1500)), 1),
    Octopus = round(pmax(200, 6000 + 3000 * cos((an - 1970) / 6) +
                         stats::rnorm(length(an), 0, 1200)), 1),
    Penaeus = round(pmax(1000, 3500 + 1000 * sin((an - 1970) / 10) +
                         stats::rnorm(length(an), 0, 700)), 1),
    stringsAsFactors = FALSE
  )
}

.sim_cpue <- function() {
  an <- 1978:2024
  mkna <- function(v, p) { v[sample.int(length(v), round(p * length(v)))] <- NA; v }
  data.frame(
    year    = as.integer(an),
    Cymbium = round(pmax(0.04, 1.2 - (an - 1978) / 60 +
                         stats::rnorm(length(an), 0, 0.25)), 3),
    Octopus = round(mkna(pmax(0.19, 1.0 + stats::rnorm(length(an), 0, 0.4)), .2), 3),
    Penaeus = round(mkna(pmax(0.001, 1.1 + stats::rnorm(length(an), 0, 0.4)), .4), 3),
    stringsAsFactors = FALSE
  )
}

## =====================================================================
## Generate + save
## =====================================================================
cymbium_bio      <- .sim_cymbium_bio()
cymbium_freq     <- .sim_cymbium_freq()
octopus_bio      <- .sim_octopus_bio()
octopus_freq     <- .sim_octopus_freq()
penaeus_bio      <- .sim_penaeus_bio()
penaeus_freq     <- .sim_penaeus_freq()
total_catches <- .sim_captures()
annual_cpue    <- .sim_cpue()

usethis::use_data(cymbium_bio, cymbium_freq, octopus_bio, octopus_freq,
                  penaeus_bio, penaeus_freq, total_catches, annual_cpue,
                  overwrite = TRUE, compress = "xz")

## ---- demonstration CSV (fictitious extracts) ----
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(
  data.frame(length = penaeus_bio$LCT, weight = penaeus_bio$weight,
             sex = penaeus_bio$sex, maturity = penaeus_bio$maturity),
  "inst/extdata/penaeus_bio_demo.csv", row.names = FALSE)
utils::write.csv(
  data.frame(length = penaeus_freq$LCT,
             year  = as.integer(format(penaeus_freq$Date, "%Y")),
             month = as.integer(format(penaeus_freq$Date, "%m")),
             sex   = penaeus_freq$sex),
  "inst/extdata/penaeus_lengths_demo.csv", row.names = FALSE)

message("Fictitious data generated: 8 .rda datasets + 2 demonstration CSVs.")
