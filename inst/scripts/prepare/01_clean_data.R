# =============================================================================
# stockflow :: 01_clean_data.R
# Phase 1 — Data cleaning & quality control
#   - Biology & length frequencies: Octopus, Penaeus, Cymbium
#   - Annual total catches (example input)
#   - Artisanal effort / CPUE (base_cpue)
#
# Exactly reproduces the cleaning logic validated in Phase 1.
# Writes clean data to data/processed/ and a QC log.
#
# Author: Kamarel Ba — usage: source("R/prepare/01_clean_data.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tidyr, stringr, readr, ggplot2
  library(readxl)
  library(janitor)
  library(lubridate)
})

# --- paths (relative to the stockflow project root) ---------------------
dir_raw  <- "."                       # the .xlsx files are at the project root
raw_path <- function(f) file.path("..", f)   # go up one level from stockflow/
proc_dir <- "data/processed"; qc_dir <- "output/qc"
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)

qc_log <- character(0)
log_qc <- function(section, msg) {
  line <- sprintf("[%s] %s", section, msg)
  message(line); qc_log <<- c(qc_log, line)
}

# --- harmonization helpers -------------------------------------------------
strip_txt <- function(x) str_squish(as.character(x))

norm_saison <- function(x) {
  k <- x |> strip_txt() |> str_to_lower() |> str_replace_all("\\s*-\\s*", "-")
  dplyr::case_when(
    k %in% c("chaude") ~ "Chaude",
    k %in% c("froide") ~ "Froide",
    str_detect(k, "chaude-froide") ~ "Transition chaude-froide",
    str_detect(k, "froide-chaude") ~ "Transition froide-chaude",
    TRUE ~ strip_txt(x)
  )
}

norm_site <- function(x) {
  # Example normalization of site labels (fictitious names).
  k <- str_to_lower(strip_txt(x))
  m <- c("debarcadere a" = "Debarcadere A", "debarcadere b" = "Debarcadere B",
         "debarcadere c" = "Debarcadere C")
  out <- unname(m[k]); ifelse(is.na(out), str_to_title(strip_txt(x)), out)
}

norm_engin <- function(x) {
  k <- str_to_lower(strip_txt(x))
  dplyr::case_when(
    k %in% c("killi", "filet killi") ~ "Killi",
    k == "fdf"          ~ "FDF",
    k == "filet fixe"   ~ "Filet fixe",
    k %in% c("trammel net","tremail") ~ "Trammel net",
    k == "lpo"          ~ "LPO",
    TRUE ~ strip_txt(x)
  )
}

norm_lieu <- function(x) {
  x |> strip_txt() |> str_replace_all("\\s*-\\s*", " ") |>
    str_to_lower() |>
    str_replace_all(c("^zone noir$"="zone noire", "^gopp$"="gop")) |>
    str_to_title()
}

# fixes numeric typos: "13..5"->13.5 ; ".19.1"->19.1 ; "1,2"->1.2
parse_num_typo <- function(v) {
  s <- as.character(v) |> str_trim() |> str_replace_all(",", ".") |>
       str_replace("^\\.", "") |> str_replace_all("\\.\\.", ".")
  # keep at most 2 segments
  s <- vapply(str_split(s, "\\."), function(p)
        if (length(p) > 2) paste0(p[1], ".", p[2]) else paste(p, collapse="."),
        character(1))
  suppressWarnings(as.numeric(s))
}

# length-weight outlier flag: log-log residuals W=aL^b, |res-med| > 4*MAD
flag_lw_outlier <- function(df, lcol, wcol, by = NULL) {
  df$flag_lw_outlier <- FALSE
  grp <- if (is.null(by)) list(seq_len(nrow(df))) else split(seq_len(nrow(df)), df[[by]])
  for (idx in grp) {
    sub <- df[idx, ]
    ok <- which(sub[[lcol]] > 0 & sub[[wcol]] > 0 &
                is.finite(sub[[lcol]]) & is.finite(sub[[wcol]]))
    if (length(ok) < 20) next
    x <- log(sub[[lcol]][ok]); y <- log(sub[[wcol]][ok])
    cf <- coef(lm(y ~ x)); res <- y - (cf[1] + cf[2]*x)
    med <- median(res); mad_ <- mad(res, constant = 1.4826) + 1e-9
    out <- abs(res - med) > 4*mad_
    df$flag_lw_outlier[idx[ok[out]]] <- TRUE
  }
  df
}

# =============================================================================
# 1. BIOLOGICAL DATA
# =============================================================================
clean_bio <- function(file, species_fixed = NULL, lcol) {
  sheets <- c("biologie", "frequence")
  res <- list()
  for (sh in sheets) {
    df <- read_excel(raw_path(file), sheet = sh) |> rename_with(str_squish)
    n0 <- nrow(df)
    if (!is.null(species_fixed)) df$species <- species_fixed
    df <- df |>
      mutate(
        landing_site = norm_site(landing_site),
        fishing_zone        = norm_lieu(fishing_zone),
        gear       = norm_engin(gear),
        season            = norm_saison(season),
        sex              = na_if(str_to_upper(strip_txt(sex)), "NAN"),
        !!lcol            := parse_num_typo(.data[[lcol]])
      )
    if ("weight" %in% names(df)) df <- mutate(df, weight = suppressWarnings(as.numeric(weight)))

    # -- Octopus weight unit correction (rows in kg, weight<10 -> g) -------
    if (grepl("octopus", file) && "weight" %in% names(df)) {
      kg <- which(df$weight < 10)
      if (length(kg)) {
        dts <- paste(unique(as.Date(df$Date[kg])), collapse=", ")
        df$weight[kg] <- df$weight[kg] * 1000
        log_qc("OCTOPUS", sprintf("%s : %d weights kg->g (dates %s)", sh, length(kg), dts))
      }
    }
    # -- length outliers -> NA (Penaeus LCT plausible cephalothorax 7-50 mm)
    if (grepl("penaeus", file)) {
      ab <- which(df[[lcol]] < 7 | df[[lcol]] > 50)
      if (length(ab)) { log_qc("PENAEUS", sprintf("%s : %d LCT outside [7,50] mm -> NA", sh, length(ab)))
                        df[[lcol]][ab] <- NA_real_ }
    }
    # -- Cymbium Maturity: IM/M/P in uppercase --------------------------------
    if (grepl("cymbium", file) && "maturity" %in% names(df))
      df$maturity <- str_to_upper(strip_txt(df$maturity))

    # -- exact duplicates: only in 'biologie' (frequence = legitimate repeated measures)
    if (sh == "biologie") {
      d <- sum(duplicated(df)); df <- distinct(df)
      log_qc(toupper(gsub("\\D","",file)), sprintf("biologie : %d exact duplicates removed", d))
    }
    df$flag_lw_outlier <- FALSE
    res[[sh]] <- df
    log_qc(basename(file), sprintf("%s : %d -> %d rows", sh, n0, nrow(df)))
  }
  res
}

oc <- clean_bio("data_octopus.xlsx", "Octopus vulgaris", "LM")
pn <- clean_bio("data_penaeus.xlsx", "Penaeus notialis", "LCT")
cy <- clean_bio("data_cymbium.xlsx", NULL,               "LCQ")

# length-weight flags (biologie sheet)
oc$biologie <- flag_lw_outlier(oc$biologie, "LM",  "weight")
pn$biologie <- flag_lw_outlier(pn$biologie, "LCT", "weight")
cy$biologie <- flag_lw_outlier(cy$biologie, "LCQ", "weight", by = "species")

# export
walk2(list(oc, pn, cy), c("octopus","penaeus","cymbium"), function(dd, nm) {
  walk(c("biologie","frequence"), function(sh) {
    p <- file.path(proc_dir, sprintf("%s_%s_clean.csv", nm, sh))
    write_csv(dd[[sh]], p); log_qc("EXPORT", sprintf("%s (%d rows)", p, nrow(dd[[sh]])))
  })
})

# =============================================================================
# 2. ANNUAL TOTAL CATCHES
# =============================================================================
tot <- read_excel(raw_path("catch_data.xlsx"), sheet = "totales")
catch_long <- tot[, 1:4]
names(catch_long) <- c("annee","Cymbium spp","Octopus vulgaris","Penaeus notialis")
catch_long <- catch_long |>
  filter(!is.na(annee)) |> mutate(annee = as.integer(annee)) |>
  pivot_longer(-annee, names_to = "espece", values_to = "capture_t") |>
  filter(!is.na(capture_t)) |> mutate(capture_t = as.numeric(capture_t)) |>
  arrange(espece, annee)
stopifnot(all(catch_long$capture_t >= 0))
write_csv(catch_long, file.path(proc_dir, "total_catches_annual.csv"))
log_qc("CATCH", sprintf("total catches: %d rows, %d-%d",
                        nrow(catch_long), min(catch_long$annee), max(catch_long$annee)))

# Cymbium DPM prorata (cymbium/glans/pepo) — marmoratum NOT included
pror <- tot[2:4, 6:7]; names(pror) <- c("espece","prorata_DPM_2017_2023")
pror <- pror |> filter(!is.na(espece)) |>
  mutate(prorata_DPM_2017_2023 = as.numeric(prorata_DPM_2017_2023))
write_csv(pror, file.path(proc_dir, "cymbium_prorata_DPM.csv"))
log_qc("CATCH", "Cymbium DPM prorata exported (WARNING: Cymbium marmoratum absent from the prorata)")

# =============================================================================
# 3. ARTISANAL EFFORT / CPUE
# =============================================================================
cpue <- read_excel(raw_path("base_cpue.xlsx")) |> rename_with(str_squish)
n0 <- nrow(cpue)
d  <- sum(duplicated(cpue)); cpue <- distinct(cpue)
log_qc("CPUE", sprintf("%d exact duplicates removed (%d -> %d)", d, n0, nrow(cpue)))
cpue <- cpue |>
  mutate(region = str_squish(region), port = str_squish(port),
         engin  = str_to_upper(str_squish(engin)),
         nom_scientifique = str_squish(nom_scientifique),
         cpue_rec = capture_port / sorties_port)   # CPUE per record
stopifnot(all(cpue$mois %in% 1:12), all(cpue$sorties_port > 0))
write_csv(cpue, file.path(proc_dir, "base_cpue_clean.csv"))
log_qc("CPUE", sprintf("cleaned: %d rows, %d-%d, gears %s",
                       nrow(cpue), min(cpue$annee), max(cpue$annee),
                       paste(sort(unique(cpue$engin)), collapse=",")))

# ----- PROVISIONAL annual index (mean of ratios weighted by catch) ------
# NB: effort aggregation method TO BE VALIDATED (see QC report).
cpue_index <- cpue |>
  group_by(espece = nom_scientifique, annee) |>
  summarise(n_rec = n(), capture_tot = sum(capture_port, na.rm = TRUE),
            cpue_moy  = mean(cpue_rec, na.rm = TRUE),
            cpue_pond = sum(capture_port * cpue_rec, na.rm = TRUE) /
                        sum(capture_port, na.rm = TRUE),
            .groups = "drop") |>
  arrange(espece, annee)
write_csv(cpue_index, file.path(proc_dir, "cpue_index_annuel_PROVISOIRE.csv"))
log_qc("CPUE", "PROVISIONAL annual index written (effort aggregation to be confirmed)")

# QC log
writeLines(qc_log, file.path(qc_dir, "journal_nettoyage_R.txt"))
message("\n=== Phase 1 (R) complete — clean data in data/processed/ ===")