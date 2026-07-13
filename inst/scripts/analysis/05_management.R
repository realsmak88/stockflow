# =============================================================================
# stockflow :: 05_management.R
# Phase 2E — Management measures
#   (1) TAC via harvest control rule (ICES-type hockey-stick HCR)
#   (2) Quota allocation by segment (gear) in proportion to historical catches
#   (3) Spatiotemporal closures (spawning & recruitment peaks)
#   (4) Technical measures: minimum landing size (TMD = L50)
# Prerequisites: 02, 03, 04 executed.
# =============================================================================
suppressPackageStartupMessages({library(tidyverse); library(lubridate)})
proc<-"data/processed"; tab<-"output/tables"; fig<-"output/figures"

brp <- read_csv(file.path(tab,"stock_assessment_BRP.csv"), show_col_types=FALSE)
mat <- read_csv(file.path(tab,"params_maturity.csv"),      show_col_types=FALSE)
lw  <- read_csv(file.path(tab,"params_length_weight.csv"), show_col_types=FALSE)

# ---- (1) TAC (HCR: ramp between Blim=0.2 and Btrigger=1.0 of Bmsy) -----------
hcr <- function(B_Bmsy, F_Fmsy, MSY, Blim=0.2, Btrig=1.0) {
  f_frac <- dplyr::case_when(B_Bmsy >= Btrig ~ 1,
                             B_Bmsy >  Blim  ~ (B_Bmsy-Blim)/(Btrig-Blim),
                             TRUE ~ 0)
  tac <- f_frac*MSY
  if (F_Fmsy > 1) tac <- min(tac, 0.9*MSY*B_Bmsy)   # precaution in case of overfishing
  tibble(f_frac = round(f_frac,2), TAC_t = round(tac,0))
}
tac <- brp |> rowwise() |>
  mutate(hcr(B_Bmsy, F_Fmsy, MSY)) |> ungroup() |>
  select(stock, MSY, B_Bmsy, F_Fmsy, f_frac, TAC_t)
write_csv(tac, file.path(tab,"TAC_par_stock.csv")); print(tac)

# ---- (2) Quotas by segment (gear) in proportion to history (last 5 years)
cp <- read_csv(file.path(proc,"base_cpue_clean.csv"), show_col_types=FALSE)
quotas <- map_dfr(tac$stock, function(sp){
  d <- cp |> filter(nom_scientifique==sp, annee >= max(annee)-4)
  sh <- d |> group_by(segment_engin=engin) |>
    summarise(capture=sum(capture_port), .groups="drop") |>
    mutate(part_historique = capture/sum(capture),
           TAC = tac$TAC_t[tac$stock==sp],
           quota_t = round(part_historique*TAC,1), stock=sp)
  sh |> select(stock, segment_engin, part_historique, quota_t) |> arrange(desc(quota_t))
})
write_csv(quotas, file.path(tab,"quotas_par_segment.csv")); print(quotas)
# NB. "By fleet/vessel owner" extension: replace 'engin' with the fleet/registration
#     key once the industrial data (logbooks) have been joined.

# ---- (3) Spatiotemporal closures --------------------------------------
bio_files <- list(`Octopus vulgaris`=c("octopus","LM","num"),
                  `Penaeus notialis`=c("penaeus","LCT","num"),
                  `Cymbium spp`     =c("cymbium","LCQ","cym"))
season_tab <- imap_dfr(bio_files, function(v, sp){
  b <- read_csv(file.path(proc, sprintf("%s_biologie_clean.csv", v[1])), show_col_types=FALSE) |>
    mutate(mois = month(as.Date(Date)))
  if (v[3]=="num") b <- b |> mutate(mature=as.numeric(maturity)>=3, ponte=as.numeric(maturity)>=4)
  else b <- b |> mutate(mm=str_to_upper(maturity), mature=mm %in% c("M","P"), ponte=mm=="P")
  L50 <- if (sp=="Cymbium spp") mean(mat$L50[str_starts(mat$stock,"Cymbium")]) else mat$L50[mat$stock==sp]
  b <- b |> mutate(juv = .data[[v[2]]] < L50)
  b |> group_by(stock=sp, mois) |>
    summarise(n=n(), p_mature=mean(mature), p_ponte=mean(ponte), p_juv=mean(juv), .groups="drop")
})
fermetures <- season_tab |> group_by(stock) |>
  summarise(mois_pic_ponte = mois[which.max(p_ponte)],
            mois_pic_juveniles = mois[which.max(p_juv)],
            p_ponte_max = max(p_ponte), p_juv_max = max(p_juv), .groups="drop")
write_csv(fermetures, file.path(tab,"fermetures_reco.csv")); print(fermetures)

# ---- (4) Technical measures: TMD = L50 ------------------------------------
tmd <- mat |> left_join(lw, by="stock") |>
  transmute(stock, L50, TMD_reco = round(L50,1), poids_a_TMD_g = round(a*L50^b,1))
write_csv(tmd, file.path(tab,"mesures_techniques_TMD.csv")); print(tmd)

message("\n=== Phase 2E complete: TAC, quotas, closures, TMD in output/tables/ ===")