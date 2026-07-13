# =============================================================================
# stockflow :: 03_cpue_standardisation.R
# Phase 2C — CPUE standardisation (abundance index)
#   Lognormal GLMM : log(CPUE) ~ year + month + region + gear  (+ random effects)
#   Output : standardised annual index (mean = 1) with CI, back-transformed (bias).
#
# Validated decision : standardisation by GLM (vs nominal CPUE).
# Prerequisite : 01_clean_data.R  (data/processed/base_cpue_clean.csv).
# =============================================================================
suppressPackageStartupMessages({library(tidyverse); library(glmmTMB); library(broom.mixed)})

proc <- "data/processed"; tab <- "output/tables"; fig <- "output/figures"

cpue <- read_csv(file.path(proc,"base_cpue_clean.csv"), show_col_types = FALSE) |>
  filter(capture_port > 0, sorties_port > 0) |>
  mutate(log_cpue = log(capture_port/sorties_port),
         fyear  = factor(annee), fmonth = factor(mois),
         region = factor(region), engin = factor(engin))

standardise <- function(df) {
  # fixed effects year/month/gear ; region as random effect if >2 levels
  fixed <- "log_cpue ~ fyear"
  if (nlevels(droplevels(df$fmonth)) > 1) fixed <- paste(fixed, "+ fmonth")
  if (nlevels(droplevels(df$engin))  > 1) fixed <- paste(fixed, "+ engin")
  reff  <- if (nlevels(droplevels(df$region)) > 2) "+ (1|region)" else ""
  m <- glmmTMB::glmmTMB(as.formula(paste(fixed, reff)), data = df, family = gaussian())

  yrs   <- levels(droplevels(df$fyear))
  base  <- df |> count(fmonth, engin, region, sort = TRUE) |> slice(1)
  newd  <- tibble(fyear = factor(yrs, levels = levels(df$fyear)),
                  fmonth = base$fmonth, engin = base$engin, region = base$region)
  pr    <- predict(m, newd, se.fit = TRUE, re.form = NA, allow.new.levels = TRUE)
  s2    <- sigma(m)^2
  tibble(annee = as.integer(yrs),
         index = exp(pr$fit + s2/2),                      # lognormal bias correction
         lci   = exp(pr$fit - 1.96*pr$se.fit + s2/2),
         uci   = exp(pr$fit + 1.96*pr$se.fit + s2/2)) |>
    mutate(index_scaled = index/mean(index))
}

res <- cpue |> group_split(nom_scientifique) |>
  set_names(map_chr(group_split(cpue, nom_scientifique), ~unique(.x$nom_scientifique)))

walk(names(res), function(sp) {
  idx <- standardise(res[[sp]])
  write_csv(idx, file.path(tab, sprintf("cpue_std_%s.csv", tolower(word(sp,1)))))
  message(sp, " : index ", min(idx$annee), "-", max(idx$annee))
})

# figure
all_idx <- imap_dfr(res, ~ standardise(.x) |> mutate(espece = .y))
p <- ggplot(all_idx, aes(annee, index_scaled)) +
  geom_ribbon(aes(ymin = lci/index*index_scaled, ymax = uci/index*index_scaled),
              alpha = .2, fill = "#2c7fb8") +
  geom_line(colour = "#12507b") + geom_point(size = 1) +
  facet_wrap(~espece, scales = "free_y") +
  labs(title = "Standardised CPUE (lognormal GLMM)", x = NULL, y = "Index (mean = 1)") +
  theme_bw()
ggsave(file.path(fig, "fig7_cpue_standardisee.png"), p, width = 14, height = 4.5, dpi = 130)
message("\n=== Phase 2C completed ===")