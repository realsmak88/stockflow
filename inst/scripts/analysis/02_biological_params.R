# =============================================================================
# stockflow :: 02_biological_params.R
# Phase 2A — Biological parameters
#   - Length-weight relationship  W = a·L^b   (log-log regression)
#   - Logistic maturity ogive -> L50, L95   (binomial GLM)
#   - Cymbium proratas (4 species) from samples -> catch split
#
# Validated decisions: mature = stages >=3 (octopus/shrimp), M+P (volutes).
# Prerequisite: 01_clean_data.R executed (data/processed/*.csv).
# =============================================================================
suppressPackageStartupMessages({library(tidyverse); library(FSA)})

proc <- "data/processed"; tab <- "output/tables"; fig <- "output/figures"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Length-weight relationship ------------------------------------------------
fit_lw <- function(df, lcol, wcol, stock) {
  d <- df |> filter(.data[[lcol]] > 0, .data[[wcol]] > 0,
                    !isTRUE(flag_lw_outlier)) |>
    transmute(L = .data[[lcol]], W = .data[[wcol]])
  m <- lm(log(W) ~ log(L), data = d)
  tibble(stock = stock, unit_L = lcol,
         a = exp(coef(m)[1]), b = coef(m)[2],
         r2 = summary(m)$r.squared, n = nrow(d),
         Lmin = min(d$L), Lmax = max(d$L))
}

# ---- 2. Maturity ogive (binomial logit GLM) -------------------------------
fit_maturity <- function(df, lcol, stock, scale = c("numeric","cymbium")) {
  scale <- match.arg(scale)
  d <- df |> mutate(L = as.numeric(.data[[lcol]]))
  if (scale == "numeric") {
    d <- d |> mutate(mat = as.integer(suppressWarnings(as.numeric(maturity)) >= 3)) |>
      filter(!is.na(mat), !is.na(L))
  } else {
    mm <- str_to_upper(str_squish(d$maturity))
    d <- d |> mutate(mat = case_when(mm %in% c("M","P") ~ 1L,
                                     mm %in% c("IM","I") ~ 0L, TRUE ~ NA_integer_)) |>
      filter(!is.na(mat), !is.na(L))
  }
  g  <- glm(mat ~ L, data = d, family = binomial(link = "logit"))
  cf <- coef(g)
  L50 <- -cf[1]/cf[2]
  L95 <- (log(0.95/0.05) - cf[1])/cf[2]
  tibble(stock = stock, L50 = L50, L95 = L95,
         n = nrow(d), prop_mature = mean(d$mat))
}

lw_tbl <- bind_rows(
  fit_lw(read_csv(file.path(proc,"octopus_biologie_clean.csv"), show_col_types = FALSE), "LM",  "weight", "Octopus vulgaris"),
  fit_lw(read_csv(file.path(proc,"penaeus_biologie_clean.csv"), show_col_types = FALSE), "LCT", "weight", "Penaeus notialis"),
  read_csv(file.path(proc,"cymbium_biologie_clean.csv"), show_col_types = FALSE) |>
    group_split(species) |>
    map_dfr(~ fit_lw(.x, "LCQ", "weight", unique(.x$species)))
)
write_csv(lw_tbl, file.path(tab, "params_length_weight.csv"))
print(lw_tbl)

mat_tbl <- bind_rows(
  fit_maturity(read_csv(file.path(proc,"octopus_biologie_clean.csv"), show_col_types = FALSE), "LM",  "Octopus vulgaris", "numeric"),
  fit_maturity(read_csv(file.path(proc,"penaeus_biologie_clean.csv"), show_col_types = FALSE), "LCT", "Penaeus notialis", "numeric"),
  read_csv(file.path(proc,"cymbium_biologie_clean.csv"), show_col_types = FALSE) |>
    group_split(species) |>
    map_dfr(~ fit_maturity(.x, "LCQ", unique(.x$species), "cymbium"))
)
write_csv(mat_tbl, file.path(tab, "params_maturity.csv"))
print(mat_tbl)

# ---- 3. Cymbium proratas (4 species) & catch split --------------------
cy_freq <- read_csv(file.path(proc,"cymbium_frequence_clean.csv"), show_col_types = FALSE)
cy_bio  <- read_csv(file.path(proc,"cymbium_biologie_clean.csv"),  show_col_types = FALSE)
lw_cy   <- lw_tbl |> filter(str_starts(stock,"Cymbium")) |> select(stock,a,b)

allc <- bind_rows(cy_freq |> select(species, LCQ),
                  cy_bio  |> select(species, LCQ)) |>
  filter(!is.na(LCQ), !is.na(species)) |>
  left_join(lw_cy, by = c("species"="stock")) |>
  mutate(w = a * LCQ^b)

prorata <- allc |> group_by(espece = species) |>
  summarise(prorata_nombre = n(), poids = sum(w, na.rm = TRUE), .groups="drop") |>
  mutate(prorata_nombre = prorata_nombre/sum(prorata_nombre),
         prorata_poids  = poids/sum(poids)) |> select(-poids)
write_csv(prorata, file.path(tab, "cymbium_proratas_echantillons.csv"))
print(prorata)

catch <- read_csv(file.path(proc,"total_catches_annual.csv"), show_col_types = FALSE)
split_cy <- catch |> filter(espece == "Cymbium spp") |>
  tidyr::crossing(prorata |> select(sp = espece, p = prorata_poids)) |>
  transmute(annee, espece = sp, capture_t = capture_t * p)
write_csv(split_cy, file.path(proc, "captures_cymbium_par_espece.csv"))

# ---- 4. Figures --------------------------------------------------------------
logistic_p <- function(L, L50, L95) 1/(1+exp(-log(19)*(L-L50)/(L95-L50)))
p_mat <- mat_tbl |> filter(!is.na(L50)) |>
  mutate(data = purrr::pmap(list(L50,L95), function(a,b)
    tibble(L = seq(0, b*1.6, length.out = 200), p = logistic_p(L,a,b)))) |>
  tidyr::unnest(data) |>
  ggplot(aes(L, p)) + geom_line(colour = "#c0392b", linewidth = 1) +
  facet_wrap(~stock, scales = "free_x") +
  labs(title = "Maturity ogives", x = "Length", y = "P(mature)") +
  theme_bw()
ggsave(file.path(fig, "fig5_ogives_maturite.png"), p_mat, width = 12, height = 7, dpi = 130)

message("\n=== Phase 2A completed: biological parameters in output/tables/ ===")