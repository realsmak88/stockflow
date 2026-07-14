## Tests du module graphique (plots.R)
## Chaque fonction doit renvoyer un objet ggplot constructible sans erreur.

## ---- jeux de donnees synthetiques legers ----------------------------
set.seed(42)
n <- 300
bio <- data.frame(
  L   = c(stats::rnorm(n, 12, 3), stats::rnorm(n, 20, 4)),
  W   = NA_real_,
  sex = sample(c("M", "F"), 2 * n, replace = TRUE),
  Mat  = sample(c("IM", "M", "P"), 2 * n, replace = TRUE),
  Date = as.Date("2025-01-01") + sample(0:364, 2 * n, replace = TRUE)
)
bio$L <- pmax(bio$L, 1)
bio$W <- 0.01 * bio$L^3 * exp(stats::rnorm(2 * n, 0, 0.1))

is_gg <- function(p) inherits(p, "ggplot")
builds <- function(p) {
  ok <- tryCatch({ ggplot2::ggplot_build(p); TRUE }, error = function(e) FALSE)
  ok
}


test_that("plot_length_frequency returns a ggplot (simple and month facet)", {
  p1 <- plot_length_frequency(bio, "L", sex_col = "sex")
  expect_true(is_gg(p1)); expect_true(builds(p1))
  p2 <- plot_length_frequency(bio, "L", date_col = "Date")
  expect_true(is_gg(p2)); expect_true(builds(p2))
  expect_error(plot_length_frequency(bio, "inexistante"), "missing")
})


test_that("plot_length_weight fits and plots the power law", {
  ## scale 'power' (natural curve) and 'log' -> one ggplot each
  p_pow <- plot_length_weight(bio, "L", "W", scale = "power")
  expect_true(is_gg(p_pow)); expect_true(builds(p_pow))
  p_log <- plot_length_weight(bio, "L", "W", scale = "log")
  expect_true(is_gg(p_log)); expect_true(builds(p_log))

  ## scale 'both' -> assembling both panels (patchwork if available)
  p_both <- plot_length_weight(bio, "L", "W", scale = "both")
  expect_true(inherits(p_both, "ggplot") || inherits(p_both, "patchwork") ||
              is.list(p_both))

  ## backward compatibility of the old log_scale argument
  expect_true(is_gg(plot_length_weight(bio, "L", "W", log_scale = TRUE)))
  expect_true(is_gg(plot_length_weight(bio, "L", "W", log_scale = FALSE)))

  ## the reported exponent b must be close to 3 (data simulated as L^3)
  m <- stats::lm(log(bio$W) ~ log(bio$L))
  expect_equal(unname(stats::coef(m)[2]), 3, tolerance = 0.1)
})


test_that("plot_length_weight colours by sex and filters M/F", {
  ## both sexes: coloured scatter + one fit per sex
  p_both <- plot_length_weight(bio, "L", "W", sex_col = "sex", scale = "power")
  expect_true(is_gg(p_both)); expect_true(builds(p_both))

  ## single sex: restricts the data and still builds
  p_f <- plot_length_weight(bio, "L", "W", sex_col = "sex", sex = "F",
                            scale = "power")
  expect_true(is_gg(p_f)); expect_true(builds(p_f))
  p_m <- plot_length_weight(bio, "L", "W", sex_col = "sex", sex = "M",
                            scale = "log")
  expect_true(is_gg(p_m)); expect_true(builds(p_m))

  ## sex != "both" without sex_col is an error
  expect_error(plot_length_weight(bio, "L", "W", sex = "F"), "sex_col")
})


test_that("plot_vbgf accepts one or several parameter sets", {
  skip_if_not_installed("TropFishR")
  p1 <- plot_vbgf(Linf = 42, K = 0.45)
  expect_true(is_gg(p1)); expect_true(builds(p1))
  p2 <- plot_vbgf(Linf = c(42, 45), K = c(0.45, 0.38),
                  labels = c("SA", "GA"))
  expect_true(is_gg(p2)); expect_true(builds(p2))
})


test_that("plot_maturity_ogive extracts a plausible L50", {
  p <- plot_maturity_ogive(bio, "L", "Mat", mature_codes = c("M", "P"))
  expect_true(is_gg(p)); expect_true(builds(p))
  expect_error(
    plot_maturity_ogive(bio, "L", "Mat", mature_codes = NULL),
    "mature_codes")
})


test_that("plot_sex_ratio works globally and by size class", {
  p1 <- plot_sex_ratio(bio, "sex")
  expect_true(is_gg(p1)); expect_true(builds(p1))
  p2 <- plot_sex_ratio(bio, "sex", "L")
  expect_true(is_gg(p2)); expect_true(builds(p2))
})


test_that("plot_condition plots Le Cren and Fulton (with/without month)", {
  p1 <- plot_condition(bio, "L", "W", index = "lecren")
  expect_true(is_gg(p1)); expect_true(builds(p1))
  p2 <- plot_condition(bio, "L", "W", index = "fulton", date_col = "Date")
  expect_true(is_gg(p2)); expect_true(builds(p2))
})


test_that("plot_length_evolution plots the temporal evolution of lengths", {
  set.seed(11)
  ne <- 500
  bt <- data.frame(
    L    = stats::rnorm(ne, 20, 4),
    sex = sample(c("F", "M"), ne, replace = TRUE),
    Date = as.Date("2024-01-01") + sample(0:729, ne, replace = TRUE)
  )
  ## the 4 period x type combinations
  expect_true(builds(plot_length_evolution(bt, "L", "Date",
                                           period = "month", type = "histogram")))
  expect_true(builds(plot_length_evolution(bt, "L", "Date",
                                           period = "month", type = "density")))
  expect_true(builds(plot_length_evolution(bt, "L", "Date",
                                           period = "year", type = "histogram")))
  expect_true(builds(plot_length_evolution(bt, "L", "Date",
                                           period = "year", type = "density")))
  ## histogram colored by sex
  expect_true(is_gg(plot_length_evolution(bt, "L", "Date", period = "month",
                                          type = "histogram", sex_col = "sex")))
  ## invalid arguments rejected
  expect_error(plot_length_evolution(bt, "L", "Date", period = "week"))
  expect_error(plot_length_evolution(bt, "L", "Date", type = "boxplot"))
})


test_that("plot_gsi computes the GSI and plots the monthly evolution", {
  set.seed(7)
  ng <- 400
  bg <- data.frame(
    Pg   = pmax(stats::rnorm(ng, 0.05, 0.02), 0.001),
    W    = stats::runif(ng, 5, 40),
    sex = sample(c("F", "M"), ng, replace = TRUE),
    Date = as.Date("2025-01-01") + sample(0:364, ng, replace = TRUE)
  )
  ## monthly evolution, restricted to females
  p1 <- plot_gsi(bg, "Pg", "W", date_col = "Date",
                 sex_col = "sex", sex_filter = "F")
  expect_true(is_gg(p1)); expect_true(builds(p1))
  ## overall distribution (without date)
  p2 <- plot_gsi(bg, "Pg", "W")
  expect_true(is_gg(p2)); expect_true(builds(p2))

  ## the anti-artifact filter discards implausible values
  bg2 <- bg; bg2$Pg[1] <- bg2$W[1] * 10   # GSI = 1000 %
  expect_true(is_gg(plot_gsi(bg2, "Pg", "W", max_gsi = 20)))

  ## error if too few usable values
  bg3 <- bg[1:2, ]
  expect_error(plot_gsi(bg3, "Pg", "W"), "Too few")
})


test_that("plot_catch_series plots wide series", {
  df <- data.frame(year = 2000:2010,
                   A = runif(11, 100, 200),
                   B = runif(11, 50, 150))
  p <- plot_catch_series(df, "year", ylab = "Captures (t)")
  expect_true(is_gg(p)); expect_true(builds(p))
})


test_that("S3 plot methods for management objects return a ggplot", {
  skip_if_not_installed("TropFishR")
  skip_if_not_installed("FSA")

  ## multi-method M
  m <- estimate_M_all(Linf = 42, K = 0.45, temp = 25, tmax = 8)
  expect_true(is_gg(plot(m))); expect_true(builds(plot(m)))

  ## fitted length-weight relationship
  lw <- fit_length_weight(data.frame(L = bio$L, W = bio$W),
                          length_col = "L", weight_col = "W")
  expect_true(is_gg(plot(lw)))

  ## TAC + quotas
  tac <- compute_tac("ratio_E", catch_current = 5000,
                     E_current = 0.7, E_target = 0.6, buffer = 0.1)
  expect_true(is_gg(plot(tac))); expect_true(builds(plot(tac)))
  q <- allocate_quota(tac, c(A = 0.5, B = 0.3, C = 0.2))
  expect_true(is_gg(plot(q))); expect_true(builds(plot(q)))

  ## management sizes
  sa <- recommend_sizes(c(10, 13, 15, 18), Linf = 42,
                        L50_maturity = 13, MK = 1.4)
  expect_true(is_gg(plot(sa))); expect_true(builds(plot(sa)))

  ## mortality decomposition
  mort <- structure(list(Z = list(Z = 0.9), M_consensus = 0.35),
                    class = "FishStockMortality")
  expect_true(is_gg(plot(mort))); expect_true(builds(plot(mort)))
})
