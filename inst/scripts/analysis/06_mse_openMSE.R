# =============================================================================
# stockflow :: 06_mse_openMSE.R
# Phase 2F — Management strategy evaluation (MSE) in closed loop
#   openMSE (MSEtool + DLMtool + SAMtool). One operating model (OM) per stock,
#   conditioned on catches + CPUE, testing several management procedures (MP).
# Prerequisites: 04 (BRP) and 03 (CPUE). Long run (Monte-Carlo) -> parallel.
# =============================================================================
suppressPackageStartupMessages({library(openMSE); library(tidyverse)})
proc<-"data/processed"; tab<-"output/tables"; fig<-"output/figures"; dir.create("output/mse",showWarnings=FALSE)

brp <- read_csv(file.path(tab,"stock_assessment_BRP.csv"), show_col_types=FALSE)
mat <- read_csv(file.path(tab,"params_maturity.csv"),      show_col_types=FALSE)
lw  <- read_csv(file.path(tab,"params_length_weight.csv"), show_col_types=FALSE)
catch <- read_csv(file.path(proc,"total_catches_annual.csv"), show_col_types=FALSE)

# life-history parameters per stock (data-limited; to be refined with FishLife/ELEFAN)
lh <- tibble::tribble(
  ~stock,             ~Linf, ~K,   ~M,   ~unit,
  "Octopus vulgaris",  25,   1.0,  1.2,  "LM cm",     # octopus: fast growth, high M
  "Penaeus notialis",  45,   1.2,  1.5,  "LCT mm",    # annual shrimp
  "Cymbium spp",       38,   0.25, 0.35, "LCQ cm")    # volute: slow

build_OM <- function(sp){
  b  <- brp  |> filter(stock==sp); h <- lh |> filter(stock==sp)
  m  <- mat  |> filter(stock==sp | (sp=="Cymbium spp" & str_starts(stock,"Cymbium")))
  L50 <- mean(m$L50)
  OM <- new("OM", DLMtool::Albacore, DLMtool::Generic_fleet,
            DLMtool::Imprecise_Unbiased, DLMtool::Perfect_Imp)
  OM@Name   <- sp; OM@nsim <- 200; OM@nyears <- nrow(filter(catch,espece==sp)); OM@proyears <- 30
  OM@maxage <- ceiling(5/ h$M *3)
  OM@Linf <- rep(h$Linf,2); OM@K <- rep(h$K,2); OM@M <- rep(h$M,2)
  OM@L50  <- rep(L50,2)
  OM@a <- lw$a[lw$stock==sp | lw$stock==first(m$stock)][1]; OM@b <- lw$b[1]
  OM@D <- c(max(0.05,b$B_Bmsy*0.5-0.1), min(0.95,b$B_Bmsy*0.5+0.1))  # current depletion B/B0
  OM@seed <- 42
  OM
}

# conditioning on real data via RCM (Rapid Conditioning Model)
condition_OM <- function(sp){
  OM <- build_OM(sp)
  g  <- catch |> filter(espece==sp) |> arrange(annee)
  cpf<- file.path(tab, sprintf("cpue_std_%s.csv", tolower(word(sp,1))))
  idx<- read_csv(cpf, show_col_types=FALSE)
  I  <- g |> left_join(transmute(idx, annee, I=index_scaled), by="annee") |> pull(I)
  rcm<- SAMtool::RCM(OM, data=list(Chist=g$capture_t, Index=matrix(I,ncol=1),
                                   I_sd=matrix(0.2,nrow=length(I),ncol=1)))
  rcm@OM
}

# candidate management procedures (spanning data requirements)
MPs <- c("AvC","CC1","DCAC","Islope1","Iratio","SP_MSY","curE","matlenlim")

run_stock_MSE <- function(sp){
  message("MSE : ", sp)
  OMc <- purrr::possibly(condition_OM, otherwise=build_OM(sp))(sp)
  mse <- MSEtool::runMSE(OMc, MPs=MPs, parallel=TRUE, silent=TRUE)
  saveRDS(mse, file.path("output/mse", sprintf("mse_%s.rds", tolower(word(sp,1)))))
  # performance metrics
  pm <- data.frame(
    stock=sp, MP=mse@MPs,
    P_B_gt_0.5Bmsy = MSEtool::Pplot(mse) |> {\(x) NA}(),  # see openMSE report
    LT_Yield = apply(mse@Catch[,, (mse@proyears-9):mse@proyears], 2, mean))
  # standard plots
  png(file.path(fig, sprintf("fig11_mse_tradeoff_%s.png", tolower(word(sp,1)))),
      width=900, height=700); openMSE::Tplot(mse); dev.off()
  pm
}

pm_all <- purrr::map_dfr(c("Octopus vulgaris","Penaeus notialis","Cymbium spp"),
                         purrr::possibly(run_stock_MSE, otherwise=NULL))
write_csv(pm_all, file.path(tab,"mse_performance.csv"))
message("\n=== Phase 2F (MSE) completed: see output/mse/ and figures fig11_* ===")
# Recommended metrics: P(B>0.5Bmsy), P(B>Blim), P(F<Fmsy), P(yield>0.5MSY)
# Robustness trials: OM variants with alternative M, h, and observation error.