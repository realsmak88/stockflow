# Generation of fictitious length frequency data
set.seed(123)

n <- 6000

lengths <-
  
  rlnorm(
    
    n,
    
    meanlog=log(21),
    
    sdlog=0.23
    
  )

lengths <- round(lengths,1)

Month <- sample(
  
  1:12,
  
  n,
  
  replace=TRUE,
  
  prob=c(
    
    0.05,
    
    0.05,
    
    0.08,
    
    0.09,
    
    0.11,
    
    0.13,
    
    0.12,
    
    0.10,
    
    0.09,
    
    0.07,
    
    0.06,
    
    0.05)
  
)

Year <- sample(
  
  2020:2024,
  
  n,
  
  replace=TRUE
  
)

FishID <- 1:n

df <-
  
  data.frame(
    
    FishID,
    
    Year,
    
    Month,
    
    Length=lengths
    
  )

write.csv(
  
  df,
  
  "data/raw/lengths.csv",
  
  row.names=FALSE
  
)

# Generation of fictitious maturity data
Length <- seq(10,35,0.5)

Pmat <- plogis(
  
  Length,
  
  18,
  
  1.8
  
)

write.csv(
  
  data.frame(
    
    Length,
    
    Maturity=Pmat
    
  ),
  
  "data/raw/maturity.csv",
  
  row.names=FALSE
)

# Generation of fictitious catch data
Year <- 2020:2024

Catch <- c(
  
  1200,
  
  1180,
  
  1155,
  
  1080,
  
  1015
  
)

write.csv(
  
  data.frame(
    
    Year,
    
    Catch
    
  ),
  
  "data/raw/catches.csv",
  
  row.names=FALSE
)

# Generation of fictitious effort data
Year <- 2020:2024

Effort <- c(
  
  420,
  
  430,
  
  450,
  
  470,
  
  495
  
)

write.csv(
  
  data.frame(
    
    Year,
    
    Effort
    
  ),
  
  "data/raw/effort.csv",
  
  row.names=FALSE
)

library(TropFishR)

source("R/import.R")
source("R/validation.R")
source("R/prepare/prepare_lengths.R")
source("R/prepare/prepare_lfq.R")
source("R/prepare/prepare_tropfish.R")
source("R/growth.R")


lfq <- prepare_tropfish(
  stock = lengths,
  bin_size = 1
)

class(lfq)

str(lfq)

names(lfq)

plot(lfq, Fname="catch", hist.sc = 1)

lfq <- TropFishR::lfqRestructure(
  lfq,
  MA = 7
)

plot(lfq, hist.sc = 0.75)

names(lfq)

range(lengths$Length)

lfq$midLengths

pw <- TropFishR::powell_wetherall(
  lfq,
  catch_columns = 1:12
)

sa <- TropFishR::ELEFAN_SA(
  lfq,
  seasonalised = FALSE
)

sa$par

ga <- TropFishR::ELEFAN_GA(
  lfq,
  seasonalised = FALSE
)

ga$par


el <- run_elefan(
    lfq,
    Linf_range=seq(60,100, lenght.out = 100),
    K_range=seq(0.01,0.3, lenght.out = 100)
)

ga <- run_elefan_ga(lfq)

sa <- run_elefan_sa(lfq)

extract_growth_parameters(el)
extract_growth_parameters(ga)
extract_growth_parameters(sa)

comparison <- compare_growth_models(
  el,
  ga,
  sa
)

comparison <- compare_growth_models(
  el,
  ga,
  sa,
  Linf_range=seq(60,100, length.out = 100),
  K_range=seq(0.01,3, length.out = 100)
)

comparison


print_growth_comparison(comparison)

best_growth_model(
  comparison,
  criterion="Rn_max"
)

best_growth_model(
  comparison,
  criterion="balanced"
)

best_growth_model(
  comparison,
  criterion="score"
)

best <- best_growth_model(
  comparison,
  criterion="Rn_max"
)

best$ranking

growth_params <- extract_growth_parameters(sa)

Linf <- growth_params$Linf
K    <- growth_params$K

cc <- TropFishR::catchCurve(
  
  lfq,
  
  catch_columns = 1:length(lfq$midLengths),
  
  calc_ogive = TRUE,
  
  plot = TRUE,
  
  Linf = sa$result$par$Linf,
  
  K = sa$result$par$K
  
)

cc <- estimate_catchcurve(
  lfq,
  ga,
  catch_columns =1:12
)

source("R/mortality.R")

mortality <- run_mortality(
  lfq,
  ga,
  temp=25
)