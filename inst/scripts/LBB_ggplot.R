# Length-based Bayesian Biomass estimator (LBB)
# Fits LBB model to length frequency data to estimate Linf, Lc, M/K, F/K 
# Derives reference points F/M, Z/K, Lopt, Lc_opt, B/B0, B/Bmsy, Y/R
# Main code developed by Rainer Froese in May-June 2017, modified in April-May 2018
# Gianpaolo Coro and Henning Winker did the JAGS coding
# Gianpaolo added the code for "Best year" in April 2019
# Deng Palomares indicated common errors and alert messages, after experience in courses, in October 2019
# Gives option in the ID file to correct for the piling-up effect, with Pile=0 no correction, Pile=1 full correction, Pile=999 degree of correction determined by fit
# Kamarel Ba modified the code to use ggplot2 for plotting in April 2026, and added error management for the input files in June 2024

# Automatic package installation
list.of.packages <- c("R2jags", "Hmisc","lattice","survival","Formula","ggplot2","crayon")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

rm(list=ls(all=TRUE)) # clear previous variables etc
options(digits=3) # displays all numbers with three significant digits as default
graphics.off() # close graphics windows from previous sessions
library(R2jags)
library(Hmisc)
library(crayon) # to display bold and italics in console

# Select stock to be analysed
Stock       <-  "Sardinella aurita COPACE" # "Bonga_West_Africa"  #"DPS_GSA22" # "Ench_cim22-24"  # "Myox_scor_22-24"

# Select file with stock ID info
ID.File     <-  "Aurita_ID.csv" # "Example_ID.csv"    

# Settings
n.sim       <- 10   # ifelse(Stock %in% c("CodRedFSim"),1,10) # number of years to be created in simulations
smooth.ts   <- T    # use three years moving average for B/B0 time series

##############################################################
#  Functions
##############################################################
#--------------------------------------------------------
# Exploited B/B0 ratio from B&H equations, for variable F
#--------------------------------------------------------
# assuming that reported lengths are the lower bounds of length classes
# get lowest exploited (>= 0.01 F) length class and class width

BH <- function(AllLength,Linf,MK,FK,GausSel,selpar1,selpar2) {
  if(GausSel==F) {
    r.Lc     <- selpar1
    r.alpha  <- selpar2 
    Lx       <- AllLength[AllLength >= Linf*(r.Lc-4.59/r.alpha)][1]
  } else if(GausSel==T) {
    r.GLmean <- selpar1
    r.SD     <- selpar2
    Lx       <- AllLength[AllLength >= Linf*(r.GLmean-3*r.SD)][1]
  }
  class.width  <- median(diff(sort(unique(AllLength))))
  FM <- FK/MK
  
  # Linf=120;Lx=22.5;r.Lc=0.2917;r.alpha=60;MK=1.5385;FK=0.7692;FM=0.5;ZK=2.3077 
  # uncomment above row for comparison of Y'R= 0.0332, B/B0=0.467 with CodLightSim
  r            <- vector() # auxilliary reduction factor
  G            <- vector() # product of reduction factors
  SL.bh        <- vector() # selection at length
  YR1.2        <- vector() # relative yield per recruit per length class
  CPUER1.2     <- vector() # relative CPUE per recruit per length class
  B1.2         <- vector() # relative unexploited biomass per recruit by length class
  L.bh         <- seq(from=Lx, to=Linf, by=class.width) # lengths to be considered
  r.L.bh       <-  L.bh / Linf # standardized lengths
  
  # calculate selection, Y'/R and CPUE'/R for every length class
  for(o in 1 : length(r.L.bh)) { 
    if(GausSel==F) {
      if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc))), # mean selection in length class
                                                1/(1+exp(-r.alpha*(r.L.bh[o+1]-r.Lc)))))
      } else SL.bh[o] <- 1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc)))
    } else if(GausSel==T) { # gill net selection 
      if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2))), # mean selection in length class
                                                exp(-((r.L.bh[o+1]-r.GLmean)^2/(2*r.SD^2)))))
      } else SL.bh[o] <- exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2)))
    } # end of calculation of selectivity loop
    
    if(o<length(r.L.bh)) {
      r[o]       <- (1-r.L.bh[o+1])^(FK*SL.bh[o])/(1-r.L.bh[o])^(FK*SL.bh[o]) 
      G[o]       <- prod(r[1:o]) }
    if(o==1) {
      YR1.2[o] <-(FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                     (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                  (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) -
        (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                                                                                (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                           (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o] 
    } else if(o==length(r.L.bh)) {
      YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                      (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                   (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] 
    } else {
      YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                      (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                   (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] -
        (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                                                                                (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                           (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o]              
    } # end of loop to calculate yield per length class
    
    CPUER1.2[o] <- YR1.2[o] / FM # CPUE/R = Y/R divided by F/M
    
    if(o<length(r.L.bh)) {
      B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                      (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)) -
                    (1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/MK)+3*(1-r.L.bh[o+1])^2/
                                          (1+2/MK)-(1-r.L.bh[o+1])^3/(1+3/MK)))*SL.bh[o]
    } else {
      B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                      (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)))*SL.bh[o]
    }
  } # end of B&H loop through length classes
  BB0   <- sum(CPUER1.2)/sum(B1.2)
  YR    <- sum(YR1.2)
  if (!is.na(BB0) && BB0 < 0.25) {
    YR <- YR * BB0 / 0.25
  }
  #if(BB0 < 0.25) YR <- YR * BB0 / 0.25 # reduce YR if recruitment and thus productivity is reduced
  return(list(BB0,YR))
  
} # end of BH function

#------------------------------------------------------------
# Function to aggregate data by year
#------------------------------------------------------------
AG  <- function(dat) { # where dat contains dat$Year, dat$Length in cm, dat$CatchNo
  
  # aggregate normalized annual LFs by weighing with square root of sample size
  # get sum of frequencies per year
  sum.Ny  <- aggregate(Freq~Year,dat,sum)$Freq  
  # get the sqrt of the sum of frequencies for every year
  sqrt.Ny <- sqrt(sum.Ny) 
  # get highest frequency in each year
  max.Ny <- aggregate(Freq~Year,dat,max)$Freq
  # get Number of Length bins in each year
  binsN <- aggregate(Freq~Year,dat,length)$Freq    
  # create vectors for sqrt.Ni and sum.Ni to weigh LF data
  sqrt.Ni = rep(sqrt.Ny,binsN)
  sum.Ni = rep(sum.Ny,binsN)
  #Do weighing
  # Divide all years by sum.Ni and multiply by sqrt.Ni
  LF.w = dat$Freq/sum.Ni*sqrt.Ni  
  # Aggregate
  LF = aggregate(LF.w, by=list(dat$Length),FUN=sum)
  # Add correct column names
  colnames(LF) <- c("Length","Freq")         
  return(LF)
} #end of aggregate function

#-----------------------------------------------------------
# Function to plot LBB-fit for a single year
#-----------------------------------------------------------
# expects lengths relative to Linf (L/Linf)
library(ggplot2)

plot.year <- function(r.L.y, r.Freq.y, r.Lopt, r.Freq.pred.y, SL1, SL2, MK, FK, Linf, main) {
  
  # 1. Preparing the data in a data.frame
  df_plot <- data.frame(
    L = r.L.y,
    Freq_Obs = r.Freq.y,
    Freq_Pred = r.Freq.pred.y
  )
  
  # Calculation of the max for positioning the texts
  max_freq <- max(r.Freq.y, na.rm = TRUE)
  
  # 2. Creating the ggplot chart
  p <- ggplot(df_plot, aes(x = L)) +
    # Observed points
   # geom_point(aes(y = Freq_Obs), size = 1, alpha = 0.6) +
    geom_point(aes(y = Freq_Obs), shape = 21, size = 2) +
    # Prediction line (the LBB model)
    geom_line(aes(y = Freq_Pred), color = "red", linewidth = 0.8) +
    # Vertical lines for Linf and Lopt
    geom_vline(xintercept = 1, color = "darkgreen", linetype = "solid") +
    geom_vline(xintercept = r.Lopt, color = "darkgreen", linetype = "solid") +
    # Text annotations on the chart
    annotate("text", x = 1, y = 1.15 * max_freq, label = "Linf", color = "darkgreen") +
    annotate("text", x = r.Lopt, y = 1.15 * max_freq, label = "Lopt", color = "darkgreen") +
    # Parameter box (Linf and Z/K)
    annotate("label", x = 0.15, y = 0.8 * max_freq, 
             label = paste0("Linf = ", format(Linf, digits = 3)), fill = "white", label.size = NA) +
    annotate("label", x = 0.15, y = 0.6 * max_freq, 
             label = paste0("Z/K = ", format(MK + FK, digits = 3)), fill = "white", label.size = NA) +
    # Formatting
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1.2 * max_freq)) +
    labs(title = main,
         x = "Length / Linf",
         y = "Relative Frequency") +
    theme_minimal() +
    theme(panel.border = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  # Displaying the chart
  print(p)
}

#-------------------------------------------------------
# Function to apply preceding 3-years moving average
#-------------------------------------------------------
ma    <- function(x){
  x.1    <-   stats::filter(x,rep(1/3,3),sides=1)
  x.1[1] <- x[1]
  x.1[2] <- (x[1]+x[2])/2
  return(x.1)
}


#############################################################
# read files with ID and with LF data to be analyzed
#############################################################
# read ID data
tryCatch({
  dat.ID         <- read.csv(ID.File, header=T,sep = "," ,stringsAsFactors=F) 
},
error=function(cond) {
  cat("ERROR: Bad structure of input CSV file - hints to check file consistency:\nCheck your CSV file by displaying it as a text file,\ni.e. right click on the file and click Open with 'Notepad' or any other text file displayer that you have on your laptop,\nCheck that there are no floating commas at the end of each line. \nIf there are floating commas, then delete those and assure that the data being saved has not been corrupted, i.e. the columns did not move between rows and that all of the data is intact.\nGo to the last line and press the carriage return if a final blank line is not present \nErase any floating commas, then rerun the software.\n")
  stop()
}
)

#ERRORS MANAGEMENT by Deng
if (dim(dat.ID)[1]==1 && dim(dat.ID)[2]==1 && regexpr(";", as.character(dat.ID))[[1]]>10){
  cat("ERROR: The CSV file is using ';' instead of ',': To solve this, go to File, Options (in Windows) and advanced settings and change the list delimiter from semi colon to a comma. In Mac, close Excel, click on Apple icon, select Language and Region, then Advanced, then change the Number separators Grouping from semi-colon to comma then press OK.\n")
  stop()
}


# restrict ID data to selected Stock
dat.ID         <- dat.ID[dat.ID$Stock==Stock,]

if (length(dat.ID$File)>=2){
  cat("ERROR: Duplicate entry for stock",Stock,". Please use different identifiers for different entries (i.e. lines in the ID file).\n",sep="")
  stop()
}

if (is.na(dat.ID$mm.user)){
  cat("ERROR: mm.user is NA while it should be TRUE or FALSE.\n")
  stop()
}   

if (is.na(dat.ID$GausSel)){
  cat("ERROR: GausSel is NA while it should be TRUE or FALSE.\n")
  stop()
}

if (is.na(dat.ID$MergeLF)){
  cat("ERROR: MergeLF is NA while it should be TRUE or FALSE.\n")
  stop()
}

if (is.na(dat.ID$Pile) || (dat.ID$Pile!=1 && dat.ID$Pile!=0 && dat.ID$Pile!=999)){
  cat("ERROR: The Pile column in the ID file cannot be left blank or with NA but which should have either of three values: Pile=0 no correction, Pile=1 full correction, Pile=999 degree of correction determined by fit (see user guide). Also, watch out for possible blank lines in the ID file.\n")
  stop()
}

if (dat.ID$mm.user==F){
  cat("REMINDER: Lengths in the ID file should be reported in cm, whereas lengths in the catch-at-length file should be always reported in mm.\n",sep="")
}else{
  cat("REMINDER: Lengths in the ID file should be reported in mm. Lengths in the catch-at-length file should be reported in mm.\n",sep="")
}

if (!file.exists(dat.ID$File)){
  cat("ERROR: Filename in File column does not correspond to filename of raw data. For example, check if the filename lacks the '.csv' extension of the catch file is misspelled or has wrong case.\n")
  stop()
}


# read LF data
dat.raw        <- read.csv(dat.ID$File, header=T, stringsAsFactors=F) 

# restrict LF data to selected stock
dat.raw        <- dat.raw[dat.raw$Stock == Stock,] 

if (dim(dat.raw)[1]==0){
  cat("ERROR: Stock ID in ID file does not correspond to the Stock ID in DAT or Catch file (misspelled/wrong case).\n")
  stop()
}

# remove NA records
dat.raw    <- dat.raw[which(is.na(dat.raw$CatchNo)==F),]

# restrict analysis to one or more gears
if(is.na(dat.ID$Gears.user[1])==FALSE) dat.raw <- dat.raw[dat.raw$Gear %in% dat.ID$Gears.user,]

# make sure data are numeric
dat.raw$Length  <- as.numeric(dat.raw$Length)
dat.raw$CatchNo <- as.numeric(dat.raw$CatchNo)
dat.raw$Year    <- as.integer(dat.raw$Year)

# if StartYear is given, restrict data to >= StartYear
if(is.na(dat.ID$StartYear)==F) dat.raw <- dat.raw[dat.raw$Year>=dat.ID$StartYear,]

# if EndYear is given, restrict data to <= EndYear
if(is.na(dat.ID$EndYear)==F) dat.raw <- dat.raw[dat.raw$Year<=dat.ID$EndYear,]

# if Years.user are given, restrict data to these years
if(is.na(dat.ID$Years.user[[1]])==F) dat.raw <- dat.raw[dat.raw$Year %in% (strsplit(dat.ID$Years.user, ","))[[1]],] # code from GP

# use largest fish as Lmax
Lmax       <- max(dat.raw$Length)/10
# use median of largest fish per year as Lmax.med
Lmax.med   <- median(as.numeric(by(dat.raw$Length[dat.raw$CatchNo>0],dat.raw$Year[dat.raw$CatchNo>0],max)))/10

# if Linf.user is given, restict data to < Linf.user
if(is.na(dat.ID$Linf.user)==F) dat.raw <- dat.raw[dat.raw$Length<(dat.ID$Linf.user*ifelse(dat.ID$mm.user==TRUE,1,10)),] 

# if Lcut.user is given, restrict data to >= Lcut.user
if(is.na(dat.ID$Lcut.user)==F) dat.raw <- dat.raw[dat.raw$Length>=(dat.ID$Lcut.user*ifelse(dat.ID$mm.user==TRUE,1,10)),]

# sort data by year and length
dat.raw    <- dat.raw[order(dat.raw$Year,dat.raw$Length),]

# check for selected year to show B/B0
if(length(dat.ID$Year.select[dat.ID$Stock==Stock]) != 0 && is.na(dat.ID$Year.select[dat.ID$Stock==Stock])==F) {
  Year.sel   <- dat.ID$Year.select[dat.ID$Stock==Stock]
} else {Year.sel <- NA}

# Put data into vectors
StartYear  <- min(dat.raw$Year)
EndYear    <- max(dat.raw$Year)
AllYear    <- dat.raw$Year
AllLength  <- dat.raw$Length
if(dat.ID$mm.user==FALSE) AllLength <- AllLength/10 
AllFreq    <- dat.raw$CatchNo 
Years      <- sort(unique(AllYear))
nYears     <- length(Years)

# if data are simulated, add noise and n.sim more years
if(substr(Stock,start=nchar(Stock)-2,stop=nchar(Stock))=="Sim") {  
  n.L.sim      <- length(AllLength)
  AllYearSim   <- AllYear
  AllLengthSim <- AllLength
  AllFreqSim   <- rlnorm(n=n.L.sim,mean=log(AllFreq),sd=0.1)
  if(!(Stock %in% c("CodfFSim","CodRecSim"))) {  # CodfFSim and CodRecSim are simulations that should run for only one year
    for(i in 1 : (n.sim-1)) {
      AllYearSim   <- append(AllYearSim,AllYear+i)
      AllLengthSim <- append(AllLengthSim,AllLength)
      AllFreqSim   <- append(AllFreqSim,rlnorm(n=n.L.sim,mean=log(AllFreq),sd=0.1))
    }
    AllYear    <- AllYearSim
    AllLength  <- AllLengthSim
    AllFreq    <- AllFreqSim
    Years      <- sort(unique(AllYear))
    nYears     <- length(Years)
    EndYear    <- Years[nYears] }
} # end of simulation loop

#-------------------------------------------------------------------------------
# 1. INITIALIZATION AND CALCULATIONS BY YEAR
#-------------------------------------------------------------------------------
output_dir <- "LBB_Results"
if(!dir.exists(output_dir)) dir.create(output_dir)

# We create the empty list here before filling it
all_data_list <- list()

for(w in 1:nYears) {
  # Annual extraction and cleaning
  df.p <- data.frame(
    Year   = AllYear[AllYear == Years[w] & AllFreq > 0],
    Length = AllLength[AllYear == Years[w] & AllFreq > 0],
    Freq   = AllFreq[AllYear == Years[w] & AllFreq > 0]
  )
  
  if(nrow(df.p) > 0) {
    # Aggregation via the package's AG function
    LF.p <- AG(dat = df.p) 
    
    # We add the year so that ggplot can facet later
    LF.p$Year <- Years[w]
    
    # We store it in the list
    all_data_list[[w]] <- LF.p
  }
}

#-------------------------------------------------------------------------------
# 2. PREPARATION AND CHARTS (Now all_data_list exists!)
#-------------------------------------------------------------------------------
df_all_years <- do.call(rbind, all_data_list)

# Check: if df_all_years is empty, we stop to avoid a ggplot error
if(is.null(df_all_years)) {
  stop("Error: No data could be aggregated. Check your CSV files.")
}

# Creation of the faceted chart (page by page)
n_pages <- ceiling(nYears / 6)

for(z in 1:n_pages) {
  years_selection <- Years[((z-1)*6 + 1) : min(z*6, nYears)]
  df_sub <- df_all_years[df_all_years$Year %in% years_selection, ]
  
  p <- ggplot(df_sub, aes(x = Length, y = Freq)) +
    geom_point(size = 2, shape = 21, color = "steelblue", alpha = 0.8) +
    geom_line(color = "steelblue", alpha = 0.6) + 
    facet_wrap(~Year, scales = "free_y", ncol = 3, nrow = 2) +
    coord_cartesian(xlim = c(0, Lmax * 1.1)) + 
    labs(title = paste("Size frequency check -", Stock),
         subtitle = paste("Page", z, "of", n_pages),
         x = paste0("Length ", ifelse(dat.ID$mm.user, "(mm)", "(cm)")),
         y = "Frequency") +
    theme_minimal()
  
  print(p)
  ggsave(filename = file.path(output_dir, paste0("LF_Check_Page_", z, ".png")), 
         plot = p, width = 12, height = 8, dpi = 300)
  
  g <- ggplot(df_sub, aes(x = Length, y = Freq)) +
    geom_bar(stat = "identity", fill = "gray80", color = "gray60", alpha = 0.5) +
    #geom_line(color = "steelblue", alpha = 0.6) + 
    facet_wrap(~Year, scales = "free_y", ncol = 3, nrow = 2) +
    coord_cartesian(xlim = c(0, Lmax * 1.1)) + 
    labs(title = paste("Size frequency check -", Stock),
         subtitle = paste("Page", z, "of", n_pages),
         x = paste0("Length ", ifelse(dat.ID$mm.user, "(mm)", "(cm)")),
         y = "Frequency") +
    theme_minimal()
  
  print(g)
  ggsave(filename = file.path(output_dir, paste0("LF_Check_Page_", z, ".png")), 
         plot = g, width = 12, height = 8, dpi = 300)
}

#--------------------------------------------------
# Print warning if MergeLF is used
#--------------------------------------------------
if(dat.ID[dat.ID$Stock==Stock]$MergeLF==TRUE) {
  cat("Attention: LFs in subsequent years are merged and the first year is identical with the second")
}

# -------------------------------------------------
# Print years and Lmax across all data for early orientation
#--------------------------------------------------
cat("\n Lmax =",Lmax,", median Lmax =",Lmax.med,"cm, for potential setting of Linf.user in ID file \n\n") 
cat(" Years in data set (for potential cut & paste into Years.user in ID file):\n", paste(Years,collapse=","),"\n")
cat("If error without hint occurs, copy years into Years.user and delete next year to be processed from string\n\n")

#----------------------------------------------------
# Create matrix to store annual estimates
#----------------------------------------------------
Ldat      <- data.frame(Stock=rep(Stock,nYears),Year=rep(NA,nYears),
                        Linf=rep(NA,nYears),
                        Linf.lcl=rep(NA,nYears),
                        Linf.ucl=rep(NA,nYears),
                        Lc=rep(NA,nYears), # for trawl selection
                        Lc.lcl=rep(NA,nYears),
                        Lc.ucl=rep(NA,nYears),
                        Lmean=rep(NA,nYears),
                        r.alpha=rep(NA,nYears),
                        r.alpha.lcl=rep(NA,nYears),
                        r.alpha.ucl=rep(NA,nYears),
                        r.GLmean=rep(NA,nYears),r.SD=rep(NA,nYears), # for gill net selection
                        MK=rep(NA,nYears),
                        MK.lcl=rep(NA,nYears),
                        MK.ucl=rep(NA,nYears),
                        FK=rep(NA,nYears),
                        FK.lcl=rep(NA,nYears),
                        FK.ucl=rep(NA,nYears),
                        ZK=rep(NA,nYears),
                        ZK.lcl=rep(NA,nYears),
                        ZK.ucl=rep(NA,nYears),
                        FM=rep(NA,nYears),
                        FM.lcl=rep(NA,nYears),
                        FM.ucl=rep(NA,nYears),
                        r.Lopt=rep(NA,nYears),
                        BB0=rep(NA,nYears),
                        BB0.lcl=rep(NA,nYears),
                        BB0.ucl=rep(NA,nYears),
                        YR=rep(NA,nYears),
                        YR.lcl=rep(NA,nYears),
                        YR.ucl=rep(NA,nYears),
                        perc.mat=rep(NA,nYears),
                        L95=rep(NA,nYears))

Lfit       <- matrix(list(),nYears,3)           

#--------------------------------------------------------------------------------------
# Use aggregated LF data for estimation of Linf (and overall Z/K)
#--------------------------------------------------------------------------------------
df        <- data.frame(AllYear,AllLength,AllFreq)
names(df) <- c("Year","Length","Freq")

LF.all    <- AG(dat=df) # function to aggregate data by year

# standardize to max Freq
LF.all$Freq = LF.all$Freq/max(LF.all$Freq) 
# remove leading empty records
LF.all     <- LF.all[which(LF.all$Freq>0)[1] : length(LF.all$Length),]
# remove trailing empty records
LF.all     <- LF.all[1 : which(LF.all$Length==max(LF.all$Length[LF.all$Freq>0])),]

# get number of records in LF.all
n.LF.all   <- length(LF.all$Length) 

# If no Linf is provided by the user (preferred), determine Linf from fully selected LF:
# Freq=Nstart*exp(ZK*(log(1-L/Linf)-log(1-Lstart/Linf)))
# Nstart is canceled out when dividing both sides by their sums
# ---------------------------------------------------------
# determine start values of selection ogive to find first fully selected length class Lstart
L10         <- LF.all$Length[which(LF.all$Freq>0.1)[1]] # use length at 10% of peak frequency as proxy for L10
L90         <- LF.all$Length[which(LF.all$Freq>0.9)[1]] # use length at 90% of peak frequency as proxy for L90
Lc.st       <- ifelse(is.na(dat.ID$Lc.user)==TRUE,(L10 + L90)/2,dat.ID$Lc.user)  # use mean of L10 and L90 as proxy for Lc, else user input
alpha.st    <- -log(1/LF.all$Freq[which(LF.all$Freq>0.1)[1]])/(L10-Lc.st) # use rearranged logistic curve to estimate slope alpha

# determine start values for Linf and Z/K 
Linf.st     <- ifelse(is.na(dat.ID$Linf.user)==F,dat.ID$Linf.user,Lmax.med) # use Linf.user or median Lmax across years as start value for Linf in nls analysis
Lmean.st    <- sum(LF.all$Length[LF.all$Length>=Lc.st]*LF.all$Freq[LF.all$Length>=Lc.st])/
  sum(LF.all$Freq[LF.all$Length>=Lc.st])
MK.st       <- ifelse(is.na(dat.ID$MK.user)==TRUE, 1.5,dat.ID$MK.user) # default 1.5
ZK.st       <- (Linf.st-Lmean.st)/(Lmean.st-Lc.st)       # the Holt equation
FK.st       <- ifelse((ZK.st-MK.st)>0,ZK.st-MK.st,0.3)   # prevent M/K being larger than Z/K

# get vectors with fully selected length classes for Linf estimation
if(is.na(dat.ID$Lstart.user)==FALSE) {Lstart <- dat.ID$Lstart.user} else {
  Lstart     <- (alpha.st*Lc.st-log(1/0.95-1))/alpha.st   # Length where selection probability is 0.95  
  # test if there are enough (>=4) length classes for estimation of aggregated Linf and ZK 
  Lstart.i   <- which(LF.all>=Lstart)[1]
  Lmax.i     <- length(LF.all$Length)
  peak.i     <- which.max(LF.all$Freq)
  if(Lstart.i<(peak.i+1)) Lstart <- LF.all$Length[peak.i+1] # make sure fully selected length starts after peak 
  if((Lmax.i-Lstart.i)<4) Lstart <- LF.all$Length[Lstart.i-1] # make sure enough length classes are available
}
# do not include Lmax to allow Linf < Lmax and to avoid error in nls when Linf-L becomes negative
L.L         <- LF.all$Length[LF.all$Length >= Lstart  & LF.all$Length < Linf.st]
L.Freq      <- LF.all$Freq[LF.all$Length>=L.L[1]& LF.all$Length < Linf.st]

if(length(L.L)<4) {
  #modification by Gianpaolo 09 07 17 
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(6,4)
  } else if(grepl("linux",tolower(Sys.info()['sysname']))) {X11(6,4)
  } else {quartz(6,4)}
  
  plot(x=LF.all$Length,y=LF.all$Freq, bty="l",main=Stock)
  lines(x=c(Lstart,Lstart),y=c(0,0.9*max(LF.all$Freq)),lty="dashed")
  text(x=Lstart,y=max(LF.all$Freq),"Lstart")
  lines(x=c(Linf.st,Linf.st),y=c(0,0.9*max(LF.all$Freq)),lty="dashed")
  text(x=Linf.st,y=max(LF.all$Freq),"Lmax")
  stop("Too few fully selected data points: set Lstart.user\n")}

# standardize frequencies by dividing by sum of observed frequencies, needed to drop NLstart from equation
sum.L.Freq  <- sum(L.Freq)
L.Freq      <- L.Freq/sum.L.Freq

# use nls() to find Linf-ZK combination with least residuals
if(is.na(dat.ID$Linf.user)==TRUE) {
  Linf.mod    <- nls(L.Freq ~ ((Linf-L.L)/(Linf-Lstart))^ZK /
                       sum(((Linf-L.L)/(Linf-Lstart))^ZK),
                     start=list(ZK=ZK.st,Linf=Linf.st),
                     lower=c(0.5*ZK.st,0.999*Linf.st), 
                     upper=c(1.5*ZK.st,1.2*Linf.st), 
                     algorithm = "port")
  
  ZK.nls       <- as.numeric(coef(Linf.mod)[1])
  ZK.nls.sd    <- as.numeric(coef(summary(Linf.mod))[,2][1])
  ZK.nls.lcl   <- ZK.nls-1.96*ZK.nls.sd
  ZK.nls.ucl   <- ZK.nls+1.96*ZK.nls.sd
  Linf.nls     <- as.numeric(coef(Linf.mod)[2])
  Linf.nls.sd  <- as.numeric(coef(summary(Linf.mod))[,2][2])
  Linf.lcl     <- Linf.nls-1.96*Linf.nls.sd
  Linf.ucl     <- Linf.nls+1.96*Linf.nls.sd
} else {  # end of loop to determine Linf and ZK.L
  # use given Linf and determine ZK.L
  # use Linf provided by user if given
  Linf.nls    <- dat.ID$Linf.user 
  Linf.nls.sd <- 0.01*dat.ID$Linf.user
  ZK.mod      <- nls(L.Freq ~ exp(ZK*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls)))/
                       sum(exp(ZK*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls)))),
                     start=list(ZK=ZK.st),
                     lower=c(0.7*ZK.st),
                     upper=c(1.3*ZK.st), 
                     algorithm = "port")
  ZK.nls       <- as.numeric(coef(ZK.mod)[1])
  ZK.nls.sd    <- as.numeric(coef(summary(ZK.mod))[,2][1])
  ZK.nls.lcl   <- ZK.nls-1.96*ZK.nls.sd
  ZK.nls.ucl   <- ZK.nls+1.96*ZK.nls.sd
  
} # end of loop if Linf is given by user



# get vector of all lengths <= prior Linf to avoid error in equation
AllFreq       <- AllFreq[AllLength <= Linf.nls]
AllYear       <- AllYear[AllLength <= Linf.nls]
AllLength     <- AllLength[AllLength <= Linf.nls]


#-----------------------------------------
# Start LF analysis by year
#-----------------------------------------
cat("Running   Jags model to fit SL and N distributions for",dat.ID$Species,"\n")
#jagsFit<-c() #modification by GP to select the best year
# To be placed BEFORE your loop over the years
jagsFit <- rep(NA, nYears)

#-------------------------------------------------------------------------------
# 1. INITIALIZATION (To be placed before the loop for(Year in Years))
#-------------------------------------------------------------------------------
library(ggplot2)
library(gridExtra)
library(crayon)

# List to store the annual plots
all_fits_plots <- list()
output_dir <- "LBB_Results"
if(!dir.exists(output_dir)) dir.create(output_dir)

#-------------------------------------------------------------------------------
# 2. MAIN LOOP (Modified to capture the plots)
#-------------------------------------------------------------------------------
i = 0 
for(Year in Years) {
  i = i + 1 # i is the index of Years, which may contain gaps 
  # if MergeLF==TRUE and if this is the second or heigher year and no simulation, aggregate LF with previous year LF
  if(dat.ID$MergeLF==TRUE & substr(Stock,start=nchar(Stock)-2,stop=nchar(Stock))!="Sim") {
    if(i==1) {AG.yr <- c(Year,Years[2])} else { # if first year, aggregate with second year
      AG.yr <- c(Years[i-1],Year) }
  } else AG.yr <- Year
  
  # aggregate data within the year (sometimes there are more than one sample per year)
  df        <- data.frame(AllYear[AllYear%in%AG.yr],AllLength[AllYear%in%AG.yr],AllFreq[AllYear%in%AG.yr])
  names(df) <- c("Year","Length","Freq")
  LF.y      <- AG(dat=df) # function to aggregate data by year and across years
  LF.y$Freq <- LF.y$Freq/sum(LF.y$Freq) # standardize frequencies
  
  # remove empty leading and trailing records
  LF.y        <- LF.y[which(LF.y$Freq>0)[1] : length(LF.y$Length),]
  LF.y        <- LF.y[1 : which.max(LF.y$Length[LF.y$Freq>0]),]
  # get vectors
  L.y         <- LF.y$Length
  r.Freq.y    <- LF.y$Freq
  
  # fill remaining zero frequencies with very small number, to avoid error
  r.Freq.y[r.Freq.y==0] <- min(r.Freq.y[r.Freq.y>0],na.rm=T)/100
  # enter data for this year into data frame
  Ldat$Year[i]     <- Year
  
  #-------------------------------------------------------------------------
  # Estimate annual parameters Lc, alpha, M/K, F/K from LF curve with trawl-type selection
  #-------------------------------------------------------------------------
  # determine priors 
  n.L         <- length(L.y)
  Linf.pr     <- Linf.nls
  Linf.sd.pr  <- ifelse(Linf.nls.sd/Linf.nls<0.01,Linf.nls.sd,0.01*Linf.nls) # restict prior CV of Linf to < 0.01
  MK.pr       <- MK.st
  MK.sd.pr    <- ifelse(is.na(dat.ID$MK.user)==TRUE,0.15,0.075)
  Pile        <- dat.ID$Pile
  
  if(dat.ID$GausSel==FALSE){ # apply trawl-like selection 
    Lc.pr        <- ifelse(is.na(dat.ID$Lc.user)==TRUE,1.02*Lc.st,dat.ID$Lc.user) # with 1.02 multiplier to account for systematic small underestimation
    Lc.sd.pr     <- ifelse(is.na(dat.ID$Lc.user)==TRUE,0.1*Lc.pr,0.05*Lc.pr) # assume narrower SD if Lc is given by user
    r.max.Freq   <- max(r.Freq.y,na.rm=T) 
    r.alpha.pr   <- -log(r.max.Freq/r.Freq.y[which(r.Freq.y>(0.1*r.max.Freq))[1]])/(L10/Linf.nls-Lc.st/Linf.nls) # relative alpha for standardized data
    r.alpha.sd.pr<- 0.025*r.alpha.pr 
    FK.pr        <- ifelse((ZK.nls-MK.st) > 0,ZK.nls-MK.st,0.3) # if Z/K <= M/K assume low F/K = 0.3 
    
    # list of data to pass to JAGS plus list of parameters to estimate   
    jags.data <- list ("r.Freq.y","L.y","n.L","Linf.pr","Linf.sd.pr","Lc.pr","Lc.sd.pr","r.alpha.pr","r.alpha.sd.pr","MK.pr","MK.sd.pr",
                       "FK.pr","Pile")
    jags.params <- c("r.alpha.d","Lc.d","SL","xN","FK.d","MK.d","Linf.d","pile.fac","Freq.pred")
    
    #---------------------------------
    # LBB JAGS model for trawl-like selection
    #---------------------------------
    sink("SLNMod.jags")
    cat("
  model {
  r.alpha.d_tau  <- pow(r.alpha.sd.pr, -2) 
  r.alpha.d      ~ dnorm(r.alpha.pr,r.alpha.d_tau) 

   Lc.d_tau  <- pow(Lc.sd.pr,-2)
   Lc.d      ~ dnorm(Lc.pr,Lc.d_tau) #       

   MK.d_tau  <-pow(MK.sd.pr, -2) # strong prior on M/K
   MK.d      ~ dnorm(MK.pr, MK.d_tau)

   Linf.tau  <- pow(Linf.sd.pr,-2) 
   Linf.d    ~ dnorm(Linf.pr,Linf.tau)
    
   FK.d       ~ dlnorm(log(FK.pr),4) # wide prior range for F/K

   SL[1]       ~ dlogis(0,1000)
   Freq.pred[1]<-0
   xN[1]       <-1

   p.low    <- ifelse(Pile==1,0.99,0)   
   p.hi     <- ifelse(Pile==0,0.01,1)
   pile.fac ~ dunif(p.low,p.hi)


   for(j in 2:n.L) {
    SL[j] <- 1/(1+exp(-r.alpha.d*(((L.y[j]+L.y[j-1])/2)/Linf.d-Lc.d/Linf.d))) # selection at mid-length of bin

    xN[j] <- xN[j-1]*((Linf.d-L.y[j])/(Linf.d-L.y[j-1]))^(MK.d+FK.d*SL[j]) # predicted numbers without pile-up
     
    cN[j] <- (xN[j-1]-xN[j])/(MK.d+FK.d*SL[j]) # predicted relative frequency with pile-up correction

    dN[j] <- cN[j]-xN[j] # difference between corrected and uncorrected frequencies
    
    uN[j] <- xN[j] + dN[j]*pile.fac # gradual application of correction with pile.fac between 0 and 1

		Freq.pred[j]<-uN[j]*SL[j] # relative frequencies of vulnerable individuals
		
    # normalize frequencies by dividing by sum of frequencies; multiply with 10 to avoid small numbers and with 1000 for effective sample size
    r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10*1000
  }	
  
  #><> LIKELIHOOD FUNCTION
  #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)
  r.Freq.y[2:n.L] ~ ddirch(r.Freq.pred[2:n.L])  
 
  } # END OF MODEL
    ",fill = TRUE)
    sink()
    
    MODEL = "SLNMod.jags"
    jagsfitSLN <- jags.parallel(data=jags.data, working.directory=NULL, inits=NULL, 
                                parameters.to.save=jags.params, 
                                model.file=paste(MODEL), 
                                n.burnin=15000, n.thin=10, n.iter=30000, n.chains=3)
    
    #jagsFit<-c(jagsFit,jagsfitSLN$BUGSoutput$pD) #modification by GP to select the best year according to the Deviance information criterion
    # Inside the loop (e.g., at index i)
    jagsFit[i] <- jagsfitSLN$BUGSoutput$mean$deviance
    
    # use median and percentiles
    Ldat$Lc[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$Lc.d)
    Ldat$Lc.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d,0.025)
    Ldat$Lc.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d,0.975)
    Ldat$Lmean[i]   <- sum(L.y[L.y>=Ldat$Lc[i]]*r.Freq.y[L.y>=Ldat$Lc[i]])/sum(r.Freq.y[L.y>=Ldat$Lc[i]])
    Ldat$r.alpha[i] <- median(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d)
    Ldat$r.alpha.lcl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d,0.025)
    Ldat$r.alpha.ucl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d,0.975)
    Ldat$MK[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
    Ldat$MK.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.025)
    Ldat$MK.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.975)
    Ldat$FK[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
    Ldat$FK.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.025)
    Ldat$FK.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.975)
    FMi             <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
    Ldat$FM[i]      <- median(FMi)
    Ldat$FM.lcl[i]  <- quantile(FMi,0.025)
    Ldat$FM.ucl[i]  <- quantile(FMi,0.975)
    ZKi             <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
    Ldat$ZK[i]      <- median(ZKi)
    Ldat$ZK.lcl[i]  <- quantile(ZKi,0.025)
    Ldat$ZK.ucl[i]  <- quantile(ZKi,0.975)
    Ldat$r.Lopt[i]  <- 3/(3+Ldat$MK[i])
    Ldat$Linf[i]    <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
    Ldat$Linf.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.025)
    Ldat$Linf.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.975)
    
  } # end of trawl-like selection 
  
  #----------------------------------------------------------------------
  # Estimate parameters GLmean, SD, F/K, M/K if selection is gillnet-like
  #---------------------------------------------------------------------- 
  if(dat.ID$GausSel==TRUE) {
    # determine priors
    # assume length at peak Freq as mean and distance to length at 80% of peak as SD of mean
    GLmean.st <- L.y[which.max(r.Freq.y)]
    # assume SD of Gaussian selection as distance between length at peak and length at 50% of peak
    Lc.pr     <- L.y[which(r.Freq.y >= (0.5*max(r.Freq.y)))][1]
    SD.st      <- max(GLmean.st-Lc.pr,0.25*GLmean.st)
    
    cat("Running Jags model to fit SL and N distributions for gillnet-like selection\n")
    
    n.L <- length(L.y)
    
    jags.data <- list ("n.L","GLmean.st","L.y","SD.st","ZK.nls","r.Freq.y","Linf.pr","Linf.sd.pr","MK.pr")
    jags.params <- c("GLmean.d","SD.d","SL","xN","FK.d","MK.d","Linf.d","Freq.pred")
    
    #---------------------------
    # JAGS model L-based with integral
    #---------------------------
    sink("SLNMod.jags")
    cat("
      model {
      GLmean.tau <- pow(0.1*GLmean.st,-2) 
      GLmean.d   ~ dnorm(GLmean.st,GLmean.tau)
      
      SD.tau    <- pow(0.2*SD.st,-2)
      SD.d      ~ dnorm(SD.st,SD.tau)
      
      MK.d_tau  <-pow(0.15,-2)
      MK.d      ~ dnorm(MK.pr,MK.d_tau)

      Linf.tau  <- pow(Linf.sd.pr,-2)
      Linf.d    ~ dnorm(Linf.pr,Linf.tau)
      
      FK        <- (ZK.nls-1.5) # ZK overestimated in gillnet selection, used as upper range
      FK.d      ~ dunif(0,FK)  

      SL[1]~ dlogis(0,1000)
      Freq.pred[1]<-0
      xN[1]<-1
      
      for(j in 2:n.L) {
        SL[j]<- exp(-((L.y[j]-GLmean.d)^2/(2*SD.d^2)))

        xN[j]<-xN[j-1]*exp((MK.d+FK.d*SL[j])*(log(1-L.y[j]/Linf.d)-log(1-L.y[j-1]/Linf.d)))
      
        cN[j] <- (xN[j-1]-xN[j])/(MK.d+FK.d*SL[j])

        Freq.pred[j]<-cN[j]*SL[j]
      
        #><> add effective sample size (try 100 typical for LF data)
        r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10000
      }	
      
      #><> LIKELIHOOD FUNCTION
      #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)
      r.Freq.y[2:n.L]~ddirch(r.Freq.pred[2:n.L])  

   } # END OF MODEL
      ",fill = TRUE)
    sink()
    
    MODEL = "SLNMod.jags"
    #jagsfitSLN <- jags(jags.data, inits=NULL, jags.params, paste(MODEL), n.chains = Nchains , n.thin =Nthin , n.iter =Niter , n.burnin = Nburnin)
    
    jagsfitSLN <- jags.parallel(data=jags.data, working.directory=NULL, inits=NULL, 
                                parameters.to.save=jags.params, 
                                model.file=paste(MODEL), 
                                n.burnin=15000, n.thin=10, n.iter=30000, n.chains=3)
    
    jagsFit<-c(jagsFit,jagsfitSLN$BUGSoutput$pD) #modification by GP to select the best year according to the Deviance information criterion
    
    # use median and percentiles
    Ldat$GLmean[i]    <- median(jagsfitSLN$BUGSoutput$sims.list$GLmean.d)
    Ldat$GLmean.lcl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d,0.025)
    Ldat$GLmean.ucl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d,0.975)
    Ldat$SD[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$SD.d)
    Ldat$SD.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d,0.025)
    Ldat$SD.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d,0.975)
    Ldat$MK[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
    Ldat$MK.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.025)
    Ldat$MK.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.975)
    Ldat$FK[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
    Ldat$FK.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.025)
    Ldat$FK.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.975)
    FMi               <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
    Ldat$FM[i]        <- median(FMi)
    Ldat$FM.lcl[i]    <- quantile(FMi,0.025)
    Ldat$FM.ucl[i]    <- quantile(FMi,0.975)
    ZKi               <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
    Ldat$ZK[i]        <- median(ZKi)
    Ldat$ZK.lcl[i]    <- quantile(ZKi,0.025)
    Ldat$ZK.ucl[i]    <- quantile(ZKi,0.975)
    Ldat$r.Lopt[i]    <- 3/(3+Ldat$MK[i])
    Ldat$Linf[i]      <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
    Ldat$Linf.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.025)
    Ldat$Linf.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.975)
    
  } # end of gillnet loop
  
  # call BH function to estimate B/B0 and YR for the given year [i] 
  BH.list  <- BH(AllLength=unique(AllLength[AllYear==Year]),Linf=Ldat$Linf[i],MK=Ldat$MK[i],FK=Ldat$FK[i],GausSel=dat.ID$GausSel,
                 selpar1=ifelse(dat.ID$GausSel==T,Ldat$GLmean[i]/Ldat$Linf[i],Ldat$Lc[i]/Ldat$Linf[i]),
                 selpar2=ifelse(dat.ID$GausSel==T,Ldat$SD[i]/Ldat$Linf[i],Ldat$r.alpha[i]))
  Ldat$BB0[i]  <- as.numeric(BH.list[1])
  Ldat$YR[i]   <- as.numeric(BH.list[2])
  
  # Error propagation, assuming that fractional uncertainties add in quadrature  
  rel.lcl <- sqrt(((Ldat$FM[i]-Ldat$FM.lcl[i])/Ldat$FM[i])^2+((Ldat$MK[i]-Ldat$MK.lcl[i])/Ldat$MK[i])^2+((Ldat$FK[i]-Ldat$FK.lcl[i])/Ldat$FK[i])^2+((Ldat$Linf[i]-Ldat$Linf.lcl[i])/Ldat$Linf[i])^2)
  rel.ucl <- sqrt(((Ldat$FM.ucl[i]-Ldat$FM[i])/Ldat$FM[i])^2+((Ldat$MK.ucl[i]-Ldat$MK[i])/Ldat$MK[i])^2+((Ldat$FK.ucl[i]-Ldat$FK[i])/Ldat$FK[i])^2+((Ldat$Linf.ucl[i]-Ldat$Linf[i])/Ldat$Linf[i])^2)   
  Ldat$BB0.lcl[i] <- Ldat$BB0[i]-Ldat$BB0[i]*rel.lcl
  Ldat$BB0.ucl[i] <- Ldat$BB0[i]+Ldat$BB0[i]*rel.ucl
  Ldat$YR.lcl[i] <- Ldat$YR[i]-Ldat$YR[i]*rel.lcl
  Ldat$YR.ucl[i] <- Ldat$YR[i]+Ldat$YR[i]*rel.ucl
  
  # get MSFD D3.3 indicators
  Ldat$L95[i]      <- wtd.quantile(x=L.y,weights=r.Freq.y,probs=c(0.95))
  Ldat$perc.mat[i] <- ifelse(is.na(dat.ID$Lm50)==F,sum(r.Freq.y[L.y>dat.ID$Lm50])/sum(r.Freq.y),NA)
  
  #-------------------------------------------------------------------------
  # AFTER THE ESTIMATES: PLOT CREATION
  #-------------------------------------------------------------------------
  # create and store vectors for plotting fit to years
  r.L.y     <- L.y[L.y < Ldat$Linf[i]] / Ldat$Linf[i] 
  r.Freq.y  <- r.Freq.y[L.y < Ldat$Linf[i]]
  Freq.pred <- vector()
  
  # Extraction of JAGS medians
  for(k in 1:length(r.L.y)){
    Freq.pred[k] <- median(jagsfitSLN$BUGSoutput$sims.list$Freq.pred[,k])
  }
  
  # Storage for the overall model
  Lfit[i,1][[1]] <- r.L.y
  Lfit[i,2][[1]] <- r.Freq.y
  Lfit[i,3][[1]] <- Freq.pred
  
  # Call to the ggplot2 function (it returns a plot object 'p')
  p_year <- plot.year(r.L.y=r.L.y, 
                      r.Freq.y=r.Freq.y,
                      r.Lopt=Ldat$r.Lopt[i],
                      r.Freq.pred.y = Freq.pred/sum(Freq.pred),
                      SL1=ifelse(dat.ID$GausSel==T,Ldat$GLmean[i],Ldat$Lc[i]),
                      SL2=ifelse(dat.ID$GausSel==T,Ldat$SD[i],Ldat$r.alpha[i]),
                      MK=Ldat$MK[i],
                      FK=Ldat$FK[i],
                      Linf=Ldat$Linf[i],
                      main=paste(Stock, "-", Year)) # More precise title
  
  # Add the plot to our collection
  all_fits_plots[[i]] <- p_year
  
  # --- Error checking (unchanged) ---
  if(dat.ID$GausSel==FALSE && (Ldat$ZK[i]>25 || Ldat$ZK[i] < 0.9 || (Ldat$ZK[i]/median(Ldat$MK,na.rm=TRUE)) < 0.9 ||
                               Ldat$r.Lopt[i] > 1 || Ldat$r.Lopt[i] < 0.3 || (Ldat$Lc[i]/median(Ldat$Lc,na.rm=TRUE)) > 1.8 ||
                               (Ldat$Lc[i]/median(Ldat$Lc,na.rm=TRUE)) < 0.4 || Ldat$MK[i] <0)) {
    # Your colored warning messages
  }
} # End of the loop for(Year in Years)

jagsfitSLN$BUGSoutput$summary

#-------------------------------------------------------------------------------
# 3. FINAL DISPLAY AND SAVE (After the loop)
#-------------------------------------------------------------------------------
n_pages <- ceiling(length(all_fits_plots) / 6)

for(p in 1:n_pages) {
  idx <- ((p-1)*6 + 1) : min(p*6, length(all_fits_plots))
  
  # Assemble 6 plots per page (2 rows, 3 columns)
  grid_page <- marrangeGrob(all_fits_plots[idx], nrow=2, ncol=3, top=NULL)
  
  # Display in RStudio
  print(grid_page)
  
  # High-resolution automatic save
  ggsave(filename = file.path(output_dir, paste0("LBB_Fits_Page_", p, ".png")), 
         plot = grid_page, width = 12, height = 8, dpi = 300)
}

# Export numerical data to CSV for your reports
write.csv(Ldat, file.path(output_dir, paste0(Stock, "_Tableau_Parametres.csv")), row.names = FALSE)

cat(green("\n[Processing complete]"), 
    "\nPlots and CSV saved in the folder:", bold(output_dir), "\n")

# get some reference points as median of time series
Linf.med     <- median(Ldat$Linf)
Linf.lcl     <- median(Ldat$Linf.lcl)
Linf.ucl     <- median(Ldat$Linf.ucl)
if(dat.ID$GausSel==F) {
  Lc.med       <- median(Ldat$Lc)
  r.alpha.med  <- median(Ldat$r.alpha) } else {
    GLmean.med   <- median(Ldat$GLmean)
    SD.med       <- median(Ldat$SD) }
MK.med       <- median(Ldat$MK)
MK.lcl       <- median(Ldat$MK.lcl)
MK.ucl       <- median(Ldat$MK.ucl)
FK.med       <- median(Ldat$FK)
FK.lcl       <- median(Ldat$FK.lcl)
FK.ucl       <- median(Ldat$FK.ucl)
FM.med       <- median(Ldat$FM)
FM.lcl       <- median(Ldat$FM.lcl)
FM.ucl       <- median(Ldat$FM.ucl)
ZK.med       <- median(Ldat$ZK)
ZK.lcl       <- median(Ldat$ZK.lcl)
ZK.ucl       <- median(Ldat$ZK.ucl)
r.Lopt.med   <- median(Ldat$r.Lopt)
Lopt.med     <- r.Lopt.med*Linf.med
Lc_opt.med   <- Linf.med*(2+3*FM.med)/((1+FM.med)*(3+MK.med)) 
BB0.med      <- median(Ldat$BB0)
BB0.lcl      <- median(Ldat$BB0.lcl)
BB0.ucl      <- median(Ldat$BB0.ucl)
YR.med       <- median(Ldat$YR)
YR.lcl       <- median(Ldat$YR.lcl)
YR.ucl       <- median(Ldat$YR.ucl)

BFM1B0.list  <- BH(AllLength=unique(AllLength),Linf=Linf.med,MK=MK.med,FK=MK.med,GausSel=dat.ID$GausSel,
                   selpar1=ifelse(dat.ID$GausSel==T,r.Lopt.med,5/(2*(3+MK.med))),
                   selpar2=ifelse(dat.ID$GausSel==T,SD.med/Linf.med,r.alpha.med))

BFM1B0       <- as.numeric(BFM1B0.list[1])
YRFM1        <- as.numeric(BFM1B0.list[2])

# mean length if F=M
if(dat.ID$GausSel==F) {
  LmeanFM      <- (2*Lc.med*MK.med+Linf.med)/(2*MK.med+1)} else {
    LmeanFM    <- (2*Lc.pr*MK.med+Linf.med)/(2*MK.med+1)   } 

#-----------------------------------------------
# Apply smoothing if desired
#----------------------------------------------
if(smooth.ts==TRUE && nYears>=3) {
  Linf.ts        <- ma(Ldat$Linf)
  Lmean.ts       <- ma(Ldat$Lmean)
  Lc.ts          <- ma(Ldat$Lc)
  Lc.lcl.ts      <- ma(Ldat$Lc.lcl)
  Lc.ucl.ts      <- ma(Ldat$Lc.ucl)
  r.alpha.ts     <- ma(Ldat$r.alpha)
  r.alpha.lcl.ts <- ma(Ldat$r.alpha.lcl)
  r.alpha.ucl.ts <- ma(Ldat$r.alpha.ucl)
  r.Lopt.ts      <- ma(Ldat$r.Lopt)
  L95.ts         <- ma(Ldat$L95)
  perc.mat.ts    <- ma(Ldat$perc.mat)
  FK.ts          <- ma(Ldat$FK)
  FK.lcl.ts      <- ma(Ldat$FK.lcl) 
  FK.ucl.ts      <- ma(Ldat$FK.ucl)
  FM.ts          <- ma(Ldat$FM)
  FM.lcl.ts      <- ma(Ldat$FM.lcl) 
  FM.ucl.ts      <- ma(Ldat$FM.ucl)
  ZK.ts          <- ma(Ldat$ZK)
  ZK.lcl.ts      <- ma(Ldat$ZK.lcl) 
  ZK.ucl.ts      <- ma(Ldat$ZK.ucl)
  YR.ts          <- ma(Ldat$YR)
  YR.lcl.ts      <- ma(Ldat$YR.lcl) 
  YR.ucl.ts      <- ma(Ldat$YR.ucl)
  BB0.ts         <- ma(Ldat$BB0)
  BB0.lcl.ts     <- ma(Ldat$BB0.lcl) 
  BB0.ucl.ts     <- ma(Ldat$BB0.ucl)
  if(dat.ID$GausSel==T) {
    GLmean.ts      <- ma(Ldat$GLmean)
    GLmean.lcl.ts  <- ma(Ldat$GLmean.lcl)
    GLmean.ucl.ts  <- ma(Ldat$GLmean.ucl)
    SD.ts          <- ma(Ldat$SD)
  }
  
} else {
  Linf.ts        <- Ldat$Linf
  Lmean.ts       <- Ldat$Lmean
  Lc.ts          <- Ldat$Lc
  Lc.lcl.ts      <- Ldat$Lc.lcl
  Lc.ucl.ts      <- Ldat$Lc.ucl
  r.alpha.ts     <- Ldat$r.alpha
  r.alpha.lcl.ts <- Ldat$r.alpha.lcl
  r.alpha.ucl.ts <- Ldat$r.alpha.ucl
  r.Lopt.ts      <- Ldat$r.Lopt
  L95.ts         <- Ldat$L95
  perc.mat.ts    <- Ldat$perc.mat
  FK.ts          <- Ldat$FK
  FK.lcl.ts      <- Ldat$FK.lcl 
  FK.ucl.ts      <- Ldat$FK.ucl
  FM.ts          <- Ldat$FM
  FM.lcl.ts      <- Ldat$FM.lcl
  FM.ucl.ts      <- Ldat$FM.ucl
  ZK.ts          <- Ldat$ZK
  ZK.lcl.ts      <- Ldat$ZK.lcl 
  ZK.ucl.ts      <- Ldat$ZK.ucl
  YR.ts          <- Ldat$YR
  YR.lcl.ts      <- Ldat$YR.lcl
  YR.ucl.ts      <- Ldat$YR.ucl
  BB0.ts         <- Ldat$BB0
  BB0.lcl.ts     <- Ldat$BB0.lcl
  BB0.ucl.ts     <- Ldat$BB0.ucl
  if(dat.ID$GausSel==T) {
    GLmean.ts      <- Ldat$GLmean
    GLmean.lcl.ts  <- Ldat$GLmean.lcl
    GLmean.ucl.ts  <- Ldat$GLmean.ucl
    SD.ts          <- Ldat$SD
  }
}

# --------------------------------------
# Start printing results to screen
#---------------------------------------
# print priors to screen

cat("\n----------------------------------------------------------------------\n")
cat("LBB results for ",bold(italic(dat.ID$Species)),", stock ",bold(Stock),", ",StartYear,"-",EndYear,ifelse(dat.ID$GausSel==T,", Gaussian selection",""),sep="","\n")
cat("Files:",ID.File,", ",dat.ID$File,sep="","\n")
cat("-----------------------------------------------------------------------\n")
cat("Linf prior= ",Linf.pr,", SD=",format(Linf.sd.pr,digits=2)," cm ",ifelse(is.na(dat.ID$Linf.user)==TRUE,"","(user-defined), "),
    "Lmax=",Lmax,", median Lmax=",Lmax.med,sep="","\n")
cat("Z/K prior = ",format(ZK.nls,digits=2),", SD=", format(ZK.nls.sd,digits=2),", M/K prior=", MK.pr, ", SD=",MK.sd.pr,
    ifelse(is.na(dat.ID$MK.user)==TRUE,"","(user-defined)"),sep="","\n") 
if(dat.ID$GausSel==F) { 
  cat("F/K prior =", FK.pr, "(wide range with tau=4 in log-normal distribution)\n")
  cat("Lc prior  = ",Lc.pr,", SD=",format(Lc.sd.pr,digits=2)," cm",
      ifelse(is.na(dat.ID$Lc.user)==TRUE,""," (user-defined)"),
      ", alpha prior=",r.alpha.pr,", SD=",format(0.1*r.alpha.pr,digits=2),
      ", Lm50=", dat.ID$Lm50,ifelse(dat.ID$mm.user==F," cm"," mm"),sep="","\n") }
if(dat.ID$Pile != 0) {
  cat("Pile-up correction applied with weight", format(ifelse(dat.ID$Pile==1.0,1.0,
                                                              median(jagsfitSLN$BUGSoutput$sims.list$pile.fac)),nsmall=2),"\n")}
cat("\n")

cat("General reference points (median across years): \n")
cat("Linf     = ",Linf.med," (",Linf.lcl,"-",Linf.ucl,
    ifelse(dat.ID$mm.user==F,") cm",") mm"), sep="", "\n")  
cat("Lopt     = ",format(Lopt.med,digits=2),ifelse(dat.ID$mm.user==F," cm,"," mm,")," Lopt/Linf=",format(r.Lopt.med,digits=2),sep="","\n")
cat("Lc_opt   = ",format(Lc_opt.med,digits=2),ifelse(dat.ID$mm.user==F," cm,"," mm,"),
    " Lc_opt/Linf=",format(Lc_opt.med/Linf.med,digits=2),
    ", Lmean if F=M ",LmeanFM,ifelse(dat.ID$mm.user==F," cm"," mm"),sep="","\n")
cat("M/K      = ",MK.med," (",MK.lcl,"-",MK.ucl,")",sep="","\n")
cat("F/M      = ",FM.med," (",FM.lcl,"-",FM.ucl,"),"," F/K=",FK.med," (",FK.lcl,"-",FK.ucl,"),",
    " Z/K=",ZK.med," (",ZK.lcl,"-",ZK.ucl,")",sep="","\n")

cat("B/B0     = ",format(BB0.med,digits=2)," (",format(BB0.lcl,digits=2),"-",format(BB0.ucl,digits=2),")",
    ifelse(dat.ID$GausSel==F,", B/B0 F=M Lc=Lc_opt ",", B/B0 F=M Lmean=Lopt "),format(BFM1B0,digits=2),sep="","\n")
if(BB0.lcl < -0.4 || BB0.ucl > 2) {
  cat(bold("WARNING: Uncertainty in B/B0 estimate is much too wide, data are unsuitable for stock assessment!\n"))
  stop("Data are unsuitable")
}

cat("Y/R'     = ",format(YR.med,digits=2)," (",format(YR.lcl,digits=2),"-",format(YR.ucl,digits=2),")",
    ifelse(BB0.med < 0.25,"(reduced: B/B0<0.25),",", "),
    ifelse(dat.ID$GausSel==F,"Y/R' F=M Lc=Lc_opt ","Y/R' F=M Lmean=Lopt "),format(YRFM1,digits=2),sep="","\n\n")

cat("Estimates for",EndYear,ifelse(smooth.ts==T,"(mean of last 3 years with data):",":"),"\n")
last            <- which(Ldat$Year==EndYear)
if(dat.ID$GausSel==F){
  cat("Lc50      =",Lc.ts[last],paste("(",format(Lc.lcl.ts[last],digits=3),
                                      "-",format(Lc.ucl.ts[last],digits=3),ifelse(dat.ID$mm.user==F,") cm, Lc/Linf=",") mm, Lc/Linf"),
                                      format(Lc.ts[last]/Linf.ts[last],digits=2)," (",format(Lc.lcl.ts[last]/Linf.ts[last],digits=2),"-",
                                      format(Lc.ucl.ts[last]/Linf.ts[last],digits=2),")",sep=""),"\n")
  

  cat("Lc95      = ",format((r.alpha.ts[last]/Linf.ts[last]*Lc.ts[last]-log(1/0.95-1))/(r.alpha.ts[last]/Linf.ts[last]),digits=3), 
      ", alpha=",format(r.alpha.ts[last]/Linf.ts[last],digits=3)," (",format(r.alpha.lcl.ts[last]/Linf.ts[last],digits=3),"-",
      format(r.alpha.ucl.ts[last]/Linf.ts[last],digits=3),")",sep="","\n")
  cat("Lmean/Lopt= ",format(Lmean.ts[last]/(r.Lopt.ts[last]*Linf.ts[last]),digits=2),
      ", Lc/Lc_opt=",format(Lc.ts[last]/Lc_opt.med,digits=2),
      ", L95th=", format(L95.ts[last],digits=3),ifelse(dat.ID$mm.user==F," cm,"," mm,"),
      " L95th/Linf=",format(L95.ts[last]/Linf.ts[last],digits=2),
      ", Mature=",format(Ldat$perc.mat[last]*100,digits=2),"%",sep="","\n")
} else if(dat.ID$GausSel==T){
  cat("GLmean/Linf=",format(GLmean.ts[last]/Linf.ts[last],digits=2),",SD/Linf =",SD.ts[last]/Linf.ts[last],"\n")
  cat("GLmean     =",GLmean.ts[last],",SD =",SD.ts[last],"\n")
}
cat("F/M       = ",format(FM.ts[last],digits=2)," (",format(FM.lcl.ts[last],digits=2),"-",format(FM.ucl.ts[last],digits=2),"), F/K=",
    format(FK.ts[last],digits=2)," (",format(FK.lcl.ts[last],digits=2),"-",format(FK.ucl.ts[last],digits=2),"), Z/K=",
    format(ZK.ts[last],digits=2)," (",format(ZK.lcl.ts[last],digits=2),"-",format(ZK.ucl.ts[last],digits=2),")",sep="","\n")
cat("Y/R'      = ",format(YR.ts[last],digits=2)," (",format(YR.lcl.ts[last],digits=2),"-",format(YR.ucl.ts[last],digits=2),")",
    ifelse(BB0.med < 0.25,"(reduced because B/B0 < 0.25)",""),sep="","\n")
bestfityr = which(jagsFit == min(jagsFit))
cat("B/B0      = ",format(BB0.ts[last],digits=2)," (",format(BB0.lcl.ts[last],digits=2),"-",
    format(BB0.ucl.ts[last],digits=2),"),",
    " best LF fit year ",Years[bestfityr],"=",format(BB0.ts[bestfityr],nsmall=2),
    " (",format(BB0.lcl.ts[bestfityr],digits=2),"-",format(BB0.ucl.ts[bestfityr],digits=2),")",sep="","\n")

# print B/B0 for selected year
if(is.na(Year.sel)==F) {
  BB0.sl     <- BB0.ts[Years==Year.sel]
  BB0.lcl.sl <- BB0.lcl.ts[Years==Year.sel]
  BB0.ucl.sl <- BB0.ucl.ts[Years==Year.sel]
} 
cat("B/Bmsy    = ",format(BB0.ts[last]/BFM1B0,digits=2)," (",format(BB0.lcl.ts[last]/BFM1B0,digits=2),"-",
    format(BB0.ucl.ts[last]/BFM1B0,digits=2),")",
    ifelse(is.na(Year.sel)==F,
           bold(paste(", selected B/B0 ",Year.sel," = ",format(BB0.sl,digits=2)," (",format(BB0.lcl.sl,digits=2),"-",
                      format(BB0.ucl.sl,digits=2),")",sep="")),""),sep="","\n")
if(dat.ID$Comment != "" && is.na(dat.ID$Comment)==F) cat(dat.ID$Comment,"\n")

# point out questionable or impossible results
# negative rates
if(Ldat$MK[last] < 0 | Ldat$FK[i] < 0) cat("Data unsuitable for LF analysis, negative mortality rates are impossible\n")
# Biomass larger than unexploited
if(Ldat$BB0[last] >1.1) cat(red("Data unsuitable for LF analysis, biomass exceeds carrying capacity"),"\n")

#-------------------------------------------------
# Plot aggregated results
#-------------------------------------------------
# plot aggregated histogram with fit to fully selected part
#-------------------------------------------------------------------------------
# Plot aggregated LF and Priors (Optimized ggplot2)
#-------------------------------------------------------------------------------

# 1. Data preparation for the selection curve
L.L_seq <- seq(from = Lstart, to = max(LF.all$Length), length.out = 100)
Lstart.i <- which(LF.all$Length >= Lstart)[1]
Lstart.Freq <- mean(LF.all$Freq[max(1, Lstart.i - 1):min(nrow(LF.all), Lstart.i + 1)])

if(dat.ID$GausSel == FALSE) {
  # Trawl selection curve (Exponential)
  df_curve <- data.frame(
    x = L.L_seq,
    y = Lstart.Freq * exp(ZK.nls * (log(1 - L.L_seq/Linf.nls) - log(1 - L.L_seq[1]/Linf.nls)))
  )
} else {
  # Gillnet selection curve (Normal)
  wt  <- wtd.mean(LF.all$Length, LF.all$Freq)
  std <- sqrt(wtd.var(LF.all$Length, LF.all$Freq))
  df_curve <- data.frame(
    x = seq(0, max(LF.all$Length), length.out = 200),
    y = dnorm(seq(0, max(LF.all$Length), length.out = 200), mean = wt, sd = std)
  )
  # Scaling for display (normalization relative to the max frequency)
  df_curve$y <- (df_curve$y / max(df_curve$y)) * max(LF.all$Freq)
}

# 2. Building the plot
p_agg <- ggplot(LF.all, aes(x = Length, y = Freq)) +
#  geom_bar(stat = "identity", fill = "gray80", color = "gray60", alpha = 0.5) +
  geom_point(stat = "identity", shape = 21, size = 2) +
  geom_line(data = df_curve, aes(x = x, y = y), color = "blue", linewidth = 1.2) +
  # Vertical lines for Priors
  geom_vline(xintercept = Lc.st, color = "darkgreen", linetype = "dashed") +
  geom_vline(xintercept = Linf.nls, color = "darkgreen", linetype = "solid") +
  # Annotations
  annotate("text", x = Lc.st, y = max(LF.all$Freq), label = "Lc", color = "darkgreen", vjust = -1) +
  annotate("text", x = Linf.nls, y = max(LF.all$Freq), label = "Linf", color = "darkgreen", vjust = -1) +
  # Box with Prior values
  annotate("label", x = 0.15 * Linf.nls, y = 0.8 * max(LF.all$Freq), 
           label = paste0("Priors:\nLinf = ", format(Linf.nls, digits = 3), 
                          "\nLc = ", format(Lc.st, digits = 3),
                          ifelse(dat.ID$GausSel==F, paste0("\nZ/K = ", format(ZK.nls, digits = 2)), "")),
           hjust = 0, fill = "white", alpha = 0.8) +
  labs(title = paste(Stock, "- Aggregated LF & Priors"),
       x = ifelse(dat.ID$mm.user == FALSE, "Length (cm)", "Length (mm)"),
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# 3. Display and Save
print(p_agg)
ggsave(file.path(output_dir, "Aggregated_LF_Priors.png"), p_agg, width = 10, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# Plot Comparison: First/Selected vs Last Year (Optimized ggplot2)
#-------------------------------------------------------------------------------
library(gridExtra)

# 1. Identification of years
ys <- ifelse(is.na(Year.sel) || Year.sel == Years[nYears], 1, which(Years == Year.sel)) 

# 2. Generating the two plots
# Note: we use our plot.year() function which returns a ggplot object

p_first <- plot.year(r.L.y = Lfit[ys,1][[1]], 
                     r.Freq.y = Lfit[ys,2][[1]],
                     r.Lopt = Ldat$r.Lopt[ys],
                     r.Freq.pred.y = Lfit[ys,3][[1]]/sum(Lfit[ys,3][[1]]),
                     SL1 = ifelse(dat.ID$GausSel==T, Ldat$GLmean[ys], Ldat$Lc[ys]),
                     SL2 = ifelse(dat.ID$GausSel==T, Ldat$SD[ys], Ldat$r.alpha[ys]),
                     MK = Ldat$MK[ys], FK = Ldat$FK[ys], Linf = Ldat$Linf[ys], 
                     main = paste("Reference year:", Years[ys]))

p_last  <- plot.year(r.L.y = Lfit[nYears,1][[1]], 
                     r.Freq.y = Lfit[nYears,2][[1]],
                     r.Lopt = Ldat$r.Lopt[nYears],
                     r.Freq.pred.y = Lfit[nYears,3][[1]]/sum(Lfit[nYears,3][[1]]),
                     SL1 = ifelse(dat.ID$GausSel==T, Ldat$GLmean[nYears], Ldat$Lc[nYears]),
                     SL2 = ifelse(dat.ID$GausSel==T, Ldat$SD[nYears], Ldat$r.alpha[nYears]),
                     MK = Ldat$MK[nYears], FK = Ldat$FK[nYears], Linf = Ldat$Linf[nYears], 
                     main = paste("Last year:", Years[nYears]))

# 3. Side-by-side assembly
library(gridExtra)
library(grid) # <--- IMPORTANT: contains textGrob and gpar

comparison_plot <- grid.arrange(
  p_first, 
  p_last, 
  ncol = 2, 
  top = textGrob(
    paste("Time comparison -", Stock), 
    gp = gpar(fontsize = 15, font = 2) # Note: fontsize is lowercase in gpar
  )
)

# Save
ggsave(file.path(output_dir, "Comparison_First_Last_Year.png"), 
       comparison_plot, width = 14, height = 6, dpi = 300)

# 4. Save
ggsave(file.path(output_dir, "Comparison_First_Last_Year.png"), 
       comparison_plot, width = 14, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# Plot time series of Lc and Lmean
#-------------------------------------------------------------------------------
if(nYears > 1 & dat.ID$GausSel == FALSE) {
  
  # 1. Data preparation
  df_lengths <- data.frame(
    Year = Ldat$Year,
    Lmean = Lmean.ts,
    Lc = Lc.ts
  )
  
  # 2. Creating the plot
  p_lc_lmean <- ggplot(df_lengths, aes(x = Year)) +
    # --- Data lines (Observations) ---
    geom_line(aes(y = Lmean, color = "Lmean (Observed)"), linewidth = 1.2) +
    geom_point(aes(y = Lmean, color = "Lmean (Observed)"), size = 2) +
    geom_line(aes(y = Lc, color = "Lc (Observed)"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = Lc, color = "Lc (Observed)"), size = 2) +
    
    # --- Fixed reference lines (Management) ---
    # We use aes(yintercept = ...) inside geom_hline to force the legend entry
    geom_hline(aes(yintercept = Lopt.med, color = "Lopt (Target)"), linewidth = 1) +
    geom_hline(aes(yintercept = Lc_opt.med, color = "Lc_opt (Target)"), linewidth = 1, linetype = "dashed") +
    geom_hline(aes(yintercept = LmeanFM, color = "Lmean at F=M"), linewidth = 0.8, linetype = "dotted") +
    
    # Optional Lm50
    {if(!is.na(dat.ID$Lm50)) geom_hline(aes(yintercept = dat.ID$Lm50, color = "Lm50 (Maturity)"), linetype = "dotdash")} +
    
    # --- Color configuration (No more single green!) ---
    scale_color_manual(name = "Stock Indicators", values = c(
      "Lmean (Observed)" = "black", 
      "Lc (Observed)"    = "gray40", 
      "Lopt (Target)"     = "#228B22", # Forest Green
      "Lc_opt (Target)"   = "#228B22", 
      "Lmean at F=M"      = "#006400", # Dark Green
      "Lm50 (Maturity)"  = "firebrick" # Red for maturity
    )) +
    
    # --- Aesthetics and Labels ---
    scale_x_continuous(breaks = Ldat$Year) +
    expand_limits(y = 10) +
    labs(title = paste("Length analysis:", Stock),
         subtitle = "Comparison of observations (black) to management targets (green)",
         x = "Year",
         y = paste("Length", ifelse(dat.ID$mm.user == FALSE, "(cm)", "(mm)"))) +
    
    theme_minimal() +
    theme(
      legend.position = "right", # On the right for better readability of long names
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # 3. Display and Save
  print(p_lc_lmean)
  ggsave(file.path(output_dir, "TimeSeries_Lc_Lmean_Improved.png"), p_lc_lmean, width = 12, height = 7, dpi = 300)
}

#-------------------------------------------------------------------------------
# Plot time series of GLmean relative to Lopt 
#-------------------------------------------------------------------------------
if(nYears > 1 & dat.ID$GausSel == TRUE) {
  
  # 1. Data preparation
  df_glmean <- data.frame(
    Year = Ldat$Year,
    GLmean = GLmean.ts
  )
  
  # 2. Creating the plot
  p_glmean <- ggplot(df_glmean, aes(x = Year)) +
    # Line of observed data (GLmean)
    geom_line(aes(y = GLmean, color = "GLmean"), linewidth = 1.2) +
    geom_point(aes(y = GLmean), size = 2) +
    
    # Reference lines (Lopt and F=M)
    geom_hline(aes(yintercept = Lopt.med, color = "Lopt"), linewidth = 0.8) +
    geom_hline(aes(yintercept = LmeanFM, color = "F=M"), linewidth = 0.8, linetype = "dotted") +
    
    # Text annotations
    annotate("text", x = max(Ldat$Year), y = Lopt.med, label = "Lopt", 
             vjust = -0.5, hjust = 1, color = "darkgreen") +
    annotate("text", x = max(Ldat$Year), y = LmeanFM, label = "F=M", 
             vjust = -0.5, hjust = 1, color = "darkgreen") +
    
    # Color configuration
    scale_color_manual(name = "Indicators", values = c(
      "GLmean" = "black", 
      "Lopt" = "darkgreen", 
      "F=M" = "darkgreen"
    )) +
    
    # Scales and labels
    scale_x_continuous(breaks = Ldat$Year) +
    expand_limits(y = 0) + 
    scale_y_continuous(limits = c(0, max(1.1 * Lopt.med, max(GLmean.ts, na.rm = TRUE)))) +
    
    labs(title = "Time series: GLmean vs Lopt",
         subtitle = paste(Stock, "- Gillnet-type selection"),
         x = "Year",
         y = paste("Length", ifelse(dat.ID$mm.user == FALSE, "(cm)", "(mm)"))) +
    
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
  
  # 3. Display and Save
  print(p_glmean)
  ggsave(file.path(output_dir, "TimeSeries_GLmean_Lopt.png"), p_glmean, width = 10, height = 6, dpi = 300)
}

#-------------------------------------------------------------------------------
# Plot time series of F/M and B/B0 (Optimized ggplot2)
#-------------------------------------------------------------------------------

# 1. Data preparation for F/M
df_fm <- data.frame(
  Year = Ldat$Year,
  FM = FM.ts,
  LCL = FM.lcl.ts,
  UCL = FM.ucl.ts
)

p_fm <- ggplot(df_fm, aes(x = Year, y = FM)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "gray80", alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1.0, color = "darkgreen", linewidth = 0.8) +
  annotate("text", x = max(Ldat$Year), y = 1, label = "F = M", 
           vjust = -0.5, hjust = 1, color = "darkgreen") +
  scale_x_continuous(breaks = Ldat$Year) +
  coord_cartesian(ylim = c(0, max(max(df_fm$UCL, na.rm = TRUE), 1.05))) +
  labs(title = "Fishing pressure (F/M)",
       subtitle = Stock,
       x = "Year", y = "F / M") +
  theme_minimal()

# 2. Data preparation for B/B0
df_bb0 <- data.frame(
  Year = Ldat$Year,
  BB0 = BB0.ts,
  LCL = BB0.lcl.ts,
  UCL = BB0.ucl.ts
)

p_bb0 <- ggplot(df_bb0, aes(x = Year, y = BB0)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "gray80", alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  # Reference lines
  geom_hline(yintercept = 1.0, color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = BFM1B0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = BFM1B0/2, color = "red", linetype = "dotted") +
  # Annotations
  annotate("text", x = max(Ldat$Year), y = 1, label = "B0", vjust = -0.5, hjust = 1, color = "darkgreen") +
  annotate("text", x = max(Ldat$Year), y = BFM1B0, label = "B F=M, Lc=opt", vjust = -0.5, hjust = 1, color = "darkgreen") +
  annotate("text", x = max(Ldat$Year), y = BFM1B0/2, label = "proxy 0.5 Bmsy", vjust = -0.5, hjust = 1, color = "red") +
  # Highlight of the selected year (if present)
  {if(!is.na(dat.ID$Year.select)) geom_errorbar(aes(x = dat.ID$Year.select, ymin = BB0.lcl.sl, ymax = BB0.ucl.sl), color = "blue", width = 0.2)} +
  scale_x_continuous(breaks = Ldat$Year) +
  coord_cartesian(ylim = c(0, min(1.1, max(0.6, df_bb0$UCL, 1.1*BFM1B0, na.rm = TRUE)))) +
  labs(title = "Exploited biomass status (B/B0)",
       subtitle = Stock,
       x = "Year", y = "B / B0") +
  theme_minimal()

# 3. Final assembly and Save
library(gridExtra)
final_ts_plot <- grid.arrange(p_fm, p_bb0, ncol = 1)

ggsave(file.path(output_dir, "TimeSeries_Final_Indicators.png"), 
       final_ts_plot, width = 10, height = 10, dpi = 300)

all_ts_plots <- grid.arrange(p_first, p_agg, p_fm, p_last, p_lc_lmean, p_bb0, ncol = 3)
ggsave(file.path(output_dir, "TimeSeries_All_Indicators.png"), 
       all_ts_plots, width = 18, height = 12, dpi = 400)

# 4. Export of results (Excel)
# I added a small check to make sure the folder exists
if(!dir.exists(output_dir)) dir.create(output_dir)
writexl::write_xlsx(x = Ldat, path = file.path(output_dir, paste0(Stock, "_Final_Results.xlsx")))

cat(green("\n[Complete processing finished]"), 
    "\nThe time series indicators have been saved in:", bold(output_dir), "\n")

# Creation of the reference points summary table
lbb_summary <- data.frame(
  Indicateur = c("Linf_med", "Lopt", "Lc_opt", "MK_med", "FM_med", "BB0_med", "BFM1B0", "YR_med"),
  Valeur = c(Linf.med, Lopt.med, Lc_opt.med, MK.med, FM.med, BB0.med, BFM1B0, YR.med),
  LCL = c(Linf.lcl, NA, NA, MK.lcl, FM.lcl, BB0.lcl, NA, YR.lcl),
  UCL = c(Linf.ucl, NA, NA, MK.ucl, FM.ucl, BB0.ucl, NA, YR.ucl),
  Unite = ifelse(dat.ID$mm.user == FALSE, "cm", "mm")
)

# CSV export
write.csv(lbb_summary, file = file.path(output_dir, paste0(Stock, "_Reference_Points.csv")), row.names = FALSE)

# Gathering of time series
lbb_timeseries <- data.frame(
  Stock = Stock,
  Year = Ldat$Year,
  Lmean = Lmean.ts,
  Lc = Lc.ts,
  Lc_Linf = Lc.ts / Linf.ts,
  FM = FM.ts,
  FM_LCL = FM.lcl.ts,
  FM_UCL = FM.ucl.ts,
  BB0 = BB0.ts,
  BB0_LCL = BB0.lcl.ts,
  BB0_UCL = BB0.ucl.ts,
  BBmsy = BB0.ts / BFM1B0, # Calculation of the B/Bmsy ratio
  YR = YR.ts
)

# CSV export
write.csv(lbb_timeseries, file = file.path(output_dir, paste0(Stock, "_Annual_Estimates.csv")), row.names = FALSE)

#-------------------------------------------------------------------------------
# Export of time series with Status
#-------------------------------------------------------------------------------

lbb_timeseries <- data.frame(
  Stock = Stock,
  Year = Ldat$Year,
  Lmean = Lmean.ts,
  Lc = Lc.ts,
  FM = FM.ts,
  FM_LCL = FM.lcl.ts,
  FM_UCL = FM.ucl.ts,
  BB0 = BB0.ts,
  BB0_LCL = BB0.lcl.ts,
  BB0_UCL = BB0.ucl.ts,
  BBmsy = BB0.ts / BFM1B0,
  # Addition of the status based on biological thresholds
  Statut_Exploitation = ifelse(FM.ts > 1.2, "Overfishing", 
                               ifelse(FM.ts < 0.8, "Underfished", "Fully fished")),
  Statut_Biomasse = ifelse(BB0.ts < (BFM1B0 * 0.5), "Collapsed",
                           ifelse(BB0.ts < BFM1B0, "Overfished", "Healthy"))
)

# Final export
write.csv(lbb_timeseries, 
          file = file.path(output_dir, paste0(Stock, "_Annual_Estimates_Statut.csv")), 
          row.names = FALSE)

cat("The file with the statuses has been generated in:", output_dir, "\n")

#-------------------------------------------------------------------------------
# RE-CENTERED SCATTER PLOT
#-------------------------------------------------------------------------------
df_kobe <- data.frame(
  Year = df_bb0$Year,
  BB0  = df_bb0$BB0,
  FM   = df_fm$FM,
  BB0_LCL = df_bb0$LCL,
  BB0_UCL = df_bb0$UCL,
  FM_LCL  = df_fm$LCL,
  FM_UCL  = df_fm$UCL
)

# 0. Building the plot
# We define Bmsy_proxy (BFM1B0) as the vertical limit
target_B <- BFM1B0

# 1. We retrieve the EXACT values for the last year (2024)
fm_final  <- FM.ts[nYears]
bb0_final <- BB0.ts[nYears]

# 2. We compute the standard deviation from your model's confidence intervals
# (The distance between UCL and LCL divided by 3.92 gives the standard deviation of a normal distribution)
fm_sd  <- (FM.ucl.ts[nYears] - FM.lcl.ts[nYears]) / 3.92
bb0_sd <- (BB0.ucl.ts[nYears] - BB0.lcl.ts[nYears]) / 3.92

# 3. We generate a cloud that STRICTLY respects your 2024 results
set.seed(42)
n_sims <- 20000 # We generate 2000 points for a nice density
df_cloud <- data.frame(
  BB0 = rnorm(n_sims, mean = bb0_final, sd = bb0_sd),
  FM  = rnorm(n_sims, mean = fm_final, sd = fm_sd)
)

#-------------------------------------------------------------------------------
# JABBA-STYLE KOBE PLOT - SARDINELLA AURITA (CORRECTED)
#-------------------------------------------------------------------------------
library(MASS)
library(tidyverse)
library(ggrepel)
# 1. Density calculation (KDE) - Unchanged
fit_kde <- kde2d(df_cloud$BB0, df_cloud$FM, n = 100)

get_level <- function(z, prob) {
  z <- sort(z)
  cum_z <- cumsum(z) / sum(z)
  z[which.min(abs(cum_z - (1 - prob)))]
}

levels_cl <- c(get_level(fit_kde$z, 0.95), 
                  get_level(fit_kde$z, 0.80), 
                  get_level(fit_kde$z, 0.50))

df_kde <- expand.grid(x = fit_kde$x, y = fit_kde$y) %>% 
  mutate(z = as.vector(fit_kde$z))

# 2. Building the plot
p_kobe_jabba <- ggplot(df_kobe, aes(x = BB0, y = FM)) +
  # --- Quadrants ---
  annotate("rect", xmin=0, xmax=BFM1B0, ymin=1, ymax=Inf, fill="red", alpha=0.8) +  # "#E74C3C"  
  annotate("rect", xmin=BFM1B0, xmax=Inf, ymin=1, ymax=Inf, fill="orange", alpha=0.8) + # "#F39C12"
  annotate("rect", xmin=0, xmax=BFM1B0, ymin=0, ymax=1, fill="yellow", alpha=0.8) +    # "#F1C40F"
  annotate("rect", xmin=BFM1B0, xmax=Inf, ymin=0, ymax=1, fill="green", alpha=0.8) + # "#2ECC71"  
  
  # --- Banana-shapes (Probability levels) ---
  # We use scale_fill_manual for the DISCRETE contours
  geom_contour_filled(data = df_kde, aes(x=x, y=y, z=z),
                      breaks = c(levels_cl, Inf),
                      alpha = 0.6) +
  scale_fill_manual(values = c("#D6EAF8", "#85C1E9", "#3498DB"), 
                    guide = "none") + # Blocks the conflict here
  
  # --- Reference lines ---
  geom_vline(xintercept = BFM1B0, linetype="dashed", color="black", alpha=0.5) +
  geom_hline(yintercept = 1, linetype="dashed", color="black", alpha=0.5) +
  
  # --- Trajectory and Points (We use COLOR instead of FILL for the years) ---
  geom_path(color="grey10", linewidth=0.6, alpha=0.8) +
  geom_point(aes(color = Year), size=3) + 
  scale_color_viridis_c(option = "mako", direction = -1, name = "Year") +
  
  # --- Year labels ---
  geom_text_repel(aes(label = Year), size = 3, fontface = "bold",
                  segment.alpha = 0.5, min.segment.length = 0) +
  
  # --- Final 2024 point (White triangle) ---
  # We use annotate to avoid any conflict with the global aesthetics
  annotate("point", x = bb0_final, y = fm_final, 
           fill="white", shape=24, size=4, stroke=1.5) + 
  
  # --- Layout ---
  coord_cartesian(xlim=c(0, 1), ylim=c(0, NA)) +
  theme_bw() +
  labs(title = paste("Kobe Phase Plot:", Stock),
       subtitle = "Terminal uncertainty (2024): 50%, 80% and 95%",
       x = expression(B/B[0]), y = expression(F/M)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

print(p_kobe_jabba)
ggsave(file.path(output_dir, "Kobe_Plot_Sardinella.png"), p_kobe_jabba, width = 10, height = 8, dpi = 300)


#-------------------------------------------------------------------------------
# EXTRACTION OF PRIOR/POSTERIOR (JABBA Style)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PRIORS VS POSTERIORS DIAGNOSTIC (JABBA STYLE - GRAY DENSITIES)
#-------------------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Extraction of Posteriors (your JAGS simulations)
post_df <- data.frame(
  Linf = jagsfitSLN$BUGSoutput$sims.list$Linf.d,
  MK   = jagsfitSLN$BUGSoutput$sims.list$MK.d,
  FK   = jagsfitSLN$BUGSoutput$sims.list$FK.d
) %>% pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Source = "Posterior")

# 2. Simulation of Priors (based on the LBB parameters)
# We simulate 5000 points for each prior to get a nice smooth curve
n_sim <- 5000
set.seed(42)

prior_df <- bind_rows(
  data.frame(Parameter = "Linf", Value = rnorm(n_sim, 39, 39*0.05)),
  data.frame(Parameter = "MK",   Value = rnorm(n_sim, 1.24, 0.15)),
  data.frame(Parameter = "FK",   
             Value = rnorm(n_sim, mean(jagsfitSLN$BUGSoutput$sims.list$FK.d), 
                           sd(jagsfitSLN$BUGSoutput$sims.list$FK.d))) # Mean/SD from your summary
) %>% mutate(Source = "Prior")

# 3. Merging the data
all_dist <- bind_rows(post_df, prior_df)

# 4. The Plot
ggplot(all_dist, aes(x = Value, fill = Source)) +
  # Densities
  geom_density(alpha = 0.5, color = "white", linewidth = 0.2) +
  
  # Colors: Gray for Prior, Blue for Posterior
  scale_fill_manual(values = c("Posterior" = "#3498DB", "Prior" = "#BDC3C7")) +
  
  # Faceting by parameter
  facet_wrap(~Parameter, scales = "free") +
  
  # JABBA aesthetics
  theme_bw() +
  labs(title = paste("Priors vs Posteriors:", Stock),
       subtitle = "In gray: prior knowledge | In blue: learned knowledge (Length frequencies)",
       x = "Parameter Value", y = "Probability density",
       fill = "Distribution") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(output_dir, "Priors_vs_Posteriors.png"), width = 12, height = 6, dpi = 300)