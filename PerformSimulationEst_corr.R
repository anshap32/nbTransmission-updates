
################################################################################
# This program runs the full set of simulations by calling SimRunEst.R
################################################################################

rm(list = ls())
options(scipen=999)

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(caret)
library(gtools)
library(pROC)
library(devtools)
library(roxygen2)
library(stringr)
library(nbTransmission)
library(bindata)

#### Batch Mode ####

setwd("/projectnb/as-thesis/project2/sim_study")


#Getting sample size from the arguements
args <- commandArgs(trailingOnly = TRUE)
sampleSize <- as.numeric(args[1])
#Setting the observation date
if(as.numeric(args[2]) == 1){
  observationDate = "infectionDate"
  dateID = "ID"
}else if(as.numeric(args[2]) == 2){
  observationDate = "sampleDate"
  dateID = "SD"
}

lowerT <- as.numeric(args[3])
#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#The number of simulations per split
nSim <- 10



#### Interactive Mode ####
# 
# iTask <- 1
# sampleSize <- 50
# observationDate <- "infectionDate"
# dateID <- "ID"
# nSim <- 1
# iteration <- 1
# lowerT <- 2


#### Getting access to functions in other programs ####

source("SimRunEst.R")

source("nbSimulation/SimOutbreak.R")
source("nbSimulation/SimulateOutbreakS.R")
source("nbSimulation/SimEvaluate.R")
source("nbSimulation/SimCovariatesNBEffect_corr.R")

source("nbTransmission/nbProbabilities_sBS.R")
source("nbTransmission/runCV_sBS.R")
source("nbTransmission/performNB_BS_full.R")

source("nbTransmission/nbProbabilities_LR.R")
source("nbTransmission/runCV_LR.R")
source("nbTransmission/performNB_LR.R")





#### Setting the oubreak parameters ####

print(sampleSize)
print(observationDate)

#Effective population size times generation time (1 day or 1/365 years)
neg <- 0.25
pi <- 1
#Sets reproductive number to off.r if off.p = 0.5
off.r <- 1.2
off.p <- 0.5
#From Yincheng's work
w.shape <- 1.05
w.scale <- 1 / (0.0014 * 365)
w.shift <- 0.25
ws.shape <- w.shape
ws.scale <- w.scale
ws.shift <- w.shift
time <- 20
#Do you want multiple outbreaks
multOutbreaks <- TRUE 

#Genotype parameters
rootseq <- NULL
length <- 300
#This gives the equivalent of 0.5 mutations/genome/year
rate <- 0.5 / length
#Thresholds for SNP distance
thresholds <- c(lowerT, 12)

#Parameters used to manually run functions
covariates <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "timeCat")
pTraining <- 1
goldStandard <- "transmission"
truth <- "transmission"
pVar <- "pScaled"



####################### Running the Simulations ######################

est <- NULL
BS <- NULL
LR <- NULL

#Running the simuliaton for the outbreaks
for (iteration in 1:nSim){
  
  #Setting the seed for each run so I can re-run specific iterations with errors
  set.seed(iTask * 1000 + iteration)
  
  tryCatch({
    
    #Running the simulation
    runID <- paste(iTask, iteration, sampleSize, sep = "_")
    
    res <- simRunEst()
    
    est <- bind_rows(est, res$allEst %>% mutate(runID = runID))
    BS <- bind_rows(BS, res$resBS %>% mutate(runID = runID))
    LR <- bind_rows(LR, res$resLR%>% mutate(runID = runID))
    
    #Printing message that the run is finished
    print(paste0("Completed run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
    
  }, error = function(e){
    #Printing message that the run is finished
    print(paste0("Error in run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
    cat("ERROR :", conditionMessage(e), "\n")})
}

#Saving dataframes with summary of results for all of the runs
saveRDS(est, file=paste0("Corr_Simulation_Results/allEsts/estimates", "_", iTask, "_lowerT", lowerT, ".rds"))
saveRDS(BS, file=paste0("Corr_Simulation_Results/BSests/estimates", "_", iTask,  "_lowerT", lowerT,".rds"))
saveRDS(LR, file=paste0("Corr_Simulation_Results/LRests/estimates", "_", iTask,  "_lowerT", lowerT,".rds"))

