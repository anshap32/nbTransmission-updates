## Sarah's nbProbabilities function updated to include my weights
## calls runCV_LR which calls performNB_LR


nbProbabilities_sBS <- function(orderedPair, indIDVar, pairIDVar,
                            goldStdVar, covariates, label = "",
                            l = 1, n = 10, m = 1, nReps = 4, 
                            nBS = 100, pSampled = 1,
                            progressBar = TRUE){
  
  orderedPair <- as.data.frame(orderedPair)
  
  
  #### checking dataset is OK for function  #### 
  
  #Creating variables with the individual indID variable
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Checking that the named variables are in the data frame
  if(!indIDVar1 %in% names(orderedPair)){
    stop(paste0(indIDVar1, " is not in the data frame."))
  }
  if(!indIDVar2 %in% names(orderedPair)){
    stop(paste0(indIDVar2, " is not in the data frame."))
  }
  if(!pairIDVar %in% names(orderedPair)){
    stop(paste0(pairIDVar, " is not in the data frame."))
  }
  if(!goldStdVar %in% names(orderedPair)){
    stop(paste0(goldStdVar, " is not in the data frame."))
  }
  
  #Checking that the covariates are in the data frame
  covarTest <- covariates %in% names(orderedPair)
  if(FALSE %in% covarTest){
    stop("At least one of the covariates is not in the data frame.")
  }
  
  #Checking that all of the covariates are factors
  covarDf <- orderedPair[, covariates]
  notFactor <- names(covarDf)[!sapply(covarDf, is.factor)]
  notFactorC <- paste0(notFactor, collapse = ", ")
  if(length(notFactor) == 1){
    stop(paste0(notFactorC, " is not a factor"))
  }else if(length(notFactor > 1)){
    stop(paste0(notFactorC, " are not factors"))
  }
  
  
  #### Setting up data frames ####
  
  #Subsetting to only relevant columns to be more efficient
  orderedPair <- orderedPair[, c(indIDVar1, indIDVar2, pairIDVar,
                                 goldStdVar, covariates)]
  
  #Creating an pairID that where order doesn't matter (smaller ID always first)
  orderedPair$pairID_uo <- ifelse(orderedPair[, indIDVar1] < orderedPair[, indIDVar2],
                                  paste(orderedPair[, indIDVar1],
                                        orderedPair[, indIDVar2], sep = "_"),
                                  paste(orderedPair[, indIDVar2],
                                        orderedPair[, indIDVar1], sep = "_"))
  
  #Finding all pairs that can be included in the training dataset
  #And subsetting into the potential links
  posTrain <- orderedPair[!is.na(orderedPair[, goldStdVar]), ]
  
  
  
  #### Iterative Estimation Procedure ####
  
  #Initializing data frames to hold results (rAll) and coefficients (cAll)
  rAll <- data.frame("p" = numeric(), pairIDVar = character(), stringsAsFactors = FALSE)
  # names(rAll) <- c("p", pairIDVar)
  cAll <-list() #data.frame("level" = character(), "odds" = numeric(), stringsAsFactors = FALSE)
  pAll <- NULL
  
  
  if(progressBar == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = nReps, style = 3)
  }
  
  for (k in 1:nReps){
    
    #Randomly choosing the "true" infector from all possible
    #Calculating probabilities using mxn cross prediction
    cvResults <- runCV_sBS(posTrain, orderedPair, indIDVar, pairIDVar, # see function below
                          goldStdVar, covariates, l, n, m, nBS, pSampled) #, 
    rAll <- dplyr::bind_rows(rAll, cvResults$rFolds)
    # cAll <- dplyr::bind_rows(cAll, cvResults$cFolds)
    cAll[[k]] <- cvResults$cFolds
    # pAll <- dplyr::bind_rows(pAll, cvResults$pFolds)
    if(progressBar == TRUE){
      utils::setTxtProgressBar(pb, k)
    }
  }
  
  
  #### Summarizing Probabilities Over iterations ####
  
  #Averaging the probabilities over all the replicates
  sumData1 <- dplyr::group_by(rAll, !!rlang::sym(pairIDVar))
  sumData2 <- dplyr::summarize(sumData1,
                               pAvg = mean(!!rlang::sym("p"), na.rm = TRUE),
                               pSD = stats::sd(!!rlang::sym("p"), na.rm = TRUE),
                               nEstimates = sum(!is.na(!!rlang::sym("p"))),
                               label = dplyr::first(label),
                               .groups = 'drop')
  sumData2 <- dplyr::ungroup(sumData2)
  
  probs <- as.data.frame(dplyr::full_join(sumData2, orderedPair, by = pairIDVar),
                         stringsAsFactors = FALSE)
  
  #Calculating the total of all probabilities per infectee
  totalP <- stats::aggregate(probs$pAvg, by = list(probs[, indIDVar2]), sum, na.rm = TRUE)
  names(totalP) <- c(indIDVar2, "pTotal")
  
  #Calculating the scaled probabilities
  probs2 <- merge(probs, totalP, by = indIDVar2)
  probs2$pScaled <- ifelse(probs2$pTotal != 0, probs2$pAvg / probs2$pTotal, 0)
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  probs2 <- probs2[order(probs2[, indIDVar2], -probs2$pScaled), ]
  probs2$pRank <- stats::ave(-probs2$pScaled, probs2[, indIDVar2], 
                             FUN = function(x){
                               rank(x, ties.method = "min") 
                             })
  
  #Only keeping columns of interest
  probs2 <- probs2[, c("label", pairIDVar, "pAvg", "pSD", "pScaled", "pRank", "nEstimates")]
  
  
  #### Summarizing Measures of Effect Over iterations ####
  cAll_flat <- bind_rows(cAll) %>%
    group_by(level) %>%
    summarise(logorMean = mean(bs_est, na.rm = TRUE), 
              VW = mean(bs_sd ^ 2, na.rm = TRUE),
              VB = stats::var(bs_est, na.rm = TRUE)) %>%
    mutate(VT = VW + VB + VB/(n*nReps), 
           logorSE = sqrt(VT), 
           logorCILB = logorMean - 1.96*logorSE, 
           logorCIUB = logorMean + 1.96*logorSE)
  




  return(list("probabilities" = probs2, "estimates" = cAll_flat)) #, "var_probs" = pAll
}





