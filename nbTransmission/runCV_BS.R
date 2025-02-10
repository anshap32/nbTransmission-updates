## Sarah's runCV function updated to include my weights
## calls performNB_LR, called by nbProbabilities

runCV_sBS <- function(posTrain, orderedPair, indIDVar, pairIDVar,
                     goldStdVar, covariates, l, n, m, 
                     nBS = 100, pSampled = 1
                     ){
  
  #Creating variables with the individual indID variable
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  if(length(unique(posTrain$pairID_uo)) != nrow(posTrain)){ # this is for people we are unsure of infector 
    
    #Finding the pairs that are concurrent (both orders included)
    concurrentID <- posTrain[duplicated(posTrain$pairID_uo), "pairID_uo"]
    concurrent <- posTrain[posTrain$pairID_uo %in% concurrentID, ]
    
    #Randomly selecting order for concurrent pairs
    concurrent2L <- linksL <- by(concurrent,
                                 INDICES = list(concurrent$pairID_uo),
                                 FUN = function(x){
                                   x[sample(nrow(x), 1), ]
                                 })
    concurrent2 <- do.call(dplyr::bind_rows, concurrent2L)
    
    #Combining choose for concurrent pairs with the rest of the pairs
    posTrain2 <- dplyr::bind_rows(posTrain[!posTrain$pairID_uo %in% concurrentID, ],
                                  concurrent2)
  }else{
    posTrain2 <- posTrain
  }
  
  #Subsetting to just the possible links
  posLinks <- posTrain2[posTrain2[, goldStdVar] == TRUE, ]
  # rm(posLinks2)
  
  #Choosing the true infector from all possibles (if multiple) (I don't think we have this in simulated data from example)
  # set.seed(k)
  linksL <- by(posLinks,
               INDICES = list(posLinks[, indIDVar2]),
               FUN = function(x){
                 x[sample(nrow(x), 1), ]
               })
  links <- do.call(dplyr::bind_rows, linksL)
  links$linked <- TRUE
  links <- links[, c(pairIDVar, indIDVar2, "linked")]
  
  #Combining the links with the non-links that do not share an infectee with the links
  trainingFull <- dplyr::full_join(orderedPair, links, by = c(pairIDVar, indIDVar2))
  trainingFull2 <- trainingFull[trainingFull[, pairIDVar] %in% links[, pairIDVar] |
                                  (trainingFull[, pairIDVar] %in% posTrain2[, pairIDVar] &
                                     !trainingFull[, indIDVar2] %in% links[, indIDVar2]), ]
  trainingFull2[is.na(trainingFull2$linked), "linked"] <- FALSE
  # rm(trainingFull)
  
  
  #Creating the cross-valindIDation folds for that part of the training dataset
  # set.seed(k)
  cv_splits <- caret::createMultiFolds(trainingFull2$linked, k = n, times = m)
  
  #Initializing data frames to hold results and coefficients
  # rFolds <- data.frame("p" = numeric(), pairIDVar = character(), stringsAsFactors = FALSE)
  # names(rFolds) <- c("p", pairIDVar)
  # cFolds <- data.frame("level" = character(), "est" = numeric(), stringsAsFactors = FALSE)
  
  rFolds <- list()
  cFolds <- list()
  
  #Running the methods for all of the CV Folds
  for (i in 1:length(cv_splits)){
    
    #Finding training dataset
    trainingPairID <- trainingFull2[, pairIDVar][cv_splits[[i]]]
    trainingRaw <- trainingFull2[trainingFull2[, pairIDVar] %in% trainingPairID, ]
    
    #Finding the infectee whose infector was found in the training dataset
    foundInfector <- trainingRaw[trainingRaw$linked == TRUE, ]
    #Finding all pairs that share an infectee with the pairs where the infector was found
    shareInfectee <- orderedPair[!orderedPair[, pairIDVar] %in% foundInfector[, pairIDVar] &
                                   orderedPair[, indIDVar2] %in% foundInfector[, indIDVar2], ]
    shareInfectee$linked <- FALSE
    
    #Creating the training datasets
    #Setting probabilities to 1 for training links and 0 for training non-links that are
    #defined as such by the gold standard (not just by sharing an infectee in the above df)
    training <- dplyr::bind_rows(trainingRaw, shareInfectee)
    training$p <- ifelse(training[, goldStdVar] == FALSE, 0,
                         ifelse(training[, "linked"] == TRUE, 1, NA))
    
    #Creating the prediction dataset
    prediction <- dplyr::full_join(orderedPair, links, by = c(pairIDVar, indIDVar2))
    prediction <- prediction[!prediction[, pairIDVar] %in% training[, pairIDVar], ]
    prediction[is.na(prediction$linked), "linked"] <- FALSE
    
    #Calculating probabilities for one split
    sim <- performNB_sBS(training, prediction, obsIDVar = pairIDVar,
                        goldStdVar = "linked", covariates, l, 
                        nBS, pSampled) #
    
    #Combining the results from fold iteration with the previous folds
    # rFolds <- dplyr::bind_rows(rFolds, sim[[1]])
    # sim[[2]]$rowOrder <- 1:nrow(sim[[2]])
    # cFolds <- dplyr::bind_rows(cFolds, sim[[2]])
    
    rFolds[[i]] <- sim[[1]]
    # sim[[2]]$rowOrder <- 1:nrow(sim[[2]])
    cFolds[[i]] <- sim[[2]]
    cFolds[[i]]$fold <- i
  }
  rFolds <- rFolds %>% bind_rows()
  cFolds <- cFolds %>% bind_rows()
  return(list("rFolds" = rFolds, "cFolds" = cFolds)) #, "pFolds" = pFolds
}
