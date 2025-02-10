## Sarah's performNB function updated to include my weights
## called by runCV which is called by nbProbabilities

performNB_sBS <- function(training, prediction, obsIDVar, goldStdVar, 
                         covariates, l = 1, 
                         nBS = 100, pSampled = 1){
  
  #### Checking variable names ####
  
  #Checking that the named variables are in the dataframe
  if(!obsIDVar %in% names(training) | !obsIDVar %in% names(prediction)){
    stop(paste0(obsIDVar, " is not in the dataframe."))
  }
  if(!goldStdVar %in% names(training)){
    stop(paste0(goldStdVar, " is not in the dataframe."))
  }
  
  #Checking that the covariates are in the dataframe
  covarTestT <- c(covariates %in% names(training), covariates %in% names(prediction))
  if(FALSE %in% covarTestT){
    stop("At least one of the covariates is not in the input dataframes.")
  }
  
  #Checking that all of the covariates are factors
  covarDf <- training[, covariates]
  notFactor <- names(covarDf)[!sapply(covarDf, is.factor)]
  notFactorC <- paste0(notFactor, collapse = ", ")
  if(FALSE %in% sapply(covarDf, is.factor)){
    stop(paste0(notFactorC, " are not factors in the training data"))
  }
  covarDf <- prediction[, covariates]
  notFactor <- names(covarDf)[!sapply(covarDf, is.factor)]
  notFactorC <- paste0(notFactor, collapse = ", ")
  if(FALSE %in% sapply(covarDf, is.factor)){
    stop(paste0(notFactorC, " are not factors in the training data"))
  }
  
  
  #Making sure there are both events and non-event observations.
  #If not return NA for probabilities and print a warning
  if(sum(training[, goldStdVar] == TRUE, na.rm = TRUE) == 0 |
     sum(training[, goldStdVar] == FALSE, na.rm = TRUE) == 0){
    
    #Setting probability of a event to 0
    probs <- prediction
    probs$p <- NA
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- dplyr::bind_rows(probs[, c("p", obsIDVar)], training[, c("p", obsIDVar)])
    coeff <- NULL
    
    warning("No events or no non-events in training set")
    
  }
  
  else{
   
    #Setting up weight table
    #(currently irrelevant but would allow for weighting in subsequent versions)
    varTable <- as.data.frame(covariates)
    varTable$variable <- as.character(covariates)
    
    # ADD IN ATTRIBUTE WEIGHT
    varTable$weight <- 1
    
    #Creating the results dataframe which is a copy of the prediction dataframe
    results <- prediction
    
    #Finding proportion of events/non-events in the training dataset (prior -- P(c))
    classTab <- prop.table(table(training[, goldStdVar]) + l)
    
    #Initializing the coefficient dataframe
    coeff <- data.frame("level" = character(), "est" = numeric(),
                        "se" = numeric(), "Pr_F" = numeric(), "Pr_T" = numeric(), stringsAsFactors = FALSE)
    
    
    #Looping through all covariates and finding frequencies in training data
    for(i in 1:length(covariates)){
      
      #Extracting covariate name
      Var <- covariates[i]
      #Determining the weight for that covariate
      W <- varTable[varTable$variable == Var, "weight"]
      #Creating a table with proportions in each level of the covariate from training data
      Tab <- prop.table(W * table(training[, Var], training[, goldStdVar]) + l, 2)
      
      #Calculating P(X|Outcome = TRUE, W) and P(X|Outcome = FALSE, W) (P(a_j|c))
      #First creating a variable with the value of this covariate
      results$covariate <- results[, Var]
      #1 if missing, otherwise it is the proportion in the second column and the row
      #associated with the factor level of the covariate value.
      results[, paste0(Var, "_T")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 2]) ^ W
      results[, paste0(Var, "_F")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 1]) ^ W
      
      # #Saving the odds: P(Outcome = TRUE|X=x, W) / P(Outcome = FALSE|X=x, W) -- THIS BECOMES THE SAME THING AS P(X=x|Outcome = TRUE) / P(X=x|Outcome = FALSE)
      # odds <- Tab[, 2] / Tab[, 1]
      # est <- log(odds/odds[1]) #P(X=x|Outcome = TRUE) / P(X=x|Outcome = FALSE)
      # num <- W * table(training[, Var], training[, goldStdVar]) + l
      # se <- NA
      # for(i in 1:(nrow(num) - 1)){
      #   numTab <- num[c(1, i+1), ]
      #   se <- c(se, sqrt(sum(1/numTab)))
      
      
      # level <- paste(Var, names(est), sep = ":")
       
      # should be a logistic regression model where we save the model fit to a list
      
      # MH_Tab <- prop.table( table(training[, Var], training[, goldStdVar]) + l, 1)
      # Pr_F <- MH_Tab[, 1]
      # Pr_T <- MH_Tab[, 2]
      
      # cTemp <- cbind.data.frame(level, est, se, Pr_F, Pr_T, stringsAsFactors = FALSE)
      # coeff <- mylogit
    }
    
    ### ADDED BY ME: LOGISTIC REGRESSION
    
    trainingT <- training %>% filter(linked == T)
    trainingF <- training %>% filter(linked == F)
    
    bs_out <- list()
    
    for (i in 1:nBS) {
      nSampled <- nrow(trainingF)*pSampled
      trainingF_samp <- sample_n(trainingF, nSampled, replace = T)
      bs_samp <- rbind(trainingF_samp, trainingT)
      
      mylogit <- glm(reformulate(covariates, goldStdVar), data = bs_samp, family = "binomial")
      tidy_fit <- tidy(mylogit)
      tidy_out <- data.frame(tidy_fit$term, tidy_fit$estimate)
      
      
      bs_out[[i]] <- tidy_out
    }
    
    bs_out1 <- bind_rows(bs_out)
    
  bs_out_df <- bs_out1 %>%
    group_by(tidy_fit.term) %>%
    summarise(bs_est = mean(tidy_fit.estimate), 
              bs_sd = sd(tidy_fit.estimate)) %>%
    rename(level = tidy_fit.term)
  
  #Calculating numerator and denominator for the probability calculation
    results$event <- apply(results[, grepl("_T", names(results))], 1, prod) * classTab[2] 
    results$nonevent <- apply(results[, grepl("_F", names(results))], 1, prod) * classTab[1]
    
    #Calculating probability of event
    probs0 <- results
    probs0$p <- probs0$event / (probs0$event + probs0$nonevent)
    
    #Combining training and prediction datasets
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- dplyr::bind_rows(probs0[, c(obsIDVar, "p")], training[, c(obsIDVar, "p")])
  }
  
  return(list("probabilities" = probs, "estimates" = bs_out_df))
}
