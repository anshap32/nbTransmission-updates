#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program adds covariates to a simulated outbreak
################################################################################


#### Description of Covariates ####

#X1 and X5 are correlated (corr = 0.8)
#For both (index cases): 60% chance of 1/a, 40% chance of 0/b
#Remaining individual level variables: 
# if infector for both is 1, 80% chance of 1 and 20% of 0
# if infector for both is 0, 80% chance of 0 and 20% of 1
# if infector has differing values for X1 and X2, 50% chance of 0 and 1



#X2 - categorical (a - 50%, b - 30%, c - 15%, d - 5%)
#Y2 - dichotomous (1 if match, 0 if not match)
#Linked pairs: 70% chance of matching, 30% change of not matching

#X3 - dichotomous (a - 70%, b - 30%)
#Y3 - categorical (1 = a-a, 2 = b-b, 3 = a-b, 4 = b-a)
#Linked pairs: if infector is a - 80% chance of a-a, 20% chance of a-b
#              if infector is b - 70% chance of b-b, 30% chance of b-a

#X4 - categorical (a - 5%, b - 5%, c - 5%, d - 20%, e - 30%,
#                  f - 10%, g - 10%, h - 5%, i - 5%, j - 5%)
#Y4 - categorical (1 if match, 2 if neighbors, 3 otherwise)
#Linked pairs: 60% chance of matching, 35% chance of neighbors, 5% chance of other

#X6 - categorical (a - 33%, b - 33%, c - 33%)
#Y6 - dichotomous (1 if match, 0 if not match)
#Linked pairs: 75% chance of matching, 25% change of not matching

#timeCat - categorical (1: <1y, 2: 1-2y, 3: 2-3y, 4: 3-4y, 5: 4-5y, 6: >5y)



simCovariates <- function(indData, pairData, observationDate = "infectionDate"){
  
  #Creating observation date variables
  indData$observationDate <- indData[, observationDate]
  pairData$observationDate.1 <- pairData[, paste0(observationDate, ".1")]
  pairData$observationDate.2 <- pairData[, paste0(observationDate, ".2")]
  
  
  #### Setting initial covariate information ####
  
  #Subseting to just the transmission pairs and removing index cases
  transPairs <- pairData %>% filter(transmission == TRUE, individualID.1 %% 1000 != 0)
  
  #Finding rows that need initial data - index cases and individuals who were not sampled
  origins <- indData %>% filter(is.na(sampleDate) | infector %% 1000 == 0)
  
  # simulate X1 and X3 jointly
  # Construct a binary correlation matrix
  m <- matrix(c(1,0.8,0.8,1), ncol=2)   
  
  ## Simulate covariate pairs for each of the index cases
  
  # set.seed(19)
  x <- rmvbin(nrow(origins), margprob = c(0.6, 0.6), bincorr = m)
  # cor(x)
  # table(x[,1], x[,2])
  
  #Set initial values for X1 and X5
  indData[indData$individualID %in% origins$individualID, "X1"] <- x[,1]
  indData[indData$individualID %in% origins$individualID, "X5"] <- x[,2]
  
  #Set initial values for X2
  X2info <- cbind.data.frame(value = c("a", "b", "c", "d"), freq = c(0.5, 0.3, 0.15, 0.05),
                             stringsAsFactors = FALSE)
  indData[indData$individualID %in% origins$individualID, "X2"] <- sample(X2info$value,
                                                                          size = nrow(origins),
                                                                          replace = TRUE, 
                                                                                                                                                  prob = X2info$freq)
  #Set initial values for X3
  X3info <- cbind.data.frame(value = c("a", "b"), freq = c(0.7, 0.3), stringsAsFactors = FALSE)
  indData[indData$individualID %in% origins$individualID, "X3"] <- sample(X3info$value,
                                                                          size = nrow(origins),
                                                                          replace = TRUE, 
                                                                          prob = X3info$freq)
  #Set initial values for X4
  X4info <- cbind.data.frame(value = letters[1:10],
                             freq = c(rep(0.05, 3), 0.2, 0.3, 0.1, 0.1, rep(0.05, 3)),
                             stringsAsFactors = FALSE)
  indData[indData$individualID %in% origins$individualID, "X4"] <- sample(X4info$value,
                                                                          size = nrow(origins),
                                                                          replace = TRUE, 
                                                                          prob = X4info$freq)
  #Set initial values for X6
  X6info <- cbind.data.frame(value = c("a", "b", "c"),
                             freq = c(0.33, 0.33, 0.33),
                             stringsAsFactors = FALSE)
  indData[indData$individualID %in% origins$individualID, "X6"] <- sample(X6info$value,
                                                                          size = nrow(origins),
                                                                          replace = TRUE, 
                                                                          prob = X6info$freq)
  
  
  #### Setting values for all other individuals ####
  
  #Loop through all the rows of the dataframe and assign covariate value of infectee based on infector
  for (i in 1:nrow(transPairs)){
    infector <- transPairs[i, "individualID.1"]
    infectee <- transPairs[i, "individualID.2"]
    
    X1val <- indData[indData$individualID == infector, "X1"]
    X5val <- indData[indData$individualID == infector, "X5"]
    
    if(X1val == X5val & X1val == 1){
      
      x <- rmvbin(1, margprob = c(0.8, 0.8), bincorr = m) 
      # cor(x); table(x[,1], x[,2]); table(x)
      
    } else if(X1val == X5val & X1val == 0){
      x <- rmvbin(1, margprob = c(0.2, 0.2), bincorr = m) 
      # cor(x);table(x[,1], x[,2]); table(x)
      
    }else {
      x <- rmvbin(1, margprob = c(0.5, 0.5), bincorr = m) 
    }
    
    indData[indData$individualID == infectee, "X1"] <- x[1]
    indData[indData$individualID == infectee, "X5"] <- x[2]
    
    X2val <- indData[indData$individualID == infector, "X2"]
    #Finding the values that are not the infector's
    other <- X2info$value[!X2info$value %in% X2val]
    #Summing the base frequenies for these options
    otherSum <- sum(X2info %>% filter(value %in% other) %>% select(freq))
    #Setting the new frequency to 0.7 for the infector's value and distributing the remaining 0.3
    #between the other options in proportion with their base frequencies
    X2info <- X2info %>% mutate(freqNew = ifelse(value == X2val, 0.7, (0.3 * freq) / otherSum))
    indData[indData$individualID == infectee, "X2"] <- sample(X2info$value, size = 1, prob = X2info$freqNew)
    
    X3val <- indData[indData$individualID == infector, "X3"]
    X3info <- X3info %>% mutate(freqNewA = c(0.8, 0.2),
                                freqNewB = c(0.7, 0.3))
    indData[indData$individualID == infectee, "X3"] <- ifelse(X3val == "a",
                                                              sample(X3info$value, size = 1, prob = X3info$freqNewA),
                                                              sample(X3info$value, size = 1, prob = X3info$freqNewB))
    
    
    #Finding the X4 value for the current individual
    X4val <- indData[indData$individualID == infector, "X4"]
    #Finding the index of that letter in the table with values and frequencies
    X4index <- which(X4info$value == X4val)
    #Finding the neighbors of that letter
    neighbors <- X4info$value[c(X4index - 1, X4index + 1)]
    #Finding the sum of the frequencies of the neighbors
    neighborSum <- sum(X4info %>% filter(value %in% neighbors) %>% select(freq))
    #Finding the other letters that aren't neighbors
    other <- X4info$value[!X4info$value %in% c(X4index, neighbors)]
    #Finding the sum of the frequencies of the other values
    otherSum <- sum(X4info %>% filter(value %in% other) %>% select(freq))
    
    #Setting new frequencies to be 60% for the same value, 35% split between neighbors at their
    #frequencies, and 5% split between other values at their frequencies
    X4info <- X4info %>% mutate(freqNew = ifelse(value == X4val, 0.6,
                                                 ifelse(value %in% neighbors, (0.35 * freq) / neighborSum,
                                                        (0.05 * freq) / otherSum)))
    indData[indData$individualID == infectee, "X4"] <- sample(X4info$value, size = 1, prob = X4info$freqNew)
    
    
    X6val <- indData[indData$individualID == infector, "X6"]
    #Finding the values that are not the infector's
    other <- X6info$value[!X6info$value %in% X6val]
    #Summing the base frequenies for these options
    otherSum <- sum(X6info %>% filter(value %in% other) %>% select(freq))
    #Setting the new frequency to 0.75 for the infector's value and distributing the remaining 0.25
    #between the other options in proportion with their base frequencies
    X6info <- X6info %>% mutate(freqNew = ifelse(value == X6val, 0.75, (0.25 * freq) / otherSum))
    indData[indData$individualID == infectee, "X6"] <- sample(X6info$value, size = 1, prob = X6info$freqNew)
  }
  
  
  
  #### Creating final datasets ####
  
  #For Y4, in order to find the neighboring values need to create numeric version of letters
  letterdf <- cbind.data.frame(letters, 1:26, stringsAsFactors = FALSE)
  names(letterdf) <- c("X4", "X4n")
  
  #Creating final individual level dataset
  #Creating a variable that denotes if the individual has complete info 
  covarInd <- (indData
               %>% mutate(complete = complete.cases(.), 
                          X1 = ifelse(X1 == 1, "a", "b"), 
                          X5 = ifelse(X5 == 1, "a", "b"))
               %>% left_join(letterdf, by = "X4")
  )
  
  #Merging individual-level covariates with the pairs dataset
  #Deriving pair-level covariates from the individual level
  #Setting reference groups for covariates
  #Creating a variable that denotes if the individual has complete info
  #Removing individual-level covariates
  covarPair <- (pairData
                %>% right_join(select(covarInd, individualID, X1, X2, X3, X4, X4n, X5, X6),
                               by = c("individualID.1" = "individualID"))
                %>% right_join(select(covarInd, individualID, X1, X2, X3, X4, X4n, X5, X6),
                               by = c("individualID.2" = "individualID"),
                               suffix = c(".1", ".2"))
                %>% mutate(Y1 = ifelse(X1.1 == X1.2, 1, 0),
                           Y2 = ifelse(X2.1 == X2.2, 1, 0),
                           Y3 = ifelse(X3.1 == "a" & X3.2 == "a", 1,
                                       ifelse(X3.1 == "b" & X3.2 == "b", 2,
                                              ifelse(X3.1 == "a" & X3.2 == "b", 3,
                                                     ifelse(X3.1 == "b" & X3.2 == "a", 4, NA)))),
                           Y4 = ifelse(X4.1 == X4.2, 1,
                                       ifelse(abs(X4n.1 - X4n.2) == 1, 2, 3)),
                           Y5 = ifelse(X5.1 == "a" & X5.2 == "a", 1,
                                       ifelse(X5.1 == "b" & X5.2 == "b", 2, 3)),
                           Y6 = ifelse(X6.1 == X6.2, 1, 0),
                           Y1 = factor(Y1, levels = c(0, 1)),
                           Y2 = factor(Y2, levels = c(0, 1)),
                           Y3 = factor(Y3, levels = c(1, 2, 3, 4)),
                           Y4 = factor(Y4, levels = c(1, 2, 3)),
                           Y5 = factor(Y5, levels = c(3, 1, 2)),
                           Y6 = factor(Y6, levels = c(0, 1)),
                           #Creating a categorical time difference variable
                           observationDiff = round(as.numeric(difftime(observationDate.2,
                                                                       observationDate.1, units = "days")), 0),
                           observationDiffY = observationDiff / 365,
                           observationDiffA = abs(observationDiff),
                           timeCat = ifelse(observationDiffA <= 365, 1,
                                            ifelse(observationDiffA > 365 & observationDiffA <= 730, 2,
                                                   ifelse(observationDiffA > 730 & observationDiffA <= 1095, 3,
                                                          ifelse(observationDiffA > 1095 & observationDiffA <= 1460, 4,
                                                                 ifelse(observationDiffA > 1460 & observationDiffA <= 1825, 5, 6))))),
                           timeCat = factor(timeCat, levels = c(1, 2, 3, 4, 5, 6)),
                           complete = complete.cases(.))
                %>% select(-X1.1, -X1.2, -X2.1, -X2.2, -X3.1, -X3.2, -X4.1, -X4.2,
                           -X4n.1, -X4n.2, -X5.1, -X5.2, -X6.1, -X6.2, -observationDiffA)
  )
  
  return(list(covarPair, covarInd))
}
