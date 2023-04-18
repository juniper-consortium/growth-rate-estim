
# Functions for running model ----

# Auxiliar function
#' Loads parameters
#' linkType:
#' prior2.sigma0: base-line value for 1 / precision of covariance function
#' prior2.lengthscale0: base-line value for the length scale of covariance function
#' theta.prior2.mean: mean of theta prior
#' theta.prior2.prec: precision matrix of theta prior
#' BB.prior.rho: prior for log precision if beta-binomial case
#' NB.prior.rho: prior for log precision if negative binomial case
#' dayWeek: if T, the day-of-the-week effect is used. If F, the noise term is used
#' dayWeek.prior.prec: prior for day-of-the-week effect or noise term
#' sizeSample: samples to estimate posterior of parameters
setParametersFn <- function(linkType,
                            prior2.sigma0 = 1,
                            prior2.lengthscale0 = 100,
                            theta.prior2.mean = c(0,0),
                            theta.prior2.prec = diag(2),
                            BB.prior.rho = list(overdispersion = list(prior = "gaussian", param = c(0, 0.5))),
                            NB.prior.rho = list(size = list(prior = "gaussian", param = c(0, 0.01))),
                            unitTime = "day",
                            dayWeek = unitTime == "day",
                            dayWeek.prior.prec = list(theta = list(prior = 'loggamma', param = c(1, 0.01))),
                            sizeSample = sizeSample){
  parameters <- list(
    params = list(
      linkType = linkType,
      prior2.sigma0 = prior2.sigma0,
      prior2.range0 = 2*prior2.lengthscale0,
      theta.prior2.mean = theta.prior2.mean,
      theta.prior2.prec = theta.prior2.prec,
      BB.prior.rho = BB.prior.rho,
      NB.prior.rho = NB.prior.rho,
      unitTime = unitTime,
      dayWeek = dayWeek,
      dayWeek.prior.prec = dayWeek.prior.prec
    ),
    config = list(sizeSample = sizeSample))
    return(parameters)
}

#' countTable: data table with one row per day (of days with available data) and these columns:
#' - numberTest: number of test on the day (>= 0)
#' - positiveResults: number of positive tests (<= numberTest and >= 0)
#' - date: date of count in R date format
#' parametersModel: output of setParametersFn()
#' saveSamples: Default F. If T, returns matrixSampleDays and sampleDerivatives,
#'              two matrices of size [days, num. samples] containing samples of the posterior of the GP and GP derivative respectively.
#' minDate (optional): minimum date to include in the model
#' maxDate (optional): maximum date to include in the model
#' derivativeFromGP (optional): if growth rate distribution is estimated from samples of GP derivative. If F, dataset must have at least 7 time points between minDate and maxDate.
runModelGrowthRate <- function(countTable, parametersModel, saveSamples = F, minDate = NULL, maxDate = NULL, derivativeFromGP = F){
  # Create auxiliar variables
  levelsWeek <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  
  # Create auxiliar table with dates
  if(is.null(minDate)) minDate <- min(countTable$date)
  if(is.null(maxDate)) maxDate <- max(countTable$date)
  minDay <- 1
  maxDay <- length(seq.Date(from = minDate, to = maxDate, by = parametersModel$params$unitTime)) #as.integer(maxDate - minDate + 1)
  numDays <- maxDay #maxDay - minDay + 1
  dateTable <- data.table(dayId = minDay:maxDay,
                          date = seq.Date(from = minDate, to = maxDate, by = parametersModel$params$unitTime))
  #if(nrow(dateTable) < 7) stop("There must be at least 7 time units between minDate and maxDate")
  if(nrow(dateTable) <= 1) stop("There must be at least 2 time units between minDate and maxDate")
  
  # ---------------------------------------------------- #
  #                      FIT MODEL                       #
  # ---------------------------------------------------- #
  # Create data with all days, including the ones with missing data
  dataForModel <- data.table(dayId = minDay:maxDay,
                             date = seq.Date(from = minDate, to = maxDate, by = parametersModel$params$unitTime),
                             numberTest = as.integer(NA),
                             positiveResults = as.integer(NA))
  setkey(dataForModel, date)
  setkey(countTable, date)
  dataForModel[countTable, ":="(numberTest = i.numberTest, positiveResults = i.positiveResults)]
  dataForModel[numberTest == 0, ":="(numberTest = NA, positiveResults = NA)]
  dataForModel[, dayWeek := weekdays(date)]
  
  # Define grid for finite differences
  boundaryVertices <- 200 # Points on the left and right
  boundaryStepSize <- 10 # Distance between boundary points
  boundaryPoints <- c(seq(min(dataForModel$dayId) - boundaryStepSize*boundaryVertices, min(dataForModel$dayId) - boundaryStepSize, by = boundaryStepSize),
                      seq(min(dataForModel$dayId), max(dataForModel$dayId), by = 1), 
                      seq(max(dataForModel$dayId) + boundaryStepSize, max(dataForModel$dayId) + boundaryStepSize*boundaryVertices, by = boundaryStepSize))
  nonBoundaryIndices <- (boundaryVertices + 1):(length(boundaryPoints) - boundaryVertices)
  
  # Create mesh in 1d (mesh1d), projection matrix (A1), model (spde1), and indexes (spde1.idx)
  
  # - Mesh and projection matrix
  mesh1d <- inla.mesh.1d(boundaryPoints, boundary = "free") # NEW 09.02.2021 boundary issues
  A1 <- inla.spde.make.A(mesh1d, dataForModel[order(dayId), dayId])
  
  # - Set priors
  # Prior for day of the week
  priorDayWeek <- parametersModel$params$dayWeek.prior.prec
  # Prior for Gaussian process parameters, following (Lindgren, 2015, v63i19)
  vGP <- 2 - 1/2
  sigma0 <- parametersModel$params$prior2.sigma0
  range0 <- parametersModel$params$prior2.range0
  # Convert into tau and kappa:
  kappa0 <- sqrt(8*vGP)/range0
  tau0 <- sqrt(gamma(vGP)/(gamma(2)*sqrt(4*pi)))/(kappa0^(vGP)*sigma0)
  basis.prior2.tau <- c(-1, vGP)
  basis.prior2.kappa <- c(0, -1)
  spde1 <- inla.spde2.matern(mesh1d,
                             B.tau = cbind(log(tau0), basis.prior2.tau[1], basis.prior2.tau[2]),
                             B.kappa = cbind(log(kappa0), basis.prior2.kappa[1], basis.prior2.kappa[2]),
                             theta.prior.mean = parametersModel$params$theta.prior2.mean,
                             theta.prior.prec = parametersModel$params$theta.prior2.prec)
  
  # - Create stack data
  spde1.idx <- inla.spde.make.index("day", n.spde = spde1$n.spde)
  # Create stacked data with two datasets (observations and prediction). For paper, I'm ignoring prediction.
  stackD <- inla.stack(data = list(positiveResults = dataForModel[order(dayId), positiveResults]),
                       A = list(A1, 1, 1),
                       effects = list(c(list(Intercept = 1), spde1.idx),
                                      list(dayWeek = factor(dataForModel[order(dayId), dayWeek], levels = levelsWeek)),
                                      list(dayId2 = dataForModel[order(dayId), dayId])),
                       tag = "est")
  # Create prediction data for the next sizeGPPred days
  sizeGPPred <- 10
  predictionGPWeekday <- dateTable[dayId == max(dataForModel$dayId), weekdays(date + 1:sizeGPPred)]
  xx <- seq(max(dataForModel$dayId) + 1, max(dataForModel$dayId) + sizeGPPred, by = 1)
  A.xx <- inla.spde.make.A(mesh1d, xx)
  stack.pred <- inla.stack(data = list(positiveResults = NA),
                           A = list(A.xx),
                           effects = list(c(list(Intercept = 1), spde1.idx)),
                           tag = "pred")
  joint.stack <- inla.stack(stackD, stack.pred)
  
  # Fit model (using INLA)
  cat("Fitting model ... ")
  if(parametersModel$params$dayWeek == T){
    formula <- positiveResults ~ -1 + f(day, model = spde1) + f(dayWeek, model = 'iid', hyper = priorDayWeek, constr = T)
  }else{
    formula <- positiveResults ~ -1 + f(day, model = spde1) + f(dayId2, model = 'iid', hyper = priorDayWeek, constr = T)
  }
  if(parametersModel$params$linkType == "BB"){
    # Betabinomial
    m.spde <- inla(formula = formula, data = inla.stack.data(joint.stack),
                   control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                   family = "betabinomial",
                   Ntrials =  c(dataForModel[order(dayId), numberTest], rep(NA, length(xx))),
                   control.compute = list(config = TRUE),
                   control.family = list(hyper = parametersModel$params$BB.prior.rho))
  }else if(parametersModel$params$linkType == "NB"){
    # Negative binomial
    m.spde <- inla(formula = formula, data = inla.stack.data(joint.stack),
                   control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                   family = "nbinomial",
                   control.compute = list(config = TRUE),
                   control.family = list(hyper = parametersModel$params$NB.prior.rho))
  }
  
  # ---------------------------------------------------- #
  #                      GET SAMPLES                     #
  # ---------------------------------------------------- #
  cat("Estimating growth rate ... ")
  
  # Get samples of posterior distribution of parameters
  sampleList <- inla.posterior.sample(parametersModel$config$sizeSample, m.spde)
  
  # Extract relevant indexes from output
  rownms <- rownames(sampleList[[1]]$latent)
  set3 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "day"))
  dayIndexInSample <- set3[nonBoundaryIndices]
  
  # Create matrix of samples (and hyperparameters if needed)
  matrixSampleDays <- sapply(sampleList, function(x) x$latent[dayIndexInSample,1])
  # matrixSampleHyper <- sapply(sampleList, function(x) exp(x$hyperpar[c("Theta1 for day", "Theta2 for day")]))
  matrixSampleHyper <- sapply(sampleList, function(x) exp(x$hyperpar)) # !!!
  
  # Compute derivative
  if(derivativeFromGP == T | nrow(dateTable) < 7){
    # Sample from derivative
    sampleDerivatives <- getGrowthFromSamples_GP(matrixSampleGP = matrixSampleDays, samplesHyperparam = matrixSampleHyper, sigma0 = sigma0, range0 = range0)
  }else{
    # Compute approximate derivative using windows (+-3)
    # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
    sampleDerivatives <- getGrowthFromSamples(matrixSampleDays = matrixSampleDays)
  }
  
  # ---------------------------------------------------- #
  #                  PRODUCE OUTPUT                      #
  # ---------------------------------------------------- #
  listPosteriors <- computePosteriors(matrixSampleDays, sampleDerivatives, parametersModel)
  
  setkey(listPosteriors$posteriorGrowth, dayId)
  setkey(dateTable, dayId)
  listPosteriors$posteriorGrowth[dateTable, ":="(date = i.date)]
  
  setkey(listPosteriors$posteriorTransfGP, dayId)
  setkey(dateTable, dayId)
  listPosteriors$posteriorTransfGP[dateTable, ":="(date = i.date)]
  setkey(listPosteriors$posteriorTransfGP, date)
  setkey(countTable, date)
  listPosteriors$posteriorTransfGP[countTable, ":="(positiveResults = i.positiveResults)]
  
  if(saveSamples == F){
    return(list(posteriorGrowth = listPosteriors$posteriorGrowth, posteriorTransfGP = listPosteriors$posteriorTransfGP))
  }else{
    return(list(posteriorGrowth = listPosteriors$posteriorGrowth, posteriorTransfGP = listPosteriors$posteriorTransfGP,
                matrixSampleDays = matrixSampleDays, sampleDerivatives = t(sampleDerivatives), matrixSampleHyper = matrixSampleHyper,
                summaryHyperpar = data.table(m.spde$summary.hyperpar, keep.rownames = T),
                m.spde = m.spde)) # TODO ?
  }
  
}

# Auxiliar functions for running model ----

# INTERNAL
computePosteriors <- function(matrixSampleDays, sampleDerivatives, parametersModel, ifINLAMarginal = FALSE, m.spde = NULL){
  numDays <- nrow(matrixSampleDays)
  
  # ---------------------------------------------------- #
  #                POSTERIOR GROWTH RATE                 #
  # ---------------------------------------------------- #
  # Compute log (or inv. logit) of posterior of Gaussian process derivative in transformed space
  cat("Computing posterior of growth rate ... ")
  if(parametersModel$params$linkType == "BB"){
    tempExpGP <- exp(matrixSampleDays)
    tempLogit <- tempExpGP/(1 + tempExpGP)
    tempList <- sampleDerivatives/(1 + tempLogit)
  }else if(parametersModel$params$linkType == "NB"){
    tempList <- sampleDerivatives
  }
  tempDoubling <- abs(log(2)/tempList)
  posteriorGrowth <- data.table(dayId = 1:numDays,
                                mean = rowMeans(tempList),
                                sd = apply(tempList, 1, sd),
                                median = apply(tempList, 1, quantile, probs = 0.5, na.rm = T),
                                q0.025 = apply(tempList, 1, quantile, probs = 0.025, na.rm = T),
                                q0.975 = apply(tempList, 1, quantile, probs = 0.975, na.rm = T),
                                q0.25 = apply(tempList, 1, quantile, probs = 0.250, na.rm = T),
                                q0.75 = apply(tempList, 1, quantile, probs = 0.750, na.rm = T),
                                prob0 = apply(tempList, 1, function(x) sum(x >= 0)/parametersModel$config$sizeSample),
                                medianDoubT = apply(tempDoubling, 1, quantile, probs = 0.5, na.rm = T),
                                q0.025DoubT = apply(tempDoubling, 1, quantile, probs = 0.025, na.rm = T),
                                q0.975DoubT = apply(tempDoubling, 1, quantile, probs = 0.975, na.rm = T))
  
  # ---------------------------------------------------- #
  #                 POSTERIOR INCIDENCE                  #
  # ---------------------------------------------------- #
  # Compute posterior of Gaussian process in real space (incidence or positivity)
  cat("Computing posterior of incidence... ")
  if(ifINLAMarginal){
    # (this version is slow and applies only to INLA)
    samplesGP <- m.spde$marginals.random$day[nonBoundaryIndices]
    if(parametersModel$params$linkType %in% c("NB")){
      transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(exp, x))
    }else if(parametersModel$params$linkType == "BB"){
      transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(function(x) exp(x)/(1 + exp(x)), x))
    }
    posteriorTransfGP <- data.table(dayId = 1:numDays,
                                    t(sapply(transformedSamples, function(x) inla.qmarginal(c(0.5, 0.025, 0.975, 0.25, 0.75), x))))
    setnames(posteriorTransfGP, c("V1", "V2", "V3", "V4", "V5"), c("median", "q0.025", "q0.975", "q0.25", "q0.75"))
  }else{
    # NEW 09.01.2023
    if(parametersModel$params$linkType %in% c("NB")){
      transformedSamples <- exp(matrixSampleDays)
    }else if(parametersModel$params$linkType == "BB"){
      transformedSamples <- exp(matrixSampleDays)/(1 + exp(matrixSampleDays))
    }
    posteriorTransfGP <- data.table(dayId = 1:numDays,
                                  median = apply(transformedSamples, 1, quantile, probs = 0.5, na.rm = T),
                                  q0.025 = apply(transformedSamples, 1, quantile, probs = 0.025, na.rm = T),
                                  q0.975 = apply(transformedSamples, 1, quantile, probs = 0.975, na.rm = T),
                                  q0.25 = apply(transformedSamples, 1, quantile, probs = 0.250, na.rm = T),
                                  q0.75 = apply(transformedSamples, 1, quantile, probs = 0.750, na.rm = T))
    
  }
  
  return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP))
}

# INTERNAL
# NOT AVAILABLE ANYMORE as I don't trust INLA sampling of hyperparameters
getGrowthFromSamples <- function(matrixSampleDays){
  # Compute approximate derivative using windows (+-3)
  # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
  numDays <- nrow(matrixSampleDays)
  sampleDerivatives <- rbind(matrix(NA, nrow = 1, ncol = ncol(matrixSampleDays)),
                             ((matrixSampleDays[3,] - matrixSampleDays[1,])/2),
                             ((matrixSampleDays[5,] - matrixSampleDays[1,])/4),
                             (matrixSampleDays[7:numDays,] - matrixSampleDays[1:(numDays - 7 + 1),])/6,
                             ((matrixSampleDays[numDays,] - matrixSampleDays[numDays - 5 + 1,])/4),
                             ((matrixSampleDays[numDays,] - matrixSampleDays[numDays - 3 + 1,])/2),
                             matrix(NA, nrow = 1, ncol = ncol(matrixSampleDays)))
  return(sampleDerivatives)
}

# INTERNAL
#' matrixSampleDays: days x samples
#' samplesHyperparam: (exp(theta 1), exp(theta2)) X samples. sigma = sigma0 exp(theta1), range = range0 exp(theta2)
#' Copied FROM R47/Sandbox_deerivative.R
getGrowthFromSamples_GP <- function(matrixSampleGP, samplesHyperparam, sigma0, range0){
  numDays <- nrow(matrixSampleGP)
  sizeSample <- ncol(matrixSampleGP)
  
  # Compute distance matrix for ordered days
  distanceMatrix <- sapply(1:numDays, function(nd) abs(nd - (1:numDays)))
  auxRelativeDistanceMatrix <- matrix(data = 1:numDays, nrow = numDays, ncol = numDays, byrow = F) - matrix(data = 1:numDays, nrow = numDays, ncol = numDays, byrow = T)
  
  # Compute auxiliar vectors
  sig2Vector <- (sigma0*samplesHyperparam["Theta1 for day",])^2 # 1/samplesHyper["tau", indexSample]
  kappaVector <- sqrt(12)/(range0*samplesHyperparam["Theta2 for day",]) # samplesHyper["kappa", indexSample]
  
  # Loop per sample of (f1, ..., fn, log.tau, log.kappa)
  sampleDerivatives <- matrix(0, nrow = sizeSample, ncol = numDays)
  for(indexSample in 1:sizeSample){
    sig2Value <- sig2Vector[indexSample]
    kappaVal <- kappaVector[indexSample]
    
    expMatrix <- exp(-kappaVal*distanceMatrix)
    deltaMatrix <- sig2Value*(1 + kappaVal*distanceMatrix)*expMatrix
    invDeltaMatrix <- chol2inv(chol(deltaMatrix)) # solve vs. chol2inv system.time(31700*system.time(solve(deltaMatrix))/60)
    fVector <- matrixSampleGP[, indexSample]
    
    # Compute derivative matrices of f: D1, diag(D) in notes respectively
    d1Matrix <- -sig2Value*kappaVal^2*auxRelativeDistanceMatrix*expMatrix
    #dMatrix <- sig2Value*kappaVal^2*diag(expMatrix)*(1 - kappaVal*diag(distanceMatrix)) # here we only compute the diagonal
    dMatrixAll <- sig2Value*kappaVal^2*expMatrix*(1 - kappaVal*distanceMatrix)
    
    meanMVN <- d1Matrix%*%invDeltaMatrix%*%fVector
    
    #iSample <- sapply(1:numDays, function(x) rnorm(n = 1, mean = meanMVN[x], sd = sqrt( dMatrix[x] - d1Matrix[x,]%*%invDeltaMatrix%*%d1Matrix[x,] )))
    iSample <- MASS::mvrnorm(n = 1, mu = d1Matrix%*%invDeltaMatrix%*%fVector, Sigma = dMatrixAll - d1Matrix%*%invDeltaMatrix%*%t(d1Matrix))
    
    sampleDerivatives[indexSample,] <- iSample
  }
  return(t(sampleDerivatives))
}

# Functions for running multiple groups ----

#' Wrapper for runModelGrowthRate(), to run the model in multiple /groupsareas/partitions/locations if same settings/priors
#' countTableAll: data table with one row per day per group and these columns:
#' - labelPartition: name of group or 'partition'
#' - positiveResults: number of positive tests (<= numberTest and >= 0)
#' - numberTest: number of test on the day (>= 0)
#' - date: date of count in R date format
#' partitionTable
#' - labelPartition: name of group or 'partition'
#' - idPartition: unique identifier of each 'partition', numeric
#' - minDate: min date for which the model is going to be fitted for that partition
#' - maxDate: max date for which the model is going to be fitted for that partition
#' parametersModel: output of setParametersFn()
#' saveSamples: Default F. If T, returns matrixSampleDays_list and sampleDerivatives_list,
#'              a list of matrices of size [days, num. samples] containing samples of the posterior of the GP and GP derivative respectively.
#' runSubsetId: vector ids as in partitionTable$idPartition. Default: all ids ordered
runModelMultipleGroups <- function(countTableAll, partitionTable, parametersModel, saveSamples = F, derivativeFromGP = F,
                                   runSubsetId = sort(partitionTable$idPartition)){
  # TODO note that try catch does not exist for unique runs
  # TODO replace "saveSamples" and other inputs for "..."
  # Run model for all areas
  output <- vector("list", length(runSubsetId))
  for(ii in 1:length(runSubsetId)){
    tryCatch({
      labelPartition_ii <- partitionTable[idPartition == runSubsetId[ii], labelPartition]
      minDateModel_ii <- partitionTable[idPartition == runSubsetId[ii], minDate]
      maxDateModel_ii <- partitionTable[idPartition == runSubsetId[ii], maxDate]
      cat(labelPartition_ii, ": ", sep = "")
      
      countTable <- countTableAll[labelPartition == labelPartition_ii & date >= minDateModel_ii & date <= maxDateModel_ii,
                                  .(date, positiveResults, numberTest)]
      output[[ii]] <- runModelGrowthRate(countTable = countTable,
                                         parametersModel = parametersModel, saveSamples = saveSamples, derivativeFromGP = derivativeFromGP,
                                         minDate = minDateModel_ii, maxDate = maxDateModel_ii)
      cat("\n")
    }, error = function(e) {cat("SKIPPED ERROR for labelPartition", labelPartition_ii, ":", conditionMessage(e), "\n")})
  }
  return(output)
}

#' Stack output of runModelMultipleGroups() into a data table.
#' Intput:
#' - output: output of runModelMultipleGroups()
#' Output:
#' - posteriorGrowth
stackOutputMultipleGroups <- function(output, partitionTable){
  # Stack output
  posteriorGrowth <- do.call("rbind", lapply(1:length(output), function(ii) output[[ii]]$posteriorGrowth[, c(.SD, idPartition = ii)]))
  posteriorTransfGP <- do.call("rbind", lapply(1:length(output), function(ii) output[[ii]]$posteriorTransfGP[, c(.SD, idPartition = ii)]))
  
  #if(saveSamples == T){
  #  matrixSampleDays_list <- lapply(1:length(output), function(ii) output[[ii]]$matrixSampleDays)
  #  sampleDerivatives_list <- lapply(1:length(output), function(ii) output[[ii]]$sampleDerivatives)
  #}
  
  setkeyv(posteriorGrowth, c("date", "idPartition"))
  setkeyv(posteriorTransfGP, c("date", "idPartition"))
  posteriorGrowth[posteriorTransfGP, ":="(median_incidence = i.median, q0.025_incidence = i.q0.025, q0.975_incidence = i.q0.975,
                                          q0.25_incidence = i.q0.25, q0.75_incidence = i.q0.75, positiveResults = i.positiveResults)]
  
  setkey(posteriorGrowth, idPartition)
  setkey(partitionTable, idPartition)
  posteriorGrowth[partitionTable, labelPartition := i.labelPartition]
  
  # Calculate probability growth rate is greater than 0. Add to posteriorGrowth
  #hasSamples <- sapply(1:length(output), function(i) !is.null(output[[i]]$matrixSampleDays))
  #if(sum(hasSamples > 0)){
  #  prob0Table <- do.call("rbind",
  #                        lapply(which(hasSamples), function(i) data.table(idPartition = i,
  #                                                                         dayId = 1:nrow((output[[i]]$sampleDerivatives)),
  #                                                                         prob0 = apply(output[[i]]$sampleDerivatives, 1, function(x) sum(x >= 0)/length(x)))))
  #  setkeyv(posteriorGrowth, c("idPartition", "dayId"))
  #  setkeyv(prob0Table, c("idPartition", "dayId"))
  #  posteriorGrowth[prob0Table, prob0 := i.prob0]
  #}
  
  return(posteriorGrowth)
  
  # Add extras to output
  #posteriorGrowth[, dayId := NULL]
  #setkey(posteriorGrowth, idPartition)
  #setkey(partitionTable, idPartition)
  #posteriorGrowth[partitionTable, ":="(labelPartition = i.labelPartition, pop = i.pop)]
  #posteriorGrowth[partitionTable, ":="(labelPartition = i.labelPartition)]
  #
  #posteriorGrowth[, ":="(median_incidencePer100TH = median_incidence*100000/pop,
  #                       q0.025_incidencePer100TH = q0.025_incidence*100000/pop,
  #                       q0.975_incidencePer100TH = q0.975_incidence*100000/pop)]
  #posteriorGrowth[, positivePer100TH := positiveResults*100000/pop]
  
  #setkey(posteriorTransfGP, idPartition)
  #setkey(partitionTable, idPartition)
  #posteriorTransfGP[partitionTable, ":="(labelPartition = i.labelPartition)]
  #posteriorTransfGP[, ratio := NA]
  
  #if(saveSamples == F){
  #return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP))
  #}else{
  # TODO
  #return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP, matrixSampleDays = matrixSampleDays, sampleDerivatives = t(sampleDerivatives)))
  #}
}

