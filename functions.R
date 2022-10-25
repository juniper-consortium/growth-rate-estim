
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
                            dayWeek = T,
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
runModelGrowthRate <- function(countTable, parametersModel, saveSamples = F, minDate = NULL, maxDate = NULL){
  # Load parameters into function environment and add auxiliar variables
  list2env(parametersModel$params, envir = environment())
  list2env(parametersModel$config, envir = environment())
  levelsWeek <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  
  # Create auxiliar table with dates
  if(is.null(minDate)) minDate <- min(countTable$date)
  if(is.null(maxDate)) maxDate <- max(countTable$date)
  minDay <- 1
  maxDay <- as.integer(maxDate - minDate + 1)
  numDays <- maxDay - minDay + 1
  dateTable <- data.table(dayId = minDay:maxDay,
                          date = seq.Date(from = minDate, to = maxDate, by = "day"))
  
  # ---------------------------------------------------- #
  #                      FIT MODEL                       #
  # ---------------------------------------------------- #
  # Create data with all days, including the ones with missing data
  dataForModel <- data.table(dayId = minDay:maxDay,
                             date = seq.Date(from = min(countTable$date), to = max(countTable$date), by = "day"),
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
  priorDayWeek <- dayWeek.prior.prec
  # Prior for Gaussian process parameters, following (Lindgren, 2015, v63i19)
  vGP <- 2 - 1/2
  sigma0 <- prior2.sigma0
  range0 <- prior2.range0
  # Convert into tau and kappa:
  kappa0 <- sqrt(8*vGP)/range0
  tau0 <- sqrt(gamma(vGP)/(gamma(2)*sqrt(4*pi)))/(kappa0^(vGP)*sigma0)
  basis.prior2.tau <- c(-1, 3/2)
  basis.prior2.kappa <- c(0, -1)
  spde1 <- inla.spde2.matern(mesh1d,
                             B.tau = cbind(log(tau0), basis.prior2.tau[1], basis.prior2.tau[2]),
                             B.kappa = cbind(log(kappa0), basis.prior2.kappa[1], basis.prior2.kappa[2]),
                             theta.prior.mean = theta.prior2.mean,
                             theta.prior.prec = theta.prior2.prec)
  
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
  if(dayWeek == T){
    formula <- positiveResults ~ -1 + f(day, model = spde1) + f(dayWeek, model = 'iid', hyper = priorDayWeek, constr = T)
  }else{
    formula <- positiveResults ~ -1 + f(day, model = spde1) + f(dayId2, model = 'iid', hyper = priorDayWeek, constr = T)
  }
  if(linkType == "BB"){
    # Betabinomial
    m.spde <- inla(formula = formula, data = inla.stack.data(joint.stack),
                   control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                   family = "betabinomial",
                   Ntrials =  c(dataForModel[order(dayId), numberTest], rep(NA, length(xx))),
                   control.compute = list(config = TRUE),
                   control.family = list(hyper = BB.prior.rho))
  }else if(linkType == "NB"){
    # Negative binomial
    m.spde <- inla(formula = formula, data = inla.stack.data(joint.stack),
                   control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                   family = "nbinomial",
                   control.compute = list(config = TRUE),
                   control.family = list(hyper = NB.prior.rho))
  }
  
  # ---------------------------------------------------- #
  #                      GET SAMPLES                     #
  # ---------------------------------------------------- #
  cat("Estimating growth rate ... ")
  
  # Get samples of posterior distribution of parameters
  sampleList <- inla.posterior.sample(sizeSample, m.spde)
  
  # Extract relevant indexes from output
  rownms <- rownames(sampleList[[1]]$latent)
  set3 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "day"))
  dayIndexInSample <- set3[nonBoundaryIndices]
  
  # Create matrix of samples (and hyperparameters if needed)
  matrixSampleDays <- sapply(sampleList, function(x) x$latent[dayIndexInSample,1])
  # matrixSampleHyper <- sapply(sampleList, function(x) exp(x$hyperpar[c("Theta1 for day", "Theta2 for day")]))
  
  # Compute approximate derivative using windows (+-3)
  # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
  sampleDerivatives <- t(rbind(matrix(NA, nrow = 1, ncol = ncol(matrixSampleDays)),
                               ((matrixSampleDays[3,] - matrixSampleDays[1,])/2),
                               ((matrixSampleDays[5,] - matrixSampleDays[1,])/4),
                               (matrixSampleDays[7:numDays,] - matrixSampleDays[1:(numDays - 7 + 1),])/6,
                               ((matrixSampleDays[numDays,] - matrixSampleDays[numDays - 5 + 1,])/4),
                               ((matrixSampleDays[numDays,] - matrixSampleDays[numDays - 3 + 1,])/2),
                               matrix(NA, nrow = 1, ncol = ncol(matrixSampleDays))))
  
  # ---------------------------------------------------- #
  #                POSTERIOR GROWTH RATE                 #
  # ---------------------------------------------------- #
  # Compute log (or inv. logit) of posterior of Gaussian process derivative in transformed space
  cat("Computing posterior of growth rate ... ")
  if(linkType == "BB"){
    tempExpGP <- t(exp(matrixSampleDays))
    tempLogit <- tempExpGP/(1 + tempExpGP)
    tempList <- sampleDerivatives/(1 + tempLogit)
  }else if(linkType == "NB"){
    tempList <- sampleDerivatives
  }
  tempDoubling <- abs(log(2)/tempList)
  posteriorGrowth <- data.table(dayId = minDay:maxDay,
                                mean = colMeans(tempList),
                                sd = apply(tempList, 2, sd),
                                median = apply(tempList, 2, quantile, probs = 0.5, na.rm = T),
                                q0.025 = apply(tempList, 2, quantile, probs = 0.025, na.rm = T),
                                q0.975 = apply(tempList, 2, quantile, probs = 0.975, na.rm = T),
                                q0.25 = apply(tempList, 2, quantile, probs = 0.250, na.rm = T),
                                q0.75 = apply(tempList, 2, quantile, probs = 0.750, na.rm = T),
                                prob0 = apply(tempList, 2, function(x) sum(x >= 0)/sizeSample),
                                medianDoubT = apply(tempDoubling, 2, quantile, probs = 0.5, na.rm = T),
                                q0.025DoubT = apply(tempDoubling, 2, quantile, probs = 0.025, na.rm = T),
                                q0.975DoubT = apply(tempDoubling, 2, quantile, probs = 0.975, na.rm = T))
  setkey(posteriorGrowth, dayId)
  setkey(dateTable, dayId)
  posteriorGrowth[dateTable, ":="(date = i.date)]
  
  # ---------------------------------------------------- #
  #                 POSTERIOR INCIDENCE                  #
  # ---------------------------------------------------- #
  # Compute posterior of Gaussian process in real space (incidence or positivity)
  # (this section is slow)
  cat("Computing posterior of incidence... ")
  samplesGP <- m.spde$marginals.random$day[nonBoundaryIndices]
  if(linkType %in% c("NB")){
    transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(exp, x))
  }else if(linkType == "BB"){
    transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(function(x) exp(x)/(1 + exp(x)), x))
  }
  posteriorTransfGP <- data.table(dayId = minDay:maxDay,
                                  t(sapply(transformedSamples, function(x) inla.qmarginal(c(0.5, 0.025, 0.975, 0.25, 0.75), x))))
  setnames(posteriorTransfGP, c("V1", "V2", "V3", "V4", "V5"), c("median", "q0.025", "q0.975", "q0.25", "q0.75"))
  setkey(posteriorTransfGP, dayId)
  setkey(dateTable, dayId)
  posteriorTransfGP[dateTable, ":="(date = i.date)]
  setkey(posteriorTransfGP, date)
  setkey(countTable, date)
  posteriorTransfGP[countTable, ":="(positiveResults = i.positiveResults)]
  
  if(saveSamples == F){
    return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP))
  }else{
    return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP, matrixSampleDays = matrixSampleDays, sampleDerivatives = t(sampleDerivatives)))
  }
  
}
