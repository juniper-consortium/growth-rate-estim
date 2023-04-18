
# Code for the paper
# "Bayesian Estimation of real-time Epidemic Growth Rates using Gaussian Processes:
#  local dynamics of SARS-CoV-2 in England" (2022)

# Example of how to run model

library(data.table)
library(INLA)
library(ggplot2)

source("src/functions.R")
source("src/MCMC_functions.R")

# ---------------------------------------------------- #
#        Run model for validation section              #
# ---------------------------------------------------- #
parametersModel <- setParametersFn(linkType = "BB",
                                   prior2.sigma0 = 1,
                                   prior2.lengthscale0 = 50,
                                   theta.prior2.mean = c(1,1),
                                   theta.prior2.prec = diag(2),
                                   BB.prior.rho = list(overdispersion = list(prior = "gaussian", param = c(0, 0.4))),
                                   NB.prior.rho = list(size = list(prior = "gaussian", param = c(0, 0.01))),
                                   dayWeek = F,
                                   dayWeek.prior.prec = list(theta = list(prior = 'loggamma', param = c(1, 10))),
                                   sizeSample = 1000)

set.seed(98765)
# ---------- Scenario 1 ---------- #
popSize <- 1000000
numTests <- 0.1*popSize
numTests0 <- 0.1*numTests
numDays <- 200
load("SamplesForValidation.RData", verbose = T)
totalCases[, date := as.Date("2020-12-31") + dayId]
totalCases[, dayWeek := weekdays(date)]
totalCases[, numberTest := ifelse(dayWeek %in% c("Saturday", "Sunday"), 0.8, 1)*numTests]

countTable <- totalCases[simulation == 999, .(date, numberTest, casesInPopulation)]
countTable[, positiveResults := rhyper(nn = nrow(countTable), m = casesInPopulation, n = popSize - casesInPopulation, k = numberTest)]

parametersModel$params$linkType <- "NB"
outputSimulation1_positives <- runModelGrowthRate(countTable = countTable, parametersModel = parametersModel)
parametersModel$params$linkType <- "BB"
outputSimulation1_proportions <- runModelGrowthRate(countTable = countTable, parametersModel = parametersModel)

# ---------- Scenario 2 ---------- #
popSize <- 1000000
numTests <- 0.1*popSize
numTests0 <- 0.1*numTests
numDays <- 200
load("SamplesForValidation.RData", verbose = T)
totalCases[, date := as.Date("2020-12-31") + dayId]
totalCases[, dayWeek := weekdays(date)]
totalCases[, numberTest := ifelse(dayWeek %in% c("Saturday", "Sunday"), 0.8, 1)*round(numTests0 + dayId*(numTests - numTests0)/numDays)]

countTable <- totalCases[simulation == 999, .(date, numberTest, casesInPopulation)]
countTable[, positiveResults := rhyper(nn = nrow(countTable), m = casesInPopulation, n = popSize - casesInPopulation, k = numberTest)]

parametersModel$params$linkType <- "NB"
outputSimulation2_positives <- runModelGrowthRate(countTable = countTable, parametersModel = parametersModel)
parametersModel$params$linkType <- "BB"
outputSimulation2_proportions <- runModelGrowthRate(countTable = countTable, parametersModel = parametersModel)

# ---------- Plot Figure 1 * ---------- #
# (*Not exactly the same due to randomness not in the seed at the time of printing)
cutP <- numDays/2
rates <- c(0.03, -0.02)
dataToPlot <- rbind(data.table(outputSimulation1_positives$posteriorGrowth, scenarioLabel = "simulation 1", method = "positives model"),
                    data.table(outputSimulation1_proportions$posteriorGrowth, scenarioLabel = "simulation 1", method = "proportions model"),
                    data.table(outputSimulation2_positives$posteriorGrowth, scenarioLabel = "simulation 2", method = "positives model"),
                    data.table(outputSimulation2_proportions$posteriorGrowth, scenarioLabel = "simulation 2", method = "proportions model"))
dataToPlot[, real := (dayId <= cutP)*rates[1] + (dayId >= cutP)*rates[2]]
dataToPlot[, scenarioLevel := factor(scenarioLabel, levels = c("simulation 1", "simulation 2"))]
ggGRall <- ggplot(dataToPlot, aes(x = dayId)) + theme_laura() + facet_grid(scenarioLevel ~ method) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.5, fill = "gray80") +
  geom_ribbon(aes(ymin = q0.25, ymax = q0.75), alpha = 0.5, fill = "gray50") +
  geom_line(aes(y = median), size = 0.2) +
  geom_line(aes(y = real), linetype = 2, colour = "red") +
  labs(x = "time step", y = "growth rate posterior distribution") +
  theme(strip.background = element_rect(fill = NA, colour = "gray")) +
  coord_cartesian(ylim = c(-0.05, 0.06)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "gray", fill = NA),
        strip.background = element_rect(fill = NA, colour = "gray")) 
ggGRall

