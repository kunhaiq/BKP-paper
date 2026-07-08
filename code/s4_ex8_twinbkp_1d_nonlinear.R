## -------------------------------------------------------------------------
## Section 4, Example 8: One-dimensional TwinBKP nonlinear example
##
## This script reproduces Example 8 in Section 4 of the manuscript.
## It revisits the nonlinear binomial response setting from Example 2, but
## increases the sample size and replaces the full BKP model with the scalable
## TwinBKP approximation.
##
## The goal is to illustrate that TwinBKP uses the same response format and
## downstream workflow as BKP while replacing full kernel aggregation with a
## twinning-based global-local approximation.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------


## Load required packages.
## BKP provides the TwinBKP fitting, printing, and plotting methods.
## The tgp package is used for Latin hypercube sampling via lhs().
library(BKP)
library(tgp)

## Set console formatting options so that printed output follows the style
## used in the manuscript's R code examples.
options(
  prompt = "R> ",
  continue = "+  ",
  width = 70,
  useFancyQuotes = FALSE
)

## Create the output directory for figures if it does not already exist.
## This makes the script runnable from a clean clone of the repository.
dir.create("code/figure", recursive = TRUE, showWarnings = FALSE)


## -------------------------------------------------------------------------
## Data generation
## -------------------------------------------------------------------------

## Set the random seed to make the simulated design points, trial sizes,
## and binomial responses reproducible.
set.seed(123)

## Define the true success probability function.
##
## This is the nonlinear probability function used in Example 2:
##
##   pi(x) = 0.5 * {1 + exp(-x^2) cos[10(1 - exp(-x))/(1 + exp(-x))]}.
##
## The function is converted to a numeric vector so that it works cleanly
## with matrix inputs returned by lhs() and with prediction grids.
true_pi_fun <- function(x) {
  x <- as.numeric(x)
  0.5 * (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x))))
}

## Specify the number of training input locations.
## Compared with Example 2, the sample size is increased to demonstrate the
## scalable TwinBKP interface.
n <- 500

## Define the one-dimensional input domain [-2, 2].
Xbounds <- matrix(c(-2, 2), nrow = 1)

## Generate n input locations using a Latin hypercube design over [-2, 2].
X <- lhs(n = n, rect = Xbounds)

## Evaluate the true nonlinear probability function at the training locations.
true_pi <- true_pi_fun(X)

## Generate heterogeneous binomial trial sizes.
## Each observation has a trial size sampled from 1, ..., 100.
m <- sample(100, n, replace = TRUE)

## Generate binomial success counts from the true probability function.
## At each input location x_i, y_i is drawn from Binomial(m_i, pi(x_i)).
y <- rbinom(n, size = m, prob = true_pi)


## -------------------------------------------------------------------------
## TwinBKP model fitting
## -------------------------------------------------------------------------

## Fit the TwinBKP model to the simulated binomial data.
## The function uses a twinning-selected global subset and prediction-specific
## local neighbours. Since theta_g and theta_l are not supplied, the global
## kernel parameter is tuned on the global subset and the local range parameter
## is set internally according to the default TwinBKP rule.
TwinBKP_model_1D_2 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)

## Print a concise summary of the fitted TwinBKP model, including the global
## and local kernel settings, selected kernel parameters, loss value, global
## subset size, local-neighbour size, and twinning settings.
print(TwinBKP_model_1D_2)


## -------------------------------------------------------------------------
## Posterior summary plot
## -------------------------------------------------------------------------

## Construct a dense prediction grid over the input domain.
## This grid is used to add the true probability curve to the fitted TwinBKP
## posterior summary plot.
Xnew <- matrix(seq(-2, 2, length = 100), ncol = 1)

## Evaluate the true probability function on the prediction grid.
true_pi <- true_pi_fun(Xnew)

## Save the fitted TwinBKP posterior summaries to a PDF file.
## The plot() method displays the posterior mean, pointwise credible interval,
## observed proportions, and the twinning-selected global subset.
pdf("code/figure/ex8.pdf", width = 9, height = 6)
plot(TwinBKP_model_1D_2)

## Add the true probability curve for comparison.
lines(Xnew, true_pi, col = "black", lwd = 2)

## Add a legend entry for the true probability curve.
legend(
  x = -2.22, y = 0.85,
  legend = "True Probability",
  col = "black",
  lwd = 2,
  bty = "n",
  inset = 0.02
)

## Close the PDF graphics device.
dev.off()