## -------------------------------------------------------------------------
## Section 4, Example 1: One-dimensional BKP example
##
## This script reproduces the one-dimensional BKP example in Section 4.
## It contains two parts:
##   1. Binomial probability-surface estimation under a logistic probability
##      function.
##   2. A binary classification illustration showing the effect of the prior
##      precision parameter r0.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------

## Load required packages.
## BKP provides the Beta Kernel Process fitting, prediction, plotting,
## and simulation functions. The tgp package is used here only for the lhs()
## function, which generates a Latin hypercube design over the input domain.
library(BKP)
library(tgp)

## Set display options so that console output is formatted similarly to
## the examples shown in the manuscript.
options(
  prompt = "R> ",
  continue = "+  ",
  width = 70,
  useFancyQuotes = FALSE
)

## Create output directories if they do not already exist. This makes the
## script runnable from a clean clone of the repository.
dir.create("code/figure", recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------------
## Part I: Binomial probability-surface estimation
## -------------------------------------------------------------------------

## Set the random seed to make the simulated design points, trial sizes,
## and binomial responses reproducible.
set.seed(123)

## Define the true success probability function.
## This is a smooth logistic curve on the input interval [-2, 2]:
##
##     pi(x) = 1 / {1 + exp(-3x)}.
##
## The BKP model is fitted using simulated binomial observations generated
## from this probability function.
true_pi_fun <- function(x) {
  1 / (1 + exp(-3 * x))
}

## Specify the number of training input locations.
n <- 7

## Define the input domain. Here the problem is one-dimensional, so Xbounds
## is a 1 by 2 matrix containing the lower and upper bounds of x.
Xbounds <- matrix(c(-2, 2), nrow = 1)

## Generate n input locations using a Latin hypercube design over [-2, 2].
## The output X is an n by 1 matrix.
X <- lhs(n = n, rect = Xbounds)

## Evaluate the true success probability at the training input locations.
true_pi <- true_pi_fun(X)

## Generate heterogeneous binomial trial sizes.
## Each observation has a trial size randomly sampled from 1, ..., 100.
m <- sample(100, n, replace = TRUE)

## Generate binomial success counts from the true probability function.
## For each location x_i, y_i is drawn from Binomial(m_i, pi(x_i)).
y <- rbinom(n, size = m, prob = true_pi)

## Fit the standard BKP model.
## Since theta is not supplied, the kernel length-scale parameter is selected
## internally by leave-one-out cross-validation using the default Brier loss.
BKP_model_1D_1 <- fit_BKP(X, y, m, Xbounds = Xbounds)

## Print a concise summary of the fitted BKP model, including the number of
## observations, input dimension, kernel type, optimized kernel parameter,
## achieved loss value, and prior type.
print(BKP_model_1D_1)

## Construct a dense prediction grid over the input domain.
## This grid is used only for drawing the true probability curve on top of
## the fitted BKP plot.
Xnew <- matrix(seq(-2, 2, length = 100), ncol = 1)

## Evaluate the true probability function on the prediction grid.
true_pi <- true_pi_fun(Xnew)

## Save the posterior mean, 95% credible interval, observed proportions,
## and true probability curve to a PDF file.
pdf("code/figure/ex1.pdf", width = 6, height = 6)

## The plot() method for a fitted BKP object displays the posterior mean
## probability curve, pointwise credible interval, and observed proportions.
plot(BKP_model_1D_1)

## Add the true underlying probability function for comparison.
lines(Xnew, true_pi, col = "black", lwd = 2)

## Add a legend entry for the true probability curve.
legend(
  x = -2.24, y = 0.89,
  legend = "True Probability",
  col = "black",
  lwd = 2,
  bty = "n",
  inset = 0.02
)

## Close the PDF graphics device.
dev.off()

## Draw posterior samples of the latent success probability function from
## the fitted BKP model. The result contains three sampled probability curves
## evaluated on Xnew.
sim <- simulate(BKP_model_1D_1, nsim = 3, Xnew = Xnew)

## Save the posterior sample curves to a PDF file.
pdf("code/figure/ex1_sim.pdf", width = 6, height = 6)

## Plot the simulated probability curves.
## Each column of sim$samples corresponds to one posterior draw.
matplot(
  Xnew, sim$samples,
  type = "l",
  lty = 1,
  lwd = 1.5,
  col = rainbow(ncol(sim$samples)),
  xlab = "x",
  ylab = "Probability",
  main = "Simulated Probability Curves"
)

## Add a legend identifying the posterior sample curves.
legend(
  "topleft",
  legend = paste("sample", 1:ncol(sim$samples)),
  col = rainbow(ncol(sim$samples)),
  lty = 1,
  lwd = 2,
  bty = "n"
)

## Close the PDF graphics device.
dev.off()


## -------------------------------------------------------------------------
## Part II: Binary classification and the effect of prior precision
## -------------------------------------------------------------------------

## Reset the random seed so that the classification example is reproducible
## independently of the first binomial example.
set.seed(123)

## Use a larger number of input locations for the classification illustration.
n <- 20

## The classification example uses the same one-dimensional input domain.
Xbounds <- matrix(c(-2, 2), nrow = 1)

## Generate training input locations using Latin hypercube sampling.
X <- lhs(n = n, rect = Xbounds)

## Evaluate the true logistic probability function at the training locations.
true_pi <- true_pi_fun(X)

## For binary classification, each observation has one Bernoulli trial.
m <- rep(1, n)

## Generate deterministic binary labels from the true probability function.
## The label is 1 when pi(x) > 0.5 and 0 otherwise. This produces a clean
## classification boundary at x = 0.
y <- as.numeric(true_pi > 0.5)

## Fit a BKP classifier using a fixed prior with a small prior precision.
## A small r0 gives the prior weak influence, allowing the fitted probability
## surface to be driven mainly by the observed binary labels.
BKP_model_1D_1_class_1 <- fit_BKP(
  X, y, m,
  Xbounds = Xbounds,
  prior = "fixed",
  r0 = 0.01,
  loss = "log_loss"
)

## Construct a dense grid for adding the true probability curve to the plot.
Xnew <- matrix(seq(-2, 2, length = 100), ncol = 1)

## Evaluate the true probability function on the grid.
true_pi <- true_pi_fun(Xnew)

## Save the fitted classification probability curve for r0 = 0.01.
pdf("code/figure/ex1class001.pdf", width = 6, height = 6)

## The plot() method displays the fitted posterior mean probability curve,
## credible interval, observed labels, and classification threshold.
plot(BKP_model_1D_1_class_1)

## Add the true logistic probability curve for reference.
lines(Xnew, true_pi, col = "black", lwd = 2)

## Close the PDF graphics device.
dev.off()

## Fit a second BKP classifier using a fixed prior with r0 = 2.
## Compared with r0 = 0.01, this prior is more influential and can lead to
## a smoother and less decisive fitted classification surface.
BKP_model_1D_1_class_2 <- fit_BKP(
  X, y, m,
  Xbounds = Xbounds,
  prior = "fixed",
  r0 = 2,
  loss = "log_loss"
)

## Save the fitted classification probability curve for r0 = 2.
pdf("code/figure/ex1class2.pdf", width = 6, height = 6)

## Plot the fitted BKP classifier under the stronger prior.
plot(BKP_model_1D_1_class_2)

## Add the true logistic probability curve for reference.
lines(Xnew, true_pi, col = "black", lwd = 2)

## Close the PDF graphics device.
dev.off()