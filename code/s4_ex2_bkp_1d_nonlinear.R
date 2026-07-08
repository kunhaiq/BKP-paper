## -------------------------------------------------------------------------
## Section 4, Example 2: One-dimensional nonlinear BKP example
##
## This script reproduces Example 2 in Section 4.
## The goal is to show that the BKP model can estimate a more challenging
## one-dimensional probability surface than the logistic curve used in
## Example 1.
##
## The true success probability function contains local oscillations and
## nonlinear variation. We generate binomial observations at Latin hypercube
## input locations, fit a BKP model, visualize the posterior mean and
## credible interval, and draw posterior sample curves from the fitted model.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------

## Load required packages.
## BKP provides the model fitting, plotting, and posterior simulation methods.
## The tgp package is used only for lhs(), which generates Latin hypercube
## design points over the input domain.
library(BKP)
library(tgp)

## Set console formatting options so that printed output is consistent with
## the R code style used in the manuscript.
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
## Data generation
## -------------------------------------------------------------------------

## Set the random seed to make the simulated design points, trial sizes,
## and binomial responses reproducible.
set.seed(123)

## Define the true success probability function for Example 2.
##
## This function is more nonlinear than the logistic curve in Example 1.
## It combines an exponential damping term with an oscillatory cosine term:
##
##   pi(x) = {1 + exp(-x^2) cos[10(1 - exp(-x))/(1 + exp(-x))]} / 2.
##
## The resulting probability surface remains inside [0, 1] but exhibits
## stronger local variation, making it a more challenging test function for
## probability-surface estimation.
true_pi_fun <- function(x) {
  (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
}

## Specify the number of training input locations.
n <- 30

## Define the one-dimensional input domain [-2, 2].
## Since the input dimension is one, Xbounds is a 1 by 2 matrix.
Xbounds <- matrix(c(-2, 2), nrow = 1)

## Generate n input locations using a Latin hypercube design over [-2, 2].
## The output X is an n by 1 matrix.
X <- lhs(n = n, rect = Xbounds)

## Evaluate the true success probability at the training input locations.
true_pi <- true_pi_fun(X)

## Generate heterogeneous binomial trial sizes.
## Each observation has a trial size randomly sampled from 1, ..., 100.
m <- sample(100, n, replace = TRUE)

## Generate binomial success counts from the true nonlinear probability
## function. At each input location x_i, y_i is drawn from
## Binomial(m_i, pi(x_i)).
y <- rbinom(n, size = m, prob = true_pi)


## -------------------------------------------------------------------------
## BKP model fitting
## -------------------------------------------------------------------------

## Fit the standard BKP model to the simulated binomial data.
## Because theta is not supplied, the kernel length-scale parameter is selected
## internally by leave-one-out cross-validation using the default Brier loss.
BKP_model_1D_2 <- fit_BKP(X, y, m, Xbounds = Xbounds)


## -------------------------------------------------------------------------
## Posterior mean and credible interval plot
## -------------------------------------------------------------------------

## Construct a dense prediction grid over the input domain.
## This grid is used to draw the true probability function on top of the
## posterior summary produced by plot().
Xnew <- matrix(seq(-2, 2, length = 100), ncol = 1)

## Evaluate the true probability function on the prediction grid.
true_pi <- true_pi_fun(Xnew)

## Save the posterior mean, pointwise 95% credible interval, observed
## proportions, and true probability curve to a PDF file.
pdf("code/figure/ex2.pdf", width = 6, height = 6)

## The plot() method for a fitted BKP object displays the posterior mean
## probability curve, pointwise credible interval, and observed proportions.
plot(BKP_model_1D_2)

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


## -------------------------------------------------------------------------
## Posterior simulation
## -------------------------------------------------------------------------

## Draw posterior samples of the latent success probability function from
## the fitted BKP model. The result contains three sampled probability curves
## evaluated on Xnew.
sim <- simulate(BKP_model_1D_2, nsim = 3, Xnew = Xnew)

## Save the posterior sample curves to a PDF file.
pdf("code/figure/ex2_sim.pdf", width = 6, height = 6)

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