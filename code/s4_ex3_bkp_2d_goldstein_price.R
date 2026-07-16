## -------------------------------------------------------------------------
## Section 4, Example 3: Two-dimensional BKP example and timing comparison
##
## This script reproduces Example 3 in Section 4 of the manuscript.
## It contains two main components:
##
##   1. A two-dimensional BKP fit based on a probability surface obtained
##      from a rescaled Goldstein--Price-type test function.
##
##   2. A computational timing comparison between BKP and logistic Gaussian
##      process models under fixed and optimized kernel hyperparameters.
##
## By default, the timing experiment is not re-run. Instead, the script reads
## the pre-computed average timing results from code/result/elapsed_time_avg.csv.
## Set run_elapsed_time = TRUE to recompute the full timing experiment.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------


## Load required packages.
## BKP provides the Beta Kernel Process fitting, prediction, and plotting
## functions. The tgp package is used for Latin hypercube sampling via lhs().
## The gplite package is used to fit logistic Gaussian process models for
## computational comparison.
library(BKP)
library(tgp)
library(gplite)

## Set console formatting options so that printed output follows the style
## used in the manuscript's R code examples.
options(
  prompt = "R> ",
  continue = "+  ",
  width = 70,
  useFancyQuotes = FALSE
)

## Create output directories if they do not already exist. This makes the
## script runnable from a clean clone of the repository.
dir.create("code/figure", recursive = TRUE, showWarnings = FALSE)
dir.create("code/result", recursive = TRUE, showWarnings = FALSE)


## -------------------------------------------------------------------------
## Part I: Two-dimensional probability-surface estimation
## -------------------------------------------------------------------------

## Set the random seed to make the simulated design points, trial sizes,
## and binomial responses reproducible.
set.seed(123)

## Define the two-dimensional true probability function.
##
## The latent surface is constructed from a rescaled version of the
## Goldstein--Price function. The input X is assumed to lie in [0, 1]^2.
## The variables x1 and x2 are first mapped from [0, 1] to [-2, 2], and
## then used to evaluate the Goldstein--Price-type latent function.
##
## The constants m and s center and scale the log-transformed surface.
## Finally, pnorm() maps the latent real-valued function to a probability
## in [0, 1].
true_pi_fun <- function(X) {
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  m <- 8.6928
  s <- 2.4269
  x1 <- 4*X[,1]- 2
  x2 <- 4*X[,2]- 2
  a <- 1 + (x1 + x2 + 1)^2 *
    (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
  b <- 30 + (2*x1- 3*x2)^2 *
    (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
  f <- log(a*b)
  f <- (f- m)/s
  return(pnorm(f))  # Transform to probability
}

## Specify the number of training input locations.
n <- 100

## Define the two-dimensional input domain [0, 1]^2.
## Each row of Xbounds gives the lower and upper bounds of one input dimension.
Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)

## Generate a Latin hypercube design with n input locations over [0, 1]^2.
X <- lhs(n = n, rect = Xbounds)

## Evaluate the true probability surface at the training input locations.
true_pi <- true_pi_fun(X)

## Generate heterogeneous binomial trial sizes.
## Each observation has a trial size sampled from 1, ..., 100.
m <- sample(100, n, replace = TRUE)

## Generate binomial success counts at each input location.
## Conditional on X_i and m_i, y_i is drawn from Binomial(m_i, pi(X_i)).
y <- rbinom(n, size = m, prob = true_pi)

## Fit the standard BKP model to the simulated two-dimensional binomial data.
## Since theta is not supplied, the kernel length-scale parameter is selected
## internally by leave-one-out cross-validation using the default Brier loss.
BKP_model_2D <- fit_BKP(X, y, m, Xbounds = Xbounds)

## Print a concise summary of the fitted BKP model, including the sample size,
## input dimension, kernel type, optimized kernel parameter, loss value, and
## prior type.
print(BKP_model_2D)

## Save the BKP posterior summary plots to a PDF file.
## For a two-dimensional BKP model, plot() produces contour-style summaries
## such as posterior mean, posterior variance, and credible interval bounds.
pdf("code/figure/ex3.pdf", width = 9, height = 8)
plot(BKP_model_2D)
dev.off()

## Plot the true probability surface on a dense grid for comparison with the
## fitted BKP posterior summaries.
pdf("code/figure/ex3_true.pdf", width = 4.5, height = 4)

## Construct a regular 200 by 200 prediction grid over [0, 1]^2.
Xnew1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 200)
Xnew2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 200)
Xnew <- expand.grid(Xnew1 = Xnew1, Xnew2 = Xnew2)

## Evaluate the true probability function on the grid.
true_pi <- true_pi_fun(Xnew)

## Store the grid coordinates and true probabilities in a data frame that
## matches the input format required by the plotting helper.
df <- data.frame(x1 = Xnew$Xnew1, x2 = Xnew$Xnew2,
                 True = true_pi)

## Use BKP's internal two-dimensional plotting helper to reproduce the
## contour-style plot of the true probability surface. The function is
## non-exported but documented in the package; it is used here only to match
## the manuscript figure.
print(BKP:::my_2D_plot_fun("True", title = "True Probability", data = df))
dev.off()

## Generate ten new test input locations and print BKP predictions at these
## locations. The predict() method returns posterior summaries such as the
## posterior mean, variance, and credible interval bounds.
Xnew <- lhs(n = 10, rect = Xbounds)
print(predict(BKP_model_2D, Xnew))


## -------------------------------------------------------------------------
## Part II: Computational timing comparison between BKP and LGP
## -------------------------------------------------------------------------

## This block controls whether to re-run the full timing experiment.
##
## When run_elapsed_time is FALSE, the script reads pre-computed average
## elapsed times from code/result/elapsed_time_avg.csv. This is the default
## mode for fast reproduction of the timing figure.
##
## When run_elapsed_time is TRUE, the script re-runs the full benchmark:
## five methods, five sample sizes, and twenty repetitions. This can take
## many hours, especially for the optimized LGP model at n = 5000.
run_elapsed_time = FALSE

## Sample sizes used in the timing comparison.
n_vals <- c(200, 500, 1000, 2000, 5000)

## File path for the pre-computed average timing results.
time_file <- "code/result/elapsed_time_avg.csv"

## If run_elapsed_time is TRUE, run the full benchmark and save a new timing
## result file. Otherwise, read the pre-computed timing results.
if(run_elapsed_time){
  
  ## Initialize matrices to store elapsed times.
  ## Rows correspond to independent repetitions, and columns correspond to
  ## the sample sizes in n_vals.
  elapsed_time_given <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_single_start <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_multi_start <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_given_gp <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_optim_gp <- matrix(NA, nrow = 20, ncol = 5)
  
  ## Loop over 20 repetitions and five sample sizes.
  ## Each iteration simulates a new dataset and records the elapsed time for
  ## BKP and LGP fitting under different hyperparameter settings.
  for (i in 1:20) {
    for (j in 1:5) {
      
      ## Use a deterministic seed for each repetition/sample-size combination
      ## to make the benchmark reproducible.
      set.seed((i - 1) * 5 + j)
      cat("i =", i, "j =", j, ":")
      
      ## Select the current sample size.
      n <- n_vals[j]
      
      ## Generate a new two-dimensional training dataset.
      Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
      X <- lhs(n = n, rect = Xbounds)
      true_pi <- true_pi_fun(X)
      m <- sample(100, n, replace = TRUE)
      y <- rbinom(n, size = m, prob = true_pi)
      
      ## Measure elapsed time for BKP with a fixed kernel length-scale theta.
      ## This isolates the cost of full-kernel BKP fitting when no
      ## hyperparameter optimization is performed.
      elapsed_time_given[i, j] <- system.time({
        fit <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 1)
      })["elapsed"]
      cat(" BKP_given =", elapsed_time_given[i, j])
      
      ## Measure elapsed time for BKP when theta is optimized using a single
      ## starting value in the LOOCV-based optimization.
      elapsed_time_single_start[i, j] <- system.time({
        fit <- fit_BKP(X, y, m, Xbounds = Xbounds, n_multi_start = 1)
      })["elapsed"]
      cat(", BKP_single_start =", elapsed_time_single_start[i, j])
      
      ## Measure elapsed time for BKP under the default multi-start
      ## optimization strategy.
      elapsed_time_multi_start[i, j] <- system.time({
        fit <- fit_BKP(X, y, m, Xbounds = Xbounds)
      })["elapsed"]
      cat(", BKP_multi_start =", elapsed_time_multi_start[i, j])
      
      ## Measure elapsed time for the logistic Gaussian process model with
      ## given hyperparameters, using gplite.
      elapsed_time_given_gp[i, j] <- system.time({
        gp <- gp_init(cf = cf_sexp(), lik = lik_binomial())
        gp <- gp_fit(gp, X, y, trials = m)
      })["elapsed"]
      cat(", GP_given =", elapsed_time_given_gp[i, j])
      
      ## Measure elapsed time for the logistic Gaussian process model with
      ## hyperparameter optimization.
      elapsed_time_optim_gp[i, j] <- system.time({
        gp <- gp_init(cf = cf_sexp(), lik = lik_binomial())
        gp <- gp_optim(gp, X, y, trials = m, verbose = FALSE)
      })["elapsed"]
      cat(", GP_optim =", elapsed_time_optim_gp[i, j], "\n")
    }
  }
  
  ## Compute the average elapsed time across the 20 repetitions for each
  ## method and each sample size.
  avg_time <- matrix(NA, nrow = 5, ncol = 5)
  colnames(avg_time) <- paste0("n=", n_vals)
  rownames(avg_time) <- c("given", "single_start", "multi_start", "given_gp", "optim_gp")
  avg_time[1, ] <- colMeans(elapsed_time_given)
  avg_time[2, ] <- colMeans(elapsed_time_single_start)
  avg_time[3, ] <- colMeans(elapsed_time_multi_start)
  avg_time[4, ] <- colMeans(elapsed_time_given_gp)
  avg_time[5, ] <- colMeans(elapsed_time_optim_gp)
  
  ## Save the newly computed timing results to a separate file.
  ## The suffix "_new.csv" avoids overwriting the pre-computed benchmark
  ## file used for default reproduction.
  time_file_new <- sub("\\.csv$", "_new.csv", time_file)
  write.csv(avg_time, file = time_file_new)
  cat("Timing results saved to:", time_file_new, "\n")
}else if(file.exists(time_file)){
  
  ## In the default mode, read the pre-computed average elapsed times.
  ## This allows the timing figure to be reproduced quickly without running
  ## the expensive benchmark.
  avg_time <- as.matrix(read.csv(time_file, row.names = 1))
  rownames(avg_time) <- c("given", "single_start", "multi_start", "given_gp", "optim_gp")
  colnames(avg_time) <- c("n=200", "n=500", "n=1000", "n=2000", "n=5000")
} else {
  
  ## Stop with an informative error if the pre-computed timing file is absent.
  ## In that case, the user must either provide the CSV file or set
  ## run_elapsed_time = TRUE to regenerate the timing results.
  stop("Pre-computed timing file not found. Please run with run_elapsed_time = TRUE first.")
}


## -------------------------------------------------------------------------
## Plot the timing comparison
## -------------------------------------------------------------------------

## Save the timing comparison figure to a PDF file.
## The figure has two panels:
##   left:  fixed-hyperparameter BKP versus fixed-hyperparameter LGP;
##   right: optimization-based BKP and LGP methods.
pdf("code/figure/elapsed_time.pdf", width = 10, height = 5)

## Use a two-panel plotting layout.
par(mfrow = c(1, 2), mar = c(4,4,3,1))

## Define colors and plotting symbols for the five timing curves.
cols <- c("blue", "skyblue", "darkblue", "red", "orange")
pchs <- c(16, 17, 16, 15, 15)

## Left panel: compare BKP and LGP when hyperparameters are fixed.
## Both axes are shown on a log scale to make computational scaling easier
## to assess across sample sizes.
plot(n_vals, avg_time["given", ], log = "xy", type = "b", col = cols[1], pch = pchs[1], lwd = 2,
     xlab = "n", ylab = "Time (seconds, log scale)",
     main = "Fixed Hyperparameter",
     ylim = range(avg_time[c("given", "given_gp"), ]))

## Add the fixed-hyperparameter LGP timing curve.
lines(n_vals, avg_time["given_gp", ], type = "b", col = cols[4], pch = pchs[4], lwd = 2)

## Add reference complexity curves.
## The O(n^2) line is anchored at the first BKP timing value, while the
## O(n^3) line is anchored at the first LGP timing value.
ref_bkp <- avg_time["given", 1] * (n_vals / n_vals[1])^2
ref_gp  <- avg_time["given_gp", 1] * (n_vals / n_vals[1])^3
lines(n_vals, ref_bkp, col = "blue", lty = 2)
lines(n_vals, ref_gp, col = "red", lty = 4)

## Add the legend for the fixed-hyperparameter comparison.
legend("topleft", bty = "n",
       legend = c("BKP", "LGP", "O(n^2)", "O(n^3)"),
       col = c(cols[c(1,4)], "blue", "red"),
       pch = c(pchs[c(1,4)], NA, NA),
       lty = c(1,1,2,4),
       lwd = c(2,2,1,1),
       bg = "white")

## Right panel: compare optimization-based methods.
## The BKP single-start and multi-start settings are compared with optimized
## LGP fitting.
plot(n_vals, avg_time["single_start", ], log = "xy", type = "b", col = cols[2], pch = pchs[2], lwd = 2,
     xlab = "n", ylab = "Time (seconds, log scale)",
     main = "Optimization-based Methods",
     ylim = range(avg_time[c("single_start", "multi_start", "optim_gp"), ]))

## Add the BKP multi-start timing curve.
lines(n_vals, avg_time["multi_start", ], type = "b", col = cols[3], pch = pchs[3], lwd = 2)

## Add the optimized LGP timing curve.
lines(n_vals, avg_time["optim_gp", ], type = "b", col = cols[5], pch = pchs[5], lwd = 2)

## Add reference complexity curves for the optimization-based comparison.
## The O(n^2) reference is anchored to BKP single-start timing, and the
## O(n^3) reference is anchored to optimized LGP timing.
ref_bkp2 <- avg_time["single_start", 1] * (n_vals / n_vals[1])^2
ref_gp2  <- avg_time["optim_gp", 1] * (n_vals / n_vals[1])^3
lines(n_vals, ref_bkp2, col = "skyblue", lty = 2)
lines(n_vals, ref_gp2, col = "orange", lty = 4)

## Add the legend for the optimization-based comparison.
legend("topleft", bty = "n",
       legend = c("BKP (single start)", "BKP (multi-start)", "LGP", "O(n^2)", "O(n^3)"),
       col = c(cols[c(2,3,5)], "skyblue", "orange"),
       pch = c(pchs[c(2,3,5)], NA, NA),
       lty = c(1,1,1,2,4),
       lwd = c(2,2,2,1,1),
       bg = "white")

## Restore the default single-panel plotting layout.
par(mfrow = c(1,1))

## Close the PDF graphics device.
dev.off()