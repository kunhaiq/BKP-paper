## -------------------------------------------------------------------------
## Section 5, Application 1: Loa loa infection prevalence mapping
##
## This script reproduces the first real-data application in Section 5 of
## the manuscript. It analyzes village-level Loa loa infection data from
## North Cameroon using aggregated binomial responses.
##
## At each surveyed location, the response records the number of infected
## individuals and the total number examined. Longitude and latitude are used
## as the two spatial predictors. The observations are randomly divided into
## training and testing sets, after which two models are compared:
##
##   1. A Beta Kernel Process (BKP) model fitted directly to the aggregated
##      binomial counts.
##
##   2. A logistic Gaussian process (LGP) model fitted with the gplite package.
##
## Predictive accuracy is evaluated on the withheld test set using the Brier
## score, with the observed infection proportion y / m treated as the empirical
## target probability. The fitted models are subsequently evaluated over a
## land-only spatial grid to visualize posterior predictive means and
## variances.
##
## All paths are relative to the repository root. Generated figures are saved
## to code/figure/.
## -------------------------------------------------------------------------


## -------------------------------------------------------------------------
## Required packages
## -------------------------------------------------------------------------

## BKP provides the Beta Kernel Process fitting and prediction methods.
## gplite provides the logistic Gaussian process comparison model.
## RiskMap provides the Loa loa dataset.
## sf and rnaturalearth are used to define the geographic study region and
## remove prediction-grid points that fall outside land areas.
## ggplot2, gridExtra, and viridis support visualization and figure assembly.
library(BKP)
library(gplite)
library(ggplot2)
library(gridExtra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RiskMap)
library(viridis)


## Set console formatting options so that printed output follows the style
## used in the manuscript's R code examples.
options(
  prompt = "R> ",
  continue = "+  ",
  width = 70,
  useFancyQuotes = FALSE
)

## Create the output directory if it does not already exist. This allows the
## script to be run independently from a clean clone of the repository.
dir.create(
  "code/figure",
  recursive = TRUE,
  showWarnings = FALSE
)

## Set the random seed to make the train-test split reproducible.
set.seed(123)


## -------------------------------------------------------------------------
## Data preparation
## -------------------------------------------------------------------------

## Load the village-level Loa loa survey data from the RiskMap package.
## Each row represents one surveyed location and contains:
##
##   - longitude and latitude;
##   - the number of infected individuals;
##   - the total number of individuals examined;
##   - additional environmental variables not used in this application.
data("loaloa", package = "RiskMap")

## Extract longitude and latitude as the spatial input matrix.
## The response y contains the number of infected individuals, while m contains
## the corresponding binomial trial sizes.
X <- as.matrix(loaloa[, c("LONGITUDE", "LATITUDE")])
rownames(X) <- NULL

y <- loaloa$NO_INF
m <- loaloa$NO_EXAM

## Randomly allocate 70% of the survey locations to the training set and retain
## the remaining 30% for out-of-sample predictive evaluation.
train_idx <- sample(
  seq_len(nrow(loaloa)),
  size = 0.7 * nrow(loaloa)
)

X_train <- X[train_idx, , drop = FALSE]
y_train <- y[train_idx]
m_train <- m[train_idx]

X_test <- X[-train_idx, , drop = FALSE]
y_test <- y[-train_idx]
m_test <- m[-train_idx]

## Compute the observed infection proportions for visualization.
p <- y / m

## Define the rectangular geographic domain used for normalization, prediction,
## and plotting. The first row contains longitude bounds and the second row
## contains latitude bounds.
Xbounds <- matrix(
  c(7.8, 15.3, 3.1, 7.0),
  ncol = 2,
  byrow = TRUE
)

## Compute empirical infection proportions in the test set. These values are
## used as the target probabilities in the Brier-score comparison.
pi_test <- y_test / m_test


## -------------------------------------------------------------------------
## Geographic base map
## -------------------------------------------------------------------------

## Obtain an sf representation of the African land polygons. This object is
## used both as a background map and as a land mask for spatial prediction.
africa <- ne_countries(
  continent = "africa",
  scale = "medium",
  returnclass = "sf"
)


## -------------------------------------------------------------------------
## Observed prevalence map
## -------------------------------------------------------------------------

## Construct a plotting data frame containing the observed infection
## proportion, trial size, geographic coordinates, and train-test membership.
df_map <- data.frame(
  lon = X[, 1],
  lat = X[, 2],
  p = p,
  n = m,
  set = ifelse(
    seq_len(nrow(loaloa)) %in% train_idx,
    "Train",
    "Test"
  )
)

## Plot the observed infection proportions.
##
## Point color represents the empirical infection proportion y / m.
## Point size represents the number of individuals examined.
## Point shape distinguishes training and testing observations.
pmap <- ggplot() +
  geom_sf(
    data = africa,
    fill = "gray95",
    color = "gray60"
  ) +
  geom_point(
    data = df_map,
    aes(
      x = lon,
      y = lat,
      color = p,
      size = n,
      shape = set
    ),
    alpha = 0.8
  ) +
  scale_color_viridis_c(
    name = "Proportion y/m"
  ) +
  scale_size_continuous(
    name = "Trial count m"
  ) +
  scale_shape_manual(
    name = "Dataset",
    values = c(
      "Train" = 16,
      "Test" = 17
    )
  ) +
  coord_sf(
    xlim = Xbounds[1, ],
    ylim = Xbounds[2, ],
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Loaloa infection proportion (y/m): Train vs Test",
    x = "Longitude",
    y = "Latitude"
  )

## Save the observed prevalence map used in the manuscript.
ggsave(
  "code/figure/Loaloa_map.pdf",
  plot = pmap,
  width = 10,
  height = 5,
  dpi = 300
)


## -------------------------------------------------------------------------
## Model fitting and test-set prediction
## -------------------------------------------------------------------------

## Initialize a list for storing test-set predictive probabilities from the
## competing models.
all_preds <- list()


## Fit the BKP model directly to the aggregated binomial counts.
##
## The inputs are the spatial coordinates, the observed infected counts, and
## the corresponding numbers examined. Default kernel, prior, and
## hyperparameter-tuning settings are used.
bkp_fit <- fit_BKP(
  X_train,
  y_train,
  m_train,
  Xbounds = Xbounds
)

## Obtain BKP posterior mean predictions at the withheld test locations.
all_preds$BKP <- predict(
  bkp_fit,
  Xnew = X_test
)$mean

## Print the fitted BKP model summary reported in the manuscript workflow.
print(summary(bkp_fit))


## Fit the logistic Gaussian process comparison model.
##
## The binomial likelihood uses y_train as the number of successes and m_train
## as the location-specific trial counts.
gp <- gp_init(
  cf = cf_sexp(),
  lik = lik_binomial()
)

gp <- gp_optim(
  gp,
  X_train,
  y_train,
  trials = m_train,
  verbose = FALSE
)

## Obtain transformed LGP predictive means on the probability scale.
all_preds$LGP <- gp_pred(
  gp,
  X_test,
  transform = TRUE
)$mean


## -------------------------------------------------------------------------
## Out-of-sample predictive comparison
## -------------------------------------------------------------------------

## Compare the predicted probabilities with the empirical infection
## proportions in the withheld test set.
##
## For each method, the Brier score is computed as the mean squared difference
## between the predicted probability and y_test / m_test.
cat("\n=============================================================\n")
cat("  Comparison: Loaloa (BKP vs LGP)\n")
cat("=============================================================\n")

for (nm in names(all_preds)) {
  brier_score <- mean(
    (all_preds[[nm]] - pi_test)^2
  )
  
  cat(
    sprintf(
      "  %-8s | Brier: %.6f\n",
      nm,
      brier_score
    )
  )
}


## -------------------------------------------------------------------------
## Land-only spatial prediction grid
## -------------------------------------------------------------------------

## Construct a dense rectangular longitude-latitude grid over the study region.
x1_seq <- seq(
  Xbounds[1, 1],
  Xbounds[1, 2],
  length.out = 300
)

x2_seq <- seq(
  Xbounds[2, 1],
  Xbounds[2, 2],
  length.out = 260
)

grid <- expand.grid(
  lon = x1_seq,
  lat = x2_seq
)

## Convert the prediction grid to an sf object using the WGS84 coordinate
## reference system.
grid_sf <- st_as_sf(
  grid,
  coords = c("lon", "lat"),
  crs = 4326
)

## Retain only grid points that intersect the African land polygons. This
## prevents predictions over ocean locations from appearing in the maps.
on_land <- lengths(
  st_intersects(grid_sf, africa)
) > 0

grid_land <- grid[
  on_land,
  ,
  drop = FALSE
]

X_grid <- as.matrix(grid_land)


## -------------------------------------------------------------------------
## BKP prediction over the spatial grid
## -------------------------------------------------------------------------

## Obtain BKP posterior means and variances over the land-only prediction grid.
pred_bkp <- predict(
  bkp_fit,
  Xnew = X_grid
)

grid_land$BKP_mean <- as.vector(
  pred_bkp$mean
)

grid_land$BKP_var <- as.vector(
  pred_bkp$variance
)


## -------------------------------------------------------------------------
## LGP prediction over the spatial grid
## -------------------------------------------------------------------------

## Evaluate LGP predictions in batches to reduce peak memory use when
## processing the dense spatial grid.
chunk_size <- 5000

gm <- rep(
  NA_real_,
  nrow(X_grid)
)

gv <- rep(
  NA_real_,
  nrow(X_grid)
)

for (k in seq(
  from = 1,
  to = nrow(X_grid),
  by = chunk_size
)) {
  idx_end <- min(
    k + chunk_size - 1,
    nrow(X_grid)
  )
  
  idx <- k:idx_end
  
  chunk <- X_grid[
    idx,
    ,
    drop = FALSE
  ]
  
  pc <- gp_pred(
    gp,
    chunk,
    transform = TRUE,
    var = TRUE
  )
  
  gm[idx] <- as.vector(
    pc$mean
  )
  
  gv[idx] <- as.vector(
    pc$var
  )
}

grid_land$LGP_mean <- gm
grid_land$LGP_var <- gv


## Store the spatial prediction results in a common structure so that the two
## models can be visualized using the same plotting workflow.
raster_preds <- list(
  BKP = list(
    mean = grid_land$BKP_mean,
    variance = grid_land$BKP_var
  ),
  LGP = list(
    mean = grid_land$LGP_mean,
    variance = grid_land$LGP_var
  )
)


## -------------------------------------------------------------------------
## Common color-scale limits
## -------------------------------------------------------------------------

## Determine a common upper limit for the predictive mean maps. Using the same
## color scale across methods permits direct visual comparison of the BKP and
## LGP prevalence surfaces.
all_pmax <- sapply(
  raster_preds,
  function(x) {
    max(
      x$mean,
      na.rm = TRUE
    )
  }
)

plim <- ceiling(
  max(all_pmax) * 10
) / 10

## Determine a common upper limit for the predictive variance maps.
all_vmax <- sapply(
  raster_preds,
  function(x) {
    max(
      x$variance,
      na.rm = TRUE
    )
  }
)

vlim <- ceiling(
  max(all_vmax) * 10
) / 10


## -------------------------------------------------------------------------
## Model-specific prevalence and uncertainty maps
## -------------------------------------------------------------------------

## Define panel labels following the ordering used in the combined manuscript
## figure:
##
##   (a) BKP predicted prevalence;
##   (b) BKP prediction uncertainty;
##   (c) LGP predicted prevalence;
##   (d) LGP prediction uncertainty.
sub_labels <- c(
  "(a)",
  "(b)",
  "(c)",
  "(d)"
)

sub_idx <- 1

## Initialize lists for storing the probability and variance plots.
all_probs <- list()
all_vars <- list()

for (method in names(raster_preds)) {
  rp <- raster_preds[[method]]
  
  ## Construct a plotting data frame containing the prediction-grid
  ## coordinates, predictive means, and predictive variances.
  df_method <- data.frame(
    lon = grid_land$lon,
    lat = grid_land$lat,
    mean = rp$mean,
    variance = rp$variance
  )
  
  ## -----------------------------------------------------------------------
  ## Predicted prevalence map
  ## -----------------------------------------------------------------------
  
  p_prob <- ggplot() +
    geom_sf(
      data = africa,
      fill = grey(0.8),
      colour = "gray60"
    ) +
    geom_tile(
      data = df_method,
      aes(
        x = lon,
        y = lat,
        fill = mean
      )
    ) +
    scale_fill_viridis_c(
      option = "C",
      name = "Predicted prevalence",
      limits = c(0, plim),
      breaks = seq(
        0,
        plim,
        0.2
      ),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = grid::unit(
          0.2,
          "cm"
        ),
        barheight = grid::unit(
          4,
          "cm"
        ),
        frame.colour = NA,
        ticks.colour = NA,
        ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_sf(
      xlim = Xbounds[1, ],
      ylim = Xbounds[2, ],
      expand = FALSE
    ) +
    labs(
      title = paste(
        sub_labels[sub_idx],
        method,
        "- Predicted prevalence"
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.key.width = grid::unit(
        0.2,
        "cm"
      ),
      legend.key.height = grid::unit(
        4,
        "cm"
      ),
      legend.ticks.length = grid::unit(
        0.2,
        "cm"
      ),
      legend.title = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        size = 10
      )
    )
  
  sub_idx <- sub_idx + 1
  
  ## Individual model panels can be saved by uncommenting the following lines.
  # ggsave(
  #   sprintf(
  #     "code/figure/Loaloa_%s_prob.pdf",
  #     tolower(method)
  #   ),
  #   p_prob,
  #   width = 6.5,
  #   height = 3.8,
  #   dpi = 300
  # )
  
  
  ## -----------------------------------------------------------------------
  ## Predictive variance map
  ## -----------------------------------------------------------------------
  
  p_var <- ggplot() +
    geom_sf(
      data = africa,
      fill = grey(0.8),
      colour = "gray60"
    ) +
    geom_tile(
      data = df_method,
      aes(
        x = lon,
        y = lat,
        fill = variance
      )
    ) +
    scale_fill_viridis_c(
      option = "C",
      name = "Variance",
      limits = c(0, vlim),
      breaks = seq(
        0,
        vlim,
        0.1
      ),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = grid::unit(
          0.2,
          "cm"
        ),
        barheight = grid::unit(
          4,
          "cm"
        ),
        frame.colour = NA,
        ticks.colour = NA,
        ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_sf(
      xlim = Xbounds[1, ],
      ylim = Xbounds[2, ],
      expand = FALSE
    ) +
    labs(
      title = paste(
        sub_labels[sub_idx],
        method,
        "- Prediction uncertainty"
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.key.width = grid::unit(
        0.2,
        "cm"
      ),
      legend.key.height = grid::unit(
        4,
        "cm"
      ),
      legend.ticks.length = grid::unit(
        0.2,
        "cm"
      ),
      legend.title = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        size = 10
      )
    )
  
  sub_idx <- sub_idx + 1
  
  ## Individual variance panels can be saved by uncommenting the following
  ## lines.
  # ggsave(
  #   sprintf(
  #     "code/figure/Loaloa_%s_var.pdf",
  #     tolower(method)
  #   ),
  #   p_var,
  #   width = 6.5,
  #   height = 3.8,
  #   dpi = 300
  # )
  
  ## Store the model-specific plots for final figure assembly.
  all_probs[[method]] <- p_prob
  all_vars[[method]] <- p_var
}


## -------------------------------------------------------------------------
## Combined manuscript figure
## -------------------------------------------------------------------------

## Arrange the BKP and LGP predictive mean and variance maps into the 2-by-2
## panel layout used in the manuscript.
cat("\n===== Generating combined map (2x2) =====\n")

p_combined <- grid.arrange(
  all_probs$BKP,
  all_vars$BKP,
  all_probs$LGP,
  all_vars$LGP,
  ncol = 2,
  nrow = 2
)

## Save the combined prevalence and uncertainty figure.
ggsave(
  "code/figure/Loaloa_combined.pdf",
  plot = p_combined,
  width = 13,
  height = 9,
  dpi = 300
)

cat("\n===== Analysis complete =====\n")