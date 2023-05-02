# elliptical-sim

Simulation study for an article about multivariate extreme quantile
region estimation. Two things are done in the simulations:

1. Performance of an elliptical extreme quantile region estimator is compared to
a competing extreme quantile region estimator based on halfspace depth [1].
Relative errors for both estimators are computed for each simulation scenario,
and are then saved in the directory `sim-data/errors/`. Additionally, directory
`sim-data/` includes the following data for reproducibility:

   - true quantile regions,
   - estimates (elliptical and depth) and
   - samples, from which the estimates are calculated.

   On the other hand, directory `summmary-data/` includes some statistics such
as minimum, median and maximum relative errors for each scenario. Also, figures
corresponding to minimum, median and maximum relative errors are included.

2. An example with clover shaped quantile regions is constructed. Here we use
elliptical extreme quantile region estimator for estimation, and thus, estimated
quantile regions have an elliptical shape, even though true quantile regions are
not elliptically shaped.

## Requirements

1. Clone or unzip the repository.
    ```
    git clone https://github.com/perej1/elliptical-sim.git
    ```

2. Install required packages by running the following R command in the project's
   root folder (R package `renv` has to be installed).
    ```
    renv::restore()
    ```

## Running the simulation

In total there are six R scripts.

- `simulate.R` - Performs simulations for selected parameters.
- `simulate-batch.R` - Performs simulations for sets of different parameters.
  Additionally, this script is responsible for argument parsing.
- `summarise.R` - Computes minimum, median and maximum relative errors.
  Additionally, the script outputs figures corresponding to minimum, median and
  maximum relative errors.
- `functions.R` - Includes all the required packages and functions.
- `test-compute_error.R` - Unit tests for the function `compute_error.R`.

Run the script `simulate-batch.R` with desired arguments. For example, below we run simulations as a whole.

```
Rscript simulate-batch.R --simulate TRUE --summarise TRUE --clover TRUE
```

For example, if you want to exclude clover example, just run the following.
```
Rscript simulate-batch.R --simulate TRUE --summarise TRUE --clover FALSE
```

## References

[1] Y. He, J. H. Einmahl, Estimation of extreme depth-based quantile regions,
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 79
(2017) 449–461.