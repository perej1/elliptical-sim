# elliptical-sim

Simulation study for an article about multivariate extreme quantile region
estimation. Current simulation settings serve three different purposes:

1. Performance of an elliptical extreme quantile region estimator is compared to
a competing extreme quantile region estimator based on halfspace depth [1].
Relative errors for both estimators are computed for each simulation scenario.
However, simulation data is not stored in Github since files are too large. On
the other hand, directory `summmary-data/` includes some summary statistics for
relative errors of each scenario. Also, illustrative figures are included.

2. An example with skewed t-distribution is constructed. Here we use elliptical
extreme quantile region estimator for estimation, and thus, estimated quantile
regions have an elliptical shape, even though true quantile regions are not
elliptically shaped. For the skewed t-distribution we repeat each scenario once.
Thus, summary statistics are not relevant here, but illustrative figures are in
the directory `summary-data/figures-skew/`.

3. Performance of the elliptical extreme quantile region estimator is examined
   in high dimensional cases. Conservative estimates of the relative errors are stored in the directory `high-dim-data/data/`. On the other hand, figures about how the relative error fluctuates when dimension changes are stored in the directory `high-dim-data/figures/`.

## Requirements

1. Access to Aalto University Triton computing cluster is required for running
   the simulations. See [this link](https://scicomp.aalto.fi/triton/) for
   details about Triton.

2. Clone or unzip the repository.
    ```
    git clone https://github.com/perej1/elliptical-sim.git
    ```

3. Install required packages by running the following R command in the project's
   root folder (R package `renv` has to be installed).
    ```
    renv::restore()
    ```

## Running the simulation

In total there are eight scripts.

- `simulate.R` - Performs simulations for selected parameters.
- `sim-batch.slurm` - Performs simulations for sets of different parameters.
- `summarise.R` - Computes summary statistics and produces figures.
- `functions.R` - Includes functions needed for simulations.
- `gen-arg.R` - Generates arguments for simulation settings.
- `simulate-high-dim.R` - Standalone script for performing simulations in high
  dimensional cases.
- `test-compute_error.R` - Unit tests for the function `compute_error.R`.
- `test-compute_error_elliptical.R` - Unit tests for the function
  `compute_error_elliptical.R`.

With the following command one can run all simulation settings specified by the
`sim-args.txt`.

```
sbatch sim-batch.slurm
```

Instead one can run just one specified scenario (for this Triton is not needed).
```
Rscript simulate.R --type tdistDeg4 --s 100 --d 3 --n 1000 --p high --k medium --seed 278
```

Also, for high dimensional simulations triton is not needed. One can simply run
the following.
```
Rcsript simulate-high-dim.R
```

## Acknowledgements

We acknowledge the computational resources provided by the Aalto Science-IT
project.

## References

[1] Y. He, J. H. Einmahl, Estimation of extreme depth-based quantile regions,
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 79
(2017) 449–461.