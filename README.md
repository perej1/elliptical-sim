# elliptical-sim

Simulated data examples for an article about multivariate extreme quantile
region estimation. Simulations produce figures of the following form.

![image](figures/fig-n_5000-k_400-p_low-gamma_1.00.png)

In above figure $Q_p$ denotes the true $(1-p)$-quantile region, $\bar Q_p$
denotes estimate based on sample quantile and $\hat Q_p$ denotes estimate based
on the extreme quantile region estimator proposed in the article.

## Prequisites

R packages needed for running the simulations can be found on file `renv.lock`.
Notice that the package *extreme* can be only installed from the corresponding
Github repository with
```
devtools::install_github("perej1/extreme")
```

## Running the simulation

Run the following on bash command line
```
bash plot-batch.sh
```

As an result a folder `figures` will be generated (if it does not already
exist). The folder includes generated figures.
