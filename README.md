# psyp23rerps
Code and data to reproduce the analyses reported by Aurnhammer et al., 2023, Psychophysiology.

# Code

The ```code/``` directory contains collection of functions implementing and running the regression analyses. Run the ```do_*``` scripts, to run the main analyses.
Regressions are computed in ```julia```, visualisation and multiple comparisons correction are implemented in ```R```. ImageMagick shell commands are used to combine graphs into the final figures as published in the article.

Assuming an alias called ```julia``` and that data are placed in the ```data/``` directory, running the following commands from the ```code/``` replicates the analyses.

```julia do_lmerRT.jl```

```julia do_rERP.jl```

```Rscript stimuli_densities.r```

```Rscript do_plot_lmerRT.r```

```Rscript do_plot_rERP.r```

```cd ../plots/Figures```
```sh ../plots/Figures/mk_figures.sh```


# Data

The data are provided as a release within this repository and are to be placed in the ```data/``` directory.

# Language and package versions

**```julia v1.7.3```**

```DataFrames v1.3.4```
```Combinatorics v1.0.2```
```CSV v0.10.4```
```Distributions v0.25.62```
```StatsBase v0.33.17```
```LinearAlgebra (standard)```
```CategoricalArrays v0.10.6```
```MixedModels v4.7.0```
```PooledArrays v1.4.2```

**```R v4.1.2```**

```data.table v1.14.2```
```ggplot2 v3.3.5```
```gridExtra v2.3```

**```zsh v5.8```**

```ImageMagick v7.1.0-62```
