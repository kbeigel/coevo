# coevo
Functions for running tests of coevolution (cophylogenetic signal) to assess coevolution in symbiotic systems.


## Install the package using `devtools::install_github()`
```{r}
devtools::install_github('kbeigel/coevo')
```

```{r}
library(coevo)
```

## Example and usage
Please see [kbeigel/Tseptentrionalis_Specificity_2025](https://github.com/kbeigel/Tseptentrionalis_Specificity_2025) for the usage of this mini-package on a real-world dataset.

## Acknowledgements
The underlying functions for this mini-package are built using methods developed by Alix Matthews (@alixmatthews) whose original project code iz available at: [alixmatthews/cophylogenetic](https://github.com/alixmatthews/cophylogenetic) and [alixmatthews/parulidae_proctophyllodidae_cophylo/04_CophylogeneticTests/ParaFit_PACo_Clean.R](https://github.com/alixmatthews/parulidae_proctophyllodidae_cophylo/blob/main/04_CophylogeneticTests/ParaFit_PACo_Clean.R). Katherine Beigel (@kbeigel) developed this mini-package to streamline analysis operations for ParaFit and PACo, improve figure creation (adding bootstraps from original phylos), and put data into table formats for easier writing and downstream manipulation.