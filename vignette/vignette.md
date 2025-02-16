# coevo: Evaluating coevolution through tests of cophylogenetic signal (ParaFit and PACo)

## Install the package using `devtools::install_github()`
```{r}
devtools::install_github('kbeigel/coevo')
```

```{r}
library(coevo)
```


## Read in host-symbiont association matrix

Rownames should be the host taxa and column names should be the symbiont taxa, e.g.:

```
        symbiont_1  symbiont_2  symbiont_3  symbiont_4  symbiont_5
host_1	        1	        0	        0	        0	        0
host_2	        0	        1	        0	        0	        0
host_3	        0	        0	        1	        0	        0
host_4	        0	        0	        0	        1	        0
host_5	        0	        0	        0	        0	        1	
```


```{r}
hs_assoc = read_associations('testdata/hs_matrix.csv')
```


## Read in host and symbiont trees
```{r}
host_phylo_file = 'testdata/host_tree.nexus'
h = read_phylo(
    filepath = host_phylo_file,
    nexus_format = 'raxml', 
    taxa_to_keep = hs_assoc$host_taxa,
    label_fix = c('SamplePrefix_', '')
)

symbiont_phylo_file = 'testdata/symbiont_tree.nexus'
s = read_phylo(
    filepath = symbiont_phylo_file,
    nexus_format = 'raxml', 
    taxa_to_keep = hs_assoc$symbiont_taxa,
    label_fix = c('\\.', '')
)
```

## Run ParaFit analysis

```{r}
parafit_results = parafit_analysis(
    host_tree = h,
    symbiont_tree = s,
    assoc_mat = hs_assoc,
    n_runs = 100,
    correction_method = 'lingoes'
)
```

## Run PACo analysis
```{r}
paco_results = paco_analysis(
    host_tree = h,
    symbiont_tree = s,
    assoc_mat = hs_assoc,
    n_perm = 100000
)
```

## Compile the results from ParaFit and PACo
```{r}
# Global results
global_results = compile_global_results(
    parafit_results, paco_results
)

# Individual link tests
link_results = compile_link_results(
    parafit_results,
    paco_results,
    p_alpha = 0.05
)
```

## Plotting the cophylo
```{r}
# Plot the individual link results on cophylo
plot_cophylo(
    h_phy = h,
    s_phy = s,
    results = link_results,
    filename = 'cophylo',
    plot_h = 8, 
    plot_w = 14,
    palette = c('gold', 'magenta3', 'cyan3')
)
```