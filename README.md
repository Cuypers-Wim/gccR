# gccR
gene co-expression conservation calculations in R

# Functionality
- Expression conservation (EC) calculations between two gene expression compendia
- Identification and analysis of functional expression classes (FECs)

# ICC Method Overview

## Input
The ICC (Iterative Comparison of gene Co-expression) method starts with two gene expression matrices, M1 and M2. These matrices typically represent expression data for two distinct species. Both M1 and M2 should include measurements for multiple genes (rows) across various conditions (columns). In our experiments, log fold changes were used, where each expression value signifies a change relative to a "control" condition.

## Expression Conservation (EC) Calculation
We employ the *Iterative Comparison of gene Co-expression (ICC)* procedure for EC calculation, and implemented these steps:

1. **Symmetrical Correlation Matrices**: Obtain symmetrical correlation matrices for each strain.
2. **Orthologous Gene Pair Alignment**: Extract submatrices comprising genes with orthologous counterparts in the other strain and reorder them. This ensures that a given row 'n' in both M1 and M2 corresponds to an orthologous gene pair.
3. **Weighted Pearson Correlation**: We next perform an iterative calculation of weighted Pearson correlation for each orthologous gene pair to derive the final EC score. The weight vector is determined by the previous set of correlation values, with negative values set to '0' to minimize the impact of highly divergent gene pairs.
4. **Calculation Details**: Utilize the `wtd.cors` function from the weights package v1.0.4 for the weighted correlation calculations. Repeat the calculations until convergence is achieved, defined as a difference of less than 0.0001 (default value, adaptable) between a given EC value and the EC value of the previous iteration.

<p align="center">
  <img src="gccr.png" alt="Figure 1">
</p>

## EC Background Distributions
Distinguishing actual biological divergence in co-expression conservation from technical noise involves estimating two distributions.

1. **Conserved Distribution**: _Calculate how low EC values could be due to technical variation in the scenario of perfect conservation_. This distribution considers variation introduced by different conditions and is calculated by splitting the largest gene expression compendium and performing the ICC procedure. Optionally, the procedure can be repeated ten times on different data splits, considering experiment grouping. **Genes with EC values lower than the minimum of this distribution have likely not conserved their co-expression profile across strains/species.**
2. **Diverged Distribution**: _Calculate the maximum EC value attainable when there is no correlation between two gene co-expression profiles_. This involves permutating the expression values of only one gene per iteration, resulting in a conservative random distribution. **Genes with EC values higher than the maximum of this distribution likely have conserved their gene co-expression profile across strains/species.**

# Installation

## From R (studio)

``` r
# Installation is currently only possible via the development version from GitHub:
devtools::install_github("Cuypers-Wim/gccR")

```

## From a linux command line using conda

Experiencing issues with R package incompatibilities or complex installation procedures? Don't worry, we've got you covered. 
Whether you're running into errors with package conflicts or seeking a streamlined approach to deploy gccR on a server environment (such as a Linux server for parallel function execution), 
this Conda environment setup will simplify the process.

First install all packages required for gccR:

``` 
conda env create -f environment.yml
conda activate env_name

```

Then run this line of code in the environment to install gccR

``` 
Rscript -e "devtools::install_github('Cuypers-Wim/gccR')"
``` 

Now you're good to go and use gccR interactive or via a script.

# Example usage

```R
library(gccR)

corM1 <- get_corM(sTyph_LT2, dropNArows = TRUE)
corM2 <- get_corM(sTyph_SL1344, dropNArows = TRUE)

# extract submatrices

corM_ortho <- extract_core_submatrix(corM1, corM2, singleCopyOrthologs)

csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, as.matrix(singleCopyOrthologs), rename = TRUE)

# compute the expression conservation scores by means of iterative comparison of co-expression

EC <- getEC(corM_ortho$csM1, csM2_ordered, 0.001)

# estimate background distributions

perfectEC1  <- perfect_EC(sTyph_LT2, csM2_ordered, conv = 0.001, maxIter = 200, threads = 1)
randomEC1 <- divergedEC(sTyph_LT2, corM_ortho$csM1, csM2_ordered, EC$ECweights)

distr_plot <- plot_distributions(EC$ECfinal, perfectEC1$EC_final, randomEC1)
distr_plot

# Obtain 

# Obtain vectors of diverged and conserved genes --------------------------

distributions_list <- list(EC = EC$ECfinal, conserved = perfectEC1$EC_final, diverged = randomEC1)
names(distributions_list)

colnames_df <- c()
min_val <- c()
max_val <- c()

for (i in 1:length(distributions_list)) {
  
  colnames_df[i] <- names(distributions_list)[i]
  
  min_val[i] <- min(distributions_list[[i]])
  max_val[i] <- max(distributions_list[[i]])
  
  
}

# inspect the results

results_df <- cbind(names = colnames_df, min = min_val, max = max_val)

# diverged genes and their EC values

EC$ECfinal[which(EC$ECfinal < min(perfectEC1$EC_final))]

# conserved genes and their EC values

EC$ECfinal[which(EC$ECfinal > max(randomEC1))]


```
# Credits
## Methodology adapted from:

- Meysman P, Sánchez-Rodríguez A, Fu Q, Marchal K, Engelen K. Expression divergence between Escherichia coli and Salmonella enterica serovar Typhimurium reflects their lifestyles. Mol Biol Evol. 2013 Jun;30(6):1302-14. doi: 10.1093/molbev/mst029. Epub 2013 Feb 20. PMID: 23427276; PMCID: PMC3649669.

- Sonego, P. et al. (2015). Comparative Analysis of Gene Expression: Uncovering Expression Conservation and Divergence Between Salmonella enterica Serovar Typhimurium Strains LT2 and 14028S. In: Mengoni, A., Galardini, M., Fondi, M. (eds) Bacterial Pangenomics. Methods in Molecular Biology, vol 1231. Humana Press, New York, NY. https://doi.org/10.1007/978-1-4939-1720-4_8

 ## The Salmonella data included in this repository was retrieved from the COLOMBOS database (colombos.net):
 
- Moretto M, Sonego P, Dierckxsens N, Brilli M, Bianco L, Ledezma-Tejeida D, Gama-Castro S, Galardini M, Romualdi C, Laukens K, Collado-Vides J, Meysman P, Engelen K. COLOMBOS v3.0: leveraging gene expression compendia for cross-species analyses. Nucleic Acids Res. 2016 Jan 4;44(D1):D620-3. doi: 10.1093/nar/gkv1251. Epub 2015 Nov 19. PMID: 26586805; PMCID: PMC4702885.
