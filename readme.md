# Simulation Study for Adaptively-sampled Spatial Capture-Recapture (AS-SCR)

## Overview

Adaptively-sampled SCR is an extension of ordinary SCR whereby the selection of survey units are conditional on observations of an index variable made on the first visit to that survey unit exceeding some *a priori* threshold. The key advantage of this method is that sites where organisms are *not* present may be avoided, while still accounting for the preferential sampling induced by the adaptive method. The potential to reduce the effort, cost, and time scale with the rarity of the species. Nevertheless, the researcher must be aware of the assumptions of this model, the primary one being fulfillment of a valid statistical relationship between the index variable and the local population size that is constant throughout the survey period. Other assumptions include mutual independence of the populations sampled by each survey implements. All of the assumptions of ordinary SCR also apply here.

The repository contains R scripts that simulate spatial capture-recapture data under three different sampling schemes (adaptively sampled, fully sampled, and simple-random sampled), and analyze the data under 5 different models representing the three sampling schemes plus configurations including index values (+ post-fixes), and ignoring index values (- post-fixes), see model configurations below.

## Model configurations

### Model 1 
*AS-SCR+* is the proposed adaptive sampling procedure that we test for bias and precision. We expect to observe unbiased parameters estimates, particularly estimated population size. 

### Model 2
*AS-SCR--* represents an preferential sampling situation and we expect it should be positively biased for density, because it samples only high-density areas without integrating information obtained from indices that fall below the index threshold. 

### Model 3
*F-SCR+* is the ideal and most costly sampling procedure, incorporating SCR and index observations at every site, and it represents the most accurate and precise baseline to which we compare *AS-SCR+*. 

### Model 4
*F-SCR--* represents a standard SCR procedure implemented at all sites without the index data, which we expect will resemble *F-SCR+* with a slight loss in precision. This is an ordinary application of SCR against which we test our new procedure.

### Model 5
*SRS-SCR--* represents application of standard SCR at a simple random sample of sites equal to the corresponding number of sites sampled by AS-SCR. It makes no use of index data. This comparison is important to include because it is a sampling method which a cost-restrained survey could adopt having little initial information regarding the local population. We expect this will be unbiased, but also that it will have the poorest precision in comparison to the other models in the set due to the relatively sparse information it gathers.

## Simulation scripts

### Population simulation

The script to perform the simulation of the populations, and simulate captures is `AS_dataset_sim_v4.7.R`. There are a number of default values so that the function to perform the simulation is `AS.simulator()` may be run without arguments. This function has no side-effects, instead opting to generate data libraries from which the analysis scripts will draw.

The population is simulated simply as a Poisson random variable at `R` total populations (can be defined in the arguments; the default is `R = 10` for easy viewing of the data). The index variable is simulated as a Poisson random variable with the mean equal to `cov.mult * N[r]`, where `r` is one of the `R` populations; `cov.mult = 3` by default.

Based upon the population matrix `N`, animal activity centers are generated on a space of dimension:

```
xlim<- c(-1, 4)
ylim<- c(-1, 5)
```

having uniform random positions.

### Capture simulation 

The simulation for the captures are more involved, and have three main sections:

* The full sampling procedure (not indicated by comments, follows immediately after activity center simulation), used by Model 3 and Model 4, where captures are simulated at all `R` sites.
* The subset sampling procedure (indicated by "BEGIN SUBSET DATA SECTION"), which simulates captures at a subset of the `R` sites, conditional on the index variable, used by Model 1 and Model 2.
* The SRS sampling procedure (indicated by "BEGIN MODEL 5 DATA SECTION"), which simulates captures at a subset of the `R` sites, equal to the number of sites sampled by the corresponding adaptive sample simulation, but selected randomly (**not** conditional on index variable)

#### Output files

The output files are .Rdata files, and are referenced next by the data augmentation script.

### Augmentation script

The analytical procedure uses Bayesian estimation of parameters within JAGS, and data augmentation is performed upon the capture histories to provide stable dimensions for the matrix of individual capture histories. A series of all-zero capture histories are appended to the observed data to represent "potential" individuals. A proportion of those existed, but went unobserved due to imperfect detection, and the rest are "structural", imaginary individuals. Data augmented analysis introduces a new paramter psi that estimates how many of these "potential" individuals are the ones that existed. 

The script `AS_augmentor_v3.1.R` takes the outputs of the `AS.simulator()` function, and augments the data, providing datasets ready to be analyzed within JAGS under the models specified. 

### Model scripts

There are five model scripts, one for each of the models described above in the Model configurations section. They are represented simply by `model1.txt`, `model2.txt`, and so on.

### Analysis scripts

Again, there are five scripts for performing the analysis of the data under each sampling scheme. These are the files following the pattern `MODEL1_rerun.R`, `MODEL2_rerun.R`, `MODEL3_rerun.R`, etc. They automatically extract the data from the .Rdata libraries, and perform the JAGS analysis using the corresponding model .txt file. The script features ordinary parallelization using `doParallel`, and will generate a log indicating what analyses have been completed in the event of unexpected shutdown. 