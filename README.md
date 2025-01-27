This repository contains all the codes used in the following paper titled

- Kanrar, R., Li, C., Ghodsi, Z., Gamalo, M. (2025+). Risk-inclusive Contextual Bandits for Early Phase Clinical Trials.

## Summary of Contents

- `/code/function`: contains all the supporting functions developed to conduct empirical experiments.
- `/code/real_data`: contains scripts to clean, visualize and analyze the real data utilized in the article.
- `/code/simulation`: contains scripts to conduct all simulation experiments described in both the main manuscript and the supplementary material. 
- `/code/visualization`: contains scripts that generates all figures except Figure 1 and Figure 2 in the main manuscript.
- `/code/metadata`: Additional .RData files created to run the simulation experiments conveniently. Files contain lists to specify hyper-parameters required in the simulation experiments.


## Initial Setup:

- Clone Github Repository:


```
git clone git@github.com:rohitkanrar/RiTS.git
cd RiTS
```
- Create additional folders to save output files:

```
mkdir output
```

- Install all required R packages:

```
source("code/r/requirements.R")
```

## A Test Case:

Here, we include steps to replicate the simulation experiment conducted in Section 3.

- Run trials:
```
source("code/simulation/replicated_trial.R")

```
- Visualize results by running the script `code/visualization/replicated_trial_figures.R`.
- Generate tables by running the script `-code/visualize/replicated_trial_tables.R`.


Please feel free to email at rohitk@iastate.edu with questions or to report any bugs. 
