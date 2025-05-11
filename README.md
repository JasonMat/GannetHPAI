# GannetHPAI
Code and data for HPAI Gannet metapopulation model

The files in this repository can be used to recompile the paper, or to re-run the analysis from the start. 

# Compiling the paper (AISpreadXX.Rmd)
There is an element of the analysis that takes place as part of the knitting of the markdown document attached here. This mostly involves visualisation of aspects of the raw data and post-processing of the sample from the joint posterior of states and parameters. The markdown file draws information from the folders `Data` and `Images`. It also requires posterior MCMC outputs from the folder `Output`, which can be found at this location `https://gla-my.sharepoint.com/:f:/g/personal/jason_matthiopoulos_glasgow_ac_uk/Ep_DahDPQlBOrvyvYZllky4BC2uKPkIs-fMsDWUSLF772Q?e=PCakkD` (Approx. 50-100Mb size). The typical install and runtime for this compilation is about an hour. Markdown chunks have the `Cache` option switched on, so follow-on compilations will be faster. 

# Re-running the analysis (Real dataXX_X(Comment).R)
This is the JAGS code for the current version of the model. It takes several days to run, and exports an MCMC coda file to the folder `Output`, which is too large to include in the GitHub repository. It can be found here: `https://gla-my.sharepoint.com/:f:/g/personal/jason_matthiopoulos_glasgow_ac_uk/Ep_DahDPQlBOrvyvYZllky4BC2uKPkIs-fMsDWUSLF772Q?e=PCakkD` (Approx. 50-100Mb size). The analysis can be Demo-ed by setting `n` (the number of colonies) to a small number (e.g. 3), but the results of such analyses are not sensible since the epidemiological dynamics of the population are affected by the complete network of colonies. 

# Software 
R version 4.4.3 (2025-02-28 ucrt) -- "Trophy Case"

JAGS 4.3. 2. 

R Libraries required for Markdown document:
library(boot)

library(coda)

library(runjags)

library(maps)

library(rlang)

library(maps)

library(ggplot2)

library(runjags)

library(coda)

library(stringr)

library(terra)

library(ggplot2)

library(rnaturalearth)

library(kableExtra)

library(png)

library(tidyr)

library(tidyterra)

library(grid)

library(scales)

library(patchwork)

library(ggimage)

library(ggspatial)

library(RCurl)

library(sf)


R Libraries required for JAGS modelling

library(runjags)

library(coda)

library(stringr)
