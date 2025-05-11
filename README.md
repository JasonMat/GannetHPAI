# GannetHPAI
Code and data for HPAI Gannet metapopulation model

The files in this repository can be used to recompile the paper, or to re-run the analysis from the start. 

# Compiling the paper (AISpreadXX.Rmd)
There is an element of the analysis that takes place as part of the knitting of the markdown document attached here. This mostly involves visualisation of aspects of the raw data and post-processing of the sample from the joint posterior of states and parameters. The markdown file draws information from the folders `Data` and `Images`. 

# Re-running the analysis (Real dataXX_X(Comment).R)
This is the JAGS code for the current version of the model. It takes several days to run, and exports 
