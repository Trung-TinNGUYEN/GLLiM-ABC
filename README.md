# GLLiM-ABC

This repository contains all the numerical experiments from "Approximate Bayesian computation with surrogate posteriors" https://hal.archives-ouvertes.fr/hal-03139256.

1) Rfunctions directory:
The whole set of R functions to implement the GLLiM-ABC procedures are gathered in the Rfunctions directory: one file per function.

2) Rscripts directory:
This directory contains one R script for each of the example in the paper. To re-run  one of this example, eg the ITD example follow the steps in ITD.R. But beware that this will create new data sets and target observation and may not then produce exactly the same results than in the paper. 
Then go the ExampleScript.R file and follow the steps there. Note that examples are described in the order of the paper but all in the same file. 

3) data_et_al4examples directory
Alternatively, to save time or to reproduce exactly the paper's results and plots, the data sets and target observations used in the paper have been saved  as R objects in directory data_et_al4examples in .Rdata files with explicit names, eg dataITDdf3.Rdata.

To reproduce the paper's examples exactly then, go directly to file ExampleScript.R
and follow the steps there. 
