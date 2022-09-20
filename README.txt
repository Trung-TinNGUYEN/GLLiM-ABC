
######################################################################################
# August 2021
#####################################################################################

1) Rfunctions directory:
  
The whole set of R functions to implement the GLLiM-ABC procedures are gathered
in the Rfunctions directory: one file per function 

2) Rscripts directory:
  
This directory contains one R script for each of the example in the paper. To re-run 
one of this example, first load all the required packages as listed in the begining of
file ExampleScript.R . Then, for eg the ITD example, follow the steps in ITD.R. But beware 
that this will create new data sets and target observation and may not then produce exactly
the same results than in the paper. Then go back to the ExampleScript.R file and follow the steps 
there describing the ITD example. Note that several examples are described in this file,
in some order that corresponds to an earlier version of the paper. 
Examples corresponding to the use of SMC-ABC are described in files ExampleScript-V2.R and 
ExampleScript-SMC.R (see update below). 

3) data_et_al4examples directory
Alternatively, to save time or to reproduce exactly the paper's results and plots,
the data sets and target observations used in the paper have been saved  as R objects 
in directory data_et_al4examples in .Rdata files with explicit names, eg dataITDdf3.Rdata.

To reproduce the paper's examples exactly then, go directly to file ExampleScript.R
and follow the steps there. 

#####################################################################################
# UPDATE: March 2022
#####################################################################################
##
## GLLiM-ABC with SMC using the winference package
##
#####################################################################################

 The ABC-SMC implementation requires package winference from
 https://github.com/pierrejacob/winference
 To install winference on macos (Catalina):
 1) download the files from github: is a zip file
 2) unzip and re-buid the package to get mac binaries doing in the term where the
 files are: R CMD build winference
 This create a tar.gz file that can now be install from R-studio

 3) However, winference requires to pre-install other packages, eg Rmixmod, and 
 others, see error messages...
 4) Requires also installation of CGAL: use homebrew, macport did not work for me

 Example of R scripts are eg in dir inst/reproduceabc, eg dir mvnormal for 
 the normal location model, ex of plots etc...

## see  the following script files: follows the winference way....

 1) Normal Location Model (NML) example: mvnormal_wsmc_gllim.R, My_mvnormal_plots.R etc....

 2) Bivariate Beta : bbm_wsmc_gllim.R , model_get_bbm.R 

 3) MA(2) : ma2150_wsmc_gllim.R , model_get_ma2.R 
 see also script MA2MSE-SMC.R for new MSE computation with SMC-ABC

 3) Hyperboloid example: itdmix_smc_gllim.R, model_get_itdmix.R

### Misc: other useful functions/scripts

 get_gllim_colors.R
 CompDistWABC.R



