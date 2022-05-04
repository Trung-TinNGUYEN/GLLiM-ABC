# F. Forbes March 2022
#####################################################################################
##
## GLLiM-ABC with SMC using the winference package
##
#####################################################################################

# The ABC-SMC implementation requires package winference from
# https://github.com/pierrejacob/winference
# To install winference on macos (Catalina):
# 1) download the files from github: is a zip file
# 2) unzip and re-buid the package to get mac binaries doing in the term where the
# files are: R CMD build winference
# This create a tar.gz file that can now be install from R-studio

# 3) However, winference requires to pre-install other packages, eg Rmixmod, and 
# others, see error messages...
# 4) Requires also installation of CGAL: use homebrew, macport did not work for me

# Example of R scripts are eg in dir inst/reproduceabc, eg dir mvnormal for 
# the normal location model, ex of plots etc...

# see script files: follows the winference way....

# 1) Normal Location Model (NML) example: mvnormal_wsmc_gllim.R, My_mvnormal_plots.R etc....

# 2) Bivariate Beta : bbm_wsmc_gllim.R , model_get_bbm.R 

# 3) MA(2) : ma2150_wsmc_gllim.R , model_get_ma2.R 
# see also script MA2MSE-SMC.R for new MSE computation with SMC-ABC

# 3) Hyperboloid example: itdmix_smc_gllim.R, model_get_itdmix.R

### Misc: other usefull functions/scripts

# get_gllim_colors.R
# CompDistWABC.R


