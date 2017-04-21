# BGA Genomic Annotations Lab Pre-Lab

## Install Software

* Install R: [https://cran.r-project.org/](https://cran.r-project.org/)
* Install RStudio: [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)
* Install packages needed for this lab (in R): `install.packages(pkgs = c('corpcor', 'magrittr', 'stringr', 'readr', 'devtools', 'tidyverse'))` or 
    * If you're running this from an R console on Sphere, type `declare -x R_LIBS_USER="/home/student20/R/library”` in bash before starting R.
* Install this package: `devtools::install_github('andrewGhazi/genomicAnnotationsLab')`

Test that you've got the package installed correctly by loading the package with `library(genomicAnnotationsLab)` and running `?makeScore` to look at the help documentation for a function we will be using.  

There's nothing to turn in for the pre-lab.

# Lab

You can see the most up-to-date version of the lab [here](http://htmlpreview.github.io/?https://github.com/andrewGhazi/genomicAnnotationsLab/blob/master/bgaGenomicAnnotationsLab.html).
