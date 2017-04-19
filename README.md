# BGA Genomic Annotations Lab Pre-Lab

## Install Software

* Install R: [https://cran.r-project.org/](https://cran.r-project.org/)
* Install RStudio: [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)
* Install packages needed for this lab (in R): `install.packages(pkgs = c('corpcor', 'magrittr', 'stringr', 'readr', 'devtools', 'tidyverse'))` or 
    * if you're doing this day-of and are in a hurry: `install.packages(pkgs = c('corpcor', 'magrittr', 'stringr', 'readr', 'devtools', 'dplyr', 'tidyr', 'ggplot2'))`
    * Or if you're running this from an R console on Sphere, just add `lib.loc = '/home/student20/R/x86_64-redhat-linux-gnu-library/3.2'` to any calls to library() for example `library(genomicAnnotationsLab, lib.loc = '/home/student20/R/x86_64-redhat-linux-gnu-library/3.2')`.
* Install this package: `devtools::install_github('andrewGhazi/genomicAnnotationsLab')`

The lab itself should have been sent out as an HTML document!

# Lab

You can see the most up-to-date version of the lab [here](http://htmlpreview.github.io/?https://github.com/andrewGhazi/genomicAnnotationsLab/blob/master/bgaGenomicAnnotationsLab.html).
