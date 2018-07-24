source("http://www.bioconductor.org/biocLite.R")
bio="C://Users//user//Documents//R//win-library//3.5//bioconductor"


#### change environment
R.home() # see where R is installed
.libPaths()

# add to the library paths (temporarily)
.libPaths(c(.libPaths(), bio)) 
.libPaths()

#biocLite()  
biocLite(lib=bio)
?biocLite


#check if all bioconductor packages are the most updated version
#biocValid() 
biocValid(lib=bio)
?biocValid

# Update previously installed Bioconductor or
# CRAN packages and their dependencies.
biocLite("BiocUpgrade") 


pkgs=row.names(installed.packages(lib=bio))
biocLite(pkgs, lib=bio)      


###### install datasets from source

install.packages(path_to_file, repos = NULL, type="source")

##### look at packages
help.start()  # go to packages --> choose a package of interest (eg GRanges)-->
              # user guide, package, vignettes --> option to click slides,PDF

biocLite("GenomicRanges")  # to install a package

library(GenomicRanges, lib.loc=bio)
??GenomicRanges
example(GenomicRanges)  # sample use for the package
help(package="GenomicRanges") # info about the package

browseVignettes(package="GenomicRanges")

### to use latest version of BioConductor
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

# if you get a message like : BiocInstaller version 3.2 is too old for R version 3.3
# do the following

# quit your R session
# start a new session
# run the command  remove.packages("BiocInstaller")
# repeat the command untill R says there is no such package
# run the command source("http://bioconductor.org/biocLite.R")
# run biocValid()



### changing libPaths

#### change environment
R.home() # see where R is installed
libPaths # see where you packages are installed


.libPaths(c(.libPaths(), dir)) # add to the library paths (temporarily)
.libPaths()


#### *** take note that all these changes are only temporarily valid ***####

#to permenantly change your libPath you must do the following

R.home()
.libPaths()

# go to the "etc" subdirectory where R is installed
# 
# open a txt file>>
# put the line

# R_LIBS=C://Users//user//~   (the new lib path you want)
# R_LIBS= include the old lib paths
# R_LIBS= included the old lib paths
# 
# save file as Renviron.site  in the "etc" subdirectory

######## ** if Biobase cannot be installed *** ################
# add to .Rprofile file

local({
  old <- getOption("defaultPackages")
  options(defaultPackages = c(old, "Biobase"))
})