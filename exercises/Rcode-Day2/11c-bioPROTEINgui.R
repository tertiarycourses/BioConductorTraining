source("http://www.bioconductor.org/biocLite.R")
bio="C://Users//user//Documents//R//win-library//3.4//bioconductor"


#### change environment
R.home() # see where R is installed
.libPaths()

# add to the library paths (temporarily)
.libPaths(c(.libPaths(), bio)) 
.libPaths()

#biocLite()  
biocLite(lib=bio)

############################## Protein GUIs ####################################

# pRolocGUI - Interactive visualisation of spatial proteomics data
 biocLite(c("pRoloc","pRolocdata"), lib=bio)
 library(pRolocdata)   # datasets
 pRolocdata()    # see list of data > to load type data(dataset name) 
 data()#> put name of dataset in () to load into environment
 # ?datasetname > for more info
 library(pRolocGUI)
 pRolocVis() # put name of dataset in ()


# diGeR - analysing 2d DIGE data (differential protein spot)
# install.packages("C://Users//user//Desktop//tertiary//BioConductor(advance)//GUIs//digeR_1.3.tar.gz", type="source", repos=NULL, lib=bio)
#library(digeR)
#digeR()

# MSGFgui - investigate MS-GF+ identification protein data
library(MSGFgui)
MSGFgui()
