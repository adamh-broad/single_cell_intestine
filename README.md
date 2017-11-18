# A single-cell survey of the intestinal epithelium

## Single-cell analysis code
This repo contains <a href="https://github.com/adamh-broad/single_cell_intestine/blob/master/Analysis.md">quick start</a> R code for analysis of 7,216 intestinal epithelial cells. The easiest way to run it is to clone the repo (download ZIP button at top) and open Analysis.Rmd in RStudio. 

There are several R (3.3 or later) packages needed to run the code, each can be install using the 'install.packages' command. For example, to install the NMF package used to render a heatmap, run the R command 'install.packages(NMF)'. The 'sva' package used for batch correction can be installed from Bioconductor:

```{r }
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
``` 

## Abstract
Intestinal epithelial cells absorb nutrients, respond to microbes, function as a barrier and help to coordinate immune responses. Here we report profiling of 53,193 individual epithelial cells from the small intestine and organoids of mice, which enabled the identification and characterization of previously unknown subtypes of intestinal epithelial cell and their gene signatures. We found unexpected diversity in hormone-secreting enteroendocrine cells and constructed the taxonomy of newly identified subtypes, and distinguished between two subtypes of tuft cell, one of which expresses the epithelial cytokine Tslp and the pan-immune marker CD45, which was not previously associated with non-haematopoietic cells. We also characterized the ways in which cell-intrinsic states and the proportions of different cell types respond to bacterial and helminth infections: Salmonella infection caused an increase in the abundance of Paneth cells and enterocytes, and broad activation of an antimicrobial program; Heligmosomoides polygyrus caused an increase in the abundance of goblet and tuft cells. Our survey highlights previously unidentified markers and programs, associates sensory molecules with cell types, and uncovers principles of gut homeostasis and response to pathogens.

Experimental workflow            |  Infection with H.polygyrus
:-------------------------:|:-------------------------:
![](https://github.com/adamh-broad/single_cell_intestine/blob/master/Fig1a.jpg)  |  ![](https://github.com/adamh-broad/single_cell_intestine/blob/master/Relmb.jpg)

## Related Resources
* <a href="https://www.nature.com/articles/nature24489">Our paper</a>
* <a href="https://portals.broadinstitute.org/single_cell/study/small-intestinal-epithelium">Single Cell Portal (Broad Institute)</a>

For questions or issues email:
ahaber -at- broadinstitute.org
