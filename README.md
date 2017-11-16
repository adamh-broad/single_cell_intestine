# A single-cell survey of the intestinal epithelium

## Abstract
Intestinal epithelial cells absorb nutrients, respond to microbes, function as a barrier and help to coordinate immune responses. Here we report profiling of 53,193 individual epithelial cells from the small intestine and organoids of mice, which enabled the identification and characterization of previously unknown subtypes of intestinal epithelial cell and their gene signatures. We found unexpected diversity in hormone-secreting enteroendocrine cells and constructed the taxonomy of newly identified subtypes, and distinguished between two subtypes of tuft cell, one of which expresses the epithelial cytokine Tslp and the pan-immune marker CD45, which was not previously associated with non-haematopoietic cells. We also characterized the ways in which cell-intrinsic states and the proportions of different cell types respond to bacterial and helminth infections: Salmonella infection caused an increase in the abundance of Paneth cells and enterocytes, and broad activation of an antimicrobial program; Heligmosomoides polygyrus caused an increase in the abundance of goblet and tuft cells. Our survey highlights previously unidentified markers and programs, associates sensory molecules with cell types, and uncovers principles of gut homeostasis and response to pathogens.

Experimental workflow            |  Infection with H.polygyrus
:-------------------------:|:-------------------------:
![](https://github.com/adamh-broad/single_cell_intestine/blob/master/Fig1a.jpg)  |  ![](https://github.com/adamh-broad/single_cell_intestine/blob/master/Relmb.jpg)

## Single-cell analysis code
<a href="https://github.com/adamh-broad/single_cell_intestine/blob/master/Analysis.md">Quick start</a> R code for analysis of 7,216 intestinal epithelial cells. It's a good idea to download Fxns.R and atlas_tsne.txt first. There are several R packages needed to run this code, each can be install using the 'install.packages' command. For example, to install the NMF package used to render a heatmap, run the R command 'install.packages(NMF)'.

## Related Resources
* <a href="https://www.nature.com/articles/nature24489">Our paper</a>
* <a href="https://portals.broadinstitute.org/single_cell/study/small-intestinal-epithelium">Single Cell Portal (Broad Institute)</a>
