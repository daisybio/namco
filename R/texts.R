###### Error messages ######

fileNotFoundError = "The specified file could not be found; please check the file-path and try again."
otuFileNotFoundError = "The specified OTU-file could not be found; please check the file-path and try again."
metaFileNotFoundError = "The specified meta-file could not be found; please check the file-path and try again."
treeFileNotFoundError = "The specified phylogenetic tree file could not be found; please check the file-path and try again."
taxaFileNotFoundError = "The specified taxonomic classification file could not be found; please check the file-path and try again."
noTaxaInOtuError = "Did not find taxonomy column in the provided OTU file"
otuNoMatchTaxaError = "The names of the OTUs in the provided OTU-file and taxonomy file did not correspond or were in a different order; please adapt your files."
treeFileLoadError = "Could not load tree file; check format."
fileEmptyError = "The specified file was empty; please choose different file and try again."
otuOrMetaMissingError = "Missing path to OTU-file or meta-file."
unequalSamplesError = "OTU-file and meta-file do not contain the same samples; please check for consistent spelling."
duplicateSessionNameError = "The chosen project name is already used; please choose a different name for your project."

inconsistentColumnsForest = "Could not find all variables which were used to build model in the columns of new sample file. Please check for consistent spelling."




###### Info Texts #######

heatmapText = "Generate a ecology-oriented heatmap with different options of distance calculation. Choose ordination method for organization of rows and columns and distance method for cell values."
heatmapOrdinationText = HTML(paste0("Types of ordination methods:<br>",
                                    "<b>NMDS:</b>  Non-metric MultiDimenstional Scaling <br>",
                                    "<b>MDS/PCoA:</b>  principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) <br>",
                                    "<b>DPCoA:</b>  Double Principle Coordinate Analysis <br>",
                                    "<b>DCA:</b>  detrended correspondence analysis  <br>",
                                    "<b>CCA:</b>  correspondence analysis, or optionally, constrained correspondence analysis <br>",
                                    "<b>RDA:</b>  redundancy analysis <br>"))


