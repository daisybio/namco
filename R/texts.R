###### Error messages ######

fileNotFoundError = "The specified file could not be found; please check the file-path and try again."
otuFileNotFoundError = "The specified OTU-file could not be found; please check the file-path and try again."
metaFileNotFoundError = "The specified meta-file could not be found; please check the file-path and try again."
treeFileNotFoundError = "The specified phylogenetic tree file could not be found; please check the file-path and try again."
taxaFileNotFoundError = "The specified taxonomic classification file could not be found; please check the file-path and try again."
biomFileNotFounrError = "The specified biom file could not be found; please check the file-path and try again."
noTaxaInOtuError = "Did not find taxonomy column in the provided OTU file"
wrongTaxaColumnError = "There is a mistake in your taxonomy column; please check that you have no white-spaces inside and not more or less than 6 ; in each row. The wrong row(s) are: "
otuNoMatchTaxaError = "The names of the OTUs in the provided OTU-file and taxonomy file did not correspond or were in a different order; please adapt your files."
treeFileLoadError = "Could not load tree file; check format."
fileEmptyError = "The specified file was empty; please choose different file and try again."
otuOrMetaMissingError = "Missing path to OTU-file or meta-file."
biomMissingError = "Missing path to biom file."
unequalSamplesError = "OTU-file and meta-file do not contain the same samples; please check for consistent spelling."
duplicateSessionNameError = "The chosen project name is already used; please choose a different name for your project."

inconsistentColumnsForest = "Could not find all variables which were used to build model in the columns of new sample file. Please check for consistent spelling."

noEqualFastqPairsError = "Did not find foreward and reverse fastq file for each sample. Please check your input files again!"
differentSampleNamesFastqError = "The names of the samples in your meta file & the sample names of your fastq-files are not the same or one has more/less then the other! Please check again"


###### Info Texts #######

phyloseqSourceText = HTML(paste0("<b>Phyloseq: </b> P. McMurdie, S. Holmes. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. 2013. <a href=\' https://doi.org/10.1371/journal.pone.0061217/\'>  https://doi.org/10.1371/journal.pone.0061217 </a>."))
rheaSourceText = HTML(paste0("<b>RHEA</b>: Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 <a href=\'https://doi.org/10.7717/peerj.2836\'>https://doi.org/10.7717/peerj.2836</a>"))

welcomeText = HTML(paste0("<h3> Welcome to <i>namco</i>, the free Microbiome Explorer</h3>",
                          "<img src=\"Logo.png\" alt=\"Logo\" width=400 height=400>"))

authorsText = HTML(paste0("<b>Authors of this tool:</b>",
                          "Alexander Dietrich (<b>Contact</b> for Issues: ga89noj@mytum.de),",
                          "  [Benjamin Ölke, ",
                          "Maximilian Zwiebel] <br> ",
                          "Supervisor: Monica Matchado, Dr. Markus List, Prof. Dr. Jan Baumbach <br>",
                          "Chair of Experimental Bioinformatics, TU München <br>"))

inputDataText = HTML(paste0("<p>Namco has 2 options to upload microbiome-data:</p>
                <p><span style='text-decoration: underline;'>1) Option 1: OTU-Table and Meta-File:</span>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>1.1) OTU-Table:</span> tab-seperated table, where rows represent OTUs amd columns represent samples. Additionally one column in the file can include the <strong>taxonomic information</strong> for the corresponding OTU of that row. This column must be named <b> taxonomy </b>. <br> The order of taxonomies is: <em>Kingdom, Phylum, Class, Order, Family, Genus and Species</em>, seperated by semicolon. If information for any level is missing the space has to be kept empty, but still, the semicolon has to be present. For an OTU with only taxonomic information of the kingdom the entry would look like this: <code>Bacteria;;;;;;</code></p><p>Namco expects un-normalized input data and allows for normalization during file upload; this can be switched off in the upload window if the user has already normalized data.</p></div>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>1.1) Taxonomic Classification File:</span> tab-seperated file, with the taxonomic classification for each OTU. One column corresponds to a taxonomic level, the name of each column has to be as described in 1). The first column-name can be empty or something like \'taxa\'; the entries in this column are the OTU names, they have to correspond to those in the OTU-table and the phylogenetic tree file (if provided). In total this file has to have 8 columns.</p></div>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>1.2) Meta-file:</span> tab-seperated table with meta information for each sample, mapping at least one categorical variable to those. The first column has to be <b>identical</b> with the column names of the OTU file and has to named <b>SampleID</b></p></div>
                <p><span style='text-decoration: underline;'>2) Option 2: BIOM file:</span></p>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>2.1) combined .biom file:</span>This formate combines OTU-data, meta-data and taxonomic data in one file. See documentation here: <a href=https://biom-format.org/>BIOM</a></p></div>
                <p><span style='text-decoration: underline;'>optional) Phylogenetic tree:</span> To access the full functionalities provided by namco in addition to the OTU table and themapping file, we require a (rooted) phylogenetic tree of representative sequences from each OTU <strong>in Newick format</strong>. Providing an OTU tree is optional, however required to calculate certain measures of beta diversity for instance.</p>"))

testdataText =  HTML(paste0("<p>The testdata was taken from the following experiment: https://onlinelibrary.wiley.com/doi/abs/10.1002/mnfr.201500775. It investigates intestinal barrier integrity in diet induced obese mice. </p>"))

rarefactionInfoText = HTML(paste0("The rarefaction curve allows for analysis of species richness in your un-normalized dataset. Usually these curves grow rapidly at first and then plateau. </br>
                                  This means, that the more common species in a sample are found at first and the more rare ones are found only after sampling many species.</p>
                                  The samples with the steepest slopes overall are highlighted in red (you can choose how many you want to highlight); these samples can then be exlucded from further analysis."))

dimReductionInfoText = HTML(paste0("Here are three methods to reduce the high dimensionality of your OTU table into a low-dimensional space (2D or 3D in this case). The provided methods are: PCA, tSNE and UMAP. </br>
                                   <b>PCA</b>(Principal Component Analysis): Maximizing the variance of the OTU abundace values between samples in lower dimension. This results in clustering similar samples next to each other</br>
                                   <b>tSNE</b>(T-distributed Stochastic Neighbor Embedding): Modelling of a high dimensional datapoint into a lower dimension, such that similar points are modelled by nearby points with high probability</br>
                                   <b>UMAP</b>(Uniform manifold approximation and projection): Similar modelling to tSNE, but assuming that data is uniformly distributed on a locally connected Riemannian manifold</p>
                                   With each method you have the option to color samples by meta groups.</br>"))

confoundingInfoText = HTML(paste0("This tab allows you to find out confounding factors for each of your meta variables. Simply choose a variable of interest and check the result table. <b>[needs phylogenetic tree file to work]</b></br>
                                  The table tells you which of the other variables is considered a confounding factor and if that result is significant (p-value < 0.05). </p>
                                  In the lower half you find the explained variation of each meta variable, meaning which meta variables account for most of the measured variance."))

alphaDivText = HTML(paste0("Alpha-diversity allows to measure the diversity of species inside the samples. Here you can choose between 3 approaches to calculate this value: </p>
                           <b>Shannon-Index:</b> assumes all species are represented in a sample and that they are randomly sampled</br>
                           <b>Simspon-Index:</b> gives more weight to common or dominant species.  In this case, a few rare species with only a few representatives will not affect the diversity.</br>
                           These indices are not linear, meaning a shannon index of x is not twice as diverse as a shannon index of 2x. To account for that, the <b>effective</b> index can be calculated, which correspond to the number of equally abundant species that would yield the
same index value. (<a href=https://esajournals.onlinelibrary.wiley.com/doi/10.1890/06-1736.1> Jost 2007</a>) </br>
                           <b>Richness:</b> simply the summed up occurrence of species per sample (which strongly depends on sequencing depth of (unnormalized) samples)</p>
                           See the detailed formulas of calculation at the bottom of the page."))

alphaDivFormulas = withMathJax(paste0("For sample j:\
                                      $$Richness_j = \\sum_{i \\in OTU} I(x_{ij}>0)$$ \
                                      $$Shannon-Index_j = \\sum_{i \\in OTU} p_{ij} \\cdot \\ln p_{ij}$$ \
                                      $$effective \\; Shannon-Index_j = \\exp(Shannon-Index_j)$$ \
                                      $$Simpson-Index_j=\\sum_{i \\in OTU} p_{ij}^2 $$ \
                                      $$effective \\; Simpson-Index_j=\\frac{1}{Simpson-Index_j}$$ \
                                      with p_{ij} is the relative abundance of OTU i in sample j"))

heatmapText = "Generate a ecology-oriented heatmap with different options of distance calculation. Choose ordination method for organization of rows and columns and distance method for cell values. <b>[needs pyhlogenetic tree file to work]</b>"
heatmapText2 = HTML(paste0(phyloseqSourceText,"<br> For a detailed explaination of the phyloseq heatmap approach see: <a href=\'https://joey711.github.io/phyloseq/plot_heatmap-examples.html\'> Phyloseq-heatmaps </a>"))
heatmapOrdinationText = HTML(paste0("Types of ordination methods:<br>",
                                    "<b>NMDS:</b>  Non-metric MultiDimenstional Scaling <br>",
                                    "<b>MDS/PCoA:</b>  principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) <br>",
                                    "<b>DPCoA:</b>  Double Principle Coordinate Analysis <br>",
                                    "<b>DCA:</b>  detrended correspondence analysis  <br>",
                                    "<b>CCA:</b>  correspondence analysis, or optionally, constrained correspondence analysis <br>",
                                    "<b>RDA:</b>  redundancy analysis <br>"))

picrust2Text = HTML(paste0("<b>PICRUSt2</b> (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) is a software for predicting functional abundances based only on marker gene sequences. <br>",
                           "\"Function\" usually refers to gene families such as KEGG orthologs and Enzyme Classification numbers, but predictions can be made for any arbitrary trait. Similarly, predictions are typically based on 16S rRNA gene sequencing data, but other marker genes can also be used. <br>"))

picrust2SourceText = HTML(paste0("<b>PICRUSt2</b>: ",
                                 "Gavin M. Douglas, Vincent J. Maffei, Jesse Zaneveld, Svetlana N. Yurgel, James R. Brown, Christopher M. Taylor, Curtis Huttenhower, Morgan G. I. Langille, <b>2020</b> <br>",
                                 "<a href=https://www.biorxiv.org/content/10.1101/672295v2> PICRUSt2: An improved and customizable approach for metagenome inference </a>"))


coOcurrenceDistrText = HTML(paste0("This shows the logarithmic distribution in the normalized OTU table (black line is currently selected cutoff)"))

coOcurrenceHeatmapText = HTML(paste0("Heatmap of cutoff-effect: dark fields are being set to 0 in co-occurrence calculation. Y are samples, X are OTUs"))

coOcurrenceCutoffTitleText = HTML(paste0("<h4><b>Cutoff Assignment:<sup>1</sup></b></h4>"))

coOcurrenceCutoffText = HTML(paste0("<sup>1</sup>: All values in the normalized OTU-table, which fall below the chosen cutoff are being set to 0. Those OTUs are then counted as <i> not present </i>. "))

coOcurrenceCountsTitleText = HTML(paste0("<h4><b> Configure Count Calculation:<sup>2</sup></b></h4>"))

coOcurrenceCountsText = HTML(paste0("<sup>2</sup>: Two ways of calculating the counts are possible: <br>",
                                  "Both methods start by counting the co-occurrences of OTUs in all possible sample pairings. Those co-occurrences are then added up to generate a count table for each OTU-pair.",
                                  "By chosing a group from the meta file, this process is executed seperatly for all samples in the group corresponing to a unique covariate. The two tables are then compared: <br>",
                                  "<b> difference</b>: For each OTU pair x and y, calculate: counts(x) - counts(y), where x is the first occuring covariate. <br>",
                                  "<b> log <sub>2</sub> fold-change</b>: For each OTU pair x and y calculate: log<sub>2</sub>(counts(x)+0.001 / counts(y)+0.001), where x is the first occuring covariate."))

spiecEasiSourceText = HTML(paste0("<b>SPIEC-EASI: </b> Z. D. Kurtz, C. Müller, E. Miraldi, D. Littman, M. Blaser and R. Bonneau. Sparse and Compositionally Robust Inference of Microbial Ecological Networks. 2015. <a href=\'https://doi.org/10.1371/journal.pcbi.1004226\'>  https://doi.org/10.1371/journal.pcbi.1004226</a>"))

#### themetagenomics-texts ####

themetagenomicsSourceText = HTML(paste0("<b>themetagenomics</b>: Stephen Woloszynek, Joshua Chang Mell, Gideon Simpson, and Gail Rosen. Exploring thematic structure in 16S rRNA marker gene surveys. 2017. bioRxiv 146126; doi:<a href=\'https://doi.org/10.1101/146126\'>https://doi.org/10.1101/146126</a>"))

themetagenomicsTextTitle = HTML(paste0("<h5>Explore clustering by functional topics in your dataset! Choose covariate of interest to measure its relationship with the samples over topics distribution from the STM. </h5> ",
                                       "<br>", themetagenomicsSourceText))

themetagenomicsTextTopics = HTML(paste0("Pick the number of functional clusters you want to split your OTUs into"))

themetagenomicsTextSigma = HTML(paste0("This sets the strength of regularization towards a diagonalized covariance matrix. Setting the value above 0 can be useful if topics are becoming too highly correlated. Default is 0"))

themetagenomicsText2= HTML(paste0('Below shows topic-to-topic correlations from the samples over topics distribution. The edges represent positive',
                                  ' correlation between two topics, with the size of the edge reflecting to the magnitude of the correlation.',
                                  ' The size of the nodes are consistent with the ordination figure, reflecting the marginal topic frequencies.'))

fastqQualityText = HTML(paste0('This plot shows the base quality at each position. <br>
                                 The grey heatmap is the frequency of quality at each base position. The green line shows the mean quality score; the orange lines show the quartiles of the quality score distribution.'))


#themetagenomicsText3 = HTML(paste0("Integrates the samples over topics p(s|k) and topics over taxa p(k|t) distributions from the STM, the topic correlations from the p(s|k) component, the covariate effects from the p(s|k) component, and their relationship with the raw taxonomic abundances. The covariate effects for each topic are shown as a scatterplot of posterior weights with error bars corresponding the global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression weights with uncertainty intervals are shown. This is analogous to the mean difference between factor levels in the posterior predictive distribution. For continuous covariates, the points again represent the mean regression weights (i.e., the posterior slope estimate of the covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of the standard deviation of the predicted topic probabilities from the posterior predictive distribution. Colors indicate whether a given point was positive (red) or negative (blue) and did not enclose 0 at a user defined uncertainty interval.",
#                                   "The ordination figure maintains the color coding just decribed. The ordination is performed on p(k|t) via either PCoA (using either Jensen-Shannon, Euclidean, Hellinger, Bray-Curtis, Jaccard, or Chi-squared distance) or t-SNE. The latter iterates through decreasing perplexity values (starting at 30) until the algorithm succeeds. The top 2 or 3 axes can be shown. The radius of the topic points corresponds to the topic frequencies marginalized over taxa.",
#                                   "The bar plot behaves in accordance with LDAvis. When no topics are chosen, the overall taxa frequencies are shown. These frequencies do not equal the abundances found in the initial abundance table. Instead, they show p(k|t) multiplied by the marginal topic distribution (in counts). To determine the initial order in which taxa are shown, these two distributions are compared via Kullback-Liebler divergence and then weighted by the overall taxa frequency. The coloration of the bars indiciates the taxonomic group the inidividual taxa belong to. The groups shown are determined based on the abundance of that group in the raw abundance table. When a topic is selected, the relative frequency of a given taxa in that topic is shown in red.",
#                                   "λ controls relevance of taxa within a topic, which in turn is used to adjust the order in which the taxa are shown when a topic is selected. Relevence is essentially a weighted sum between the probability of taxa in a given topic and the probability of taxa in a given topic relative to the overall frequency of that taxa. Adjusting λ influences the relative weighting such that",
#                                   "r = λ x log p(t|k) + λ x log p(t|k)/p(x)",
#                                   "The correlation graph shows the topic correlations from p(s|k) ~ MVN(mu,sigma). Again, the coloration described above is conserved. The size of the nodes reflects the magnitude of the covariate posterior regression weight, whereas the width of the edges represents the value of the positive correlation between the connected nodes. By default, the graph estimates are determined using the the huge package, which first performs a nonparanormal transformation of p(s|k), followed by a Meinhuasen and Buhlman procedure. Alternatively, by choosing the simple method, the correlations are simply a thresholded MAP estimate of p(s|k)."))

#### combined texts ####


sourcesText = HTML(paste0("This tool was built using source-code from: <br> ",
                          rheaSourceText ,"<br>",
                          themetagenomicsSourceText,"<br>",
                          "<br>",
                          "The following packages were used for big parts of the calculations: <br>",
                          phyloseqSourceText, "<br>",
                          spiecEasiSourceText,"<br>",
                          picrust2SourceText,"<br>",
                          "<br> A full list of used packages will be provieded here soon..."))


