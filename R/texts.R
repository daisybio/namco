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

welcomeText = HTML(paste0("<h3> Welcome to <i>namco</i>, the free Microbiome Explorer</h3>",
                          "<img src=\"Logo.png\" alt=\"Logo\" width=400 height=400>"))

authorsText = HTML(paste0("<b>Authors of this tool:</b>",
                          "Alexander Dietrich (Contact for Issues: ga89noj@mytum.de)",
                          "Benjamin Ölke, ",
                          "Maximilian Zwiebel <br> ",
                          "Supervisor: Monica Matchado, Dr. Markus List, Prof. Dr. Jan Baumbach <br>",
                          "Experimental Chair of Bioinformatics, TU München <br>"))

sourcesText = HTML(paste0("This tool was built using source-code from <br> ",
                          "<b>RHEA</b>: Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836 <br>",
                          "<b>themetagenomics</b>: Stephen Woloszynek, Joshua Chang Mell, Gideon Simpson, and Gail Rosen. Exploring thematic structure in 16S rRNA marker gene surveys. 2017. bioRxiv 146126; doi: https://doi.org/10.1101/146126 <br>",
                          "A full list of used packages will be provieded here soon..."))

inputDataText = HTML(paste0("<p>Namco has 2 obligatory input files &amp; and one optional file:</p>
                <p><span style='text-decoration: underline;'>1) OTU-Table:</span> tab-seperated table, where rows represent OTUs amd columns represent samples. Additionally one column in the file can include the <strong>taxonomic information</strong> for the corresponding OTU of that row. This column must be named <b> taxonomy </b>. <br> The order of taxonomies is: <em>Kingdom, Phylum, Class, Order, Family, Genus and Species</em>, seperated by semicolon. If information for any level is missing the space has to be kept empty, but still, the semicolon has to be present. For an OTU with only taxonomic information of the kingdom the entry would look like this: <code>Bacteria;;;;;;</code></p><p>Namco expects un-normalized input data and allows for normalization during file upload; this can be switched off in the upload window if the user has already normalized data.</p>
                <p><span style='text-decoration: underline;'>1.1) Taxonomic Classification File:</span> tab-seperated file, with the taxonomic classification for each OTU. One column corresponds to a taxonomic level, the name of each column has to be as described in 1). The first column-name can be empty or something like \'taxa\'; the entries in this column are the OTU names, they have to correspond to those in the OTU-table and the phylogenetic tree file (if provided). In total this file has to have 8 columns.</p>
                <p><span style='text-decoration: underline;'>2) Meta-file:</span> tab-seperated table with meta information for each sample, mapping at least one categorical variable to those. The first column has to be <b>identical</b> with the column names of the OTU file and has to named <b>SampleID</b></p>
                <p><span style='text-decoration: underline;'>3) Phylogenetic tree:</span> To access the full functionalities provided by namco in addition to the OTU table and themapping file, we require a (rooted) phylogenetic tree of representative sequences from each OTU <strong>in Newick format</strong>. Providing an OTU tree is optional, however required to calculate certain measures of beta diversity for instance.</p>"))

testdataText =  HTML(paste0("<p>The testdata was taken from the following experiment: https://onlinelibrary.wiley.com/doi/abs/10.1002/mnfr.201500775. It investigates intestinal barrier integrity in diet induced obese mice. </p>"))

heatmapText = "Generate a ecology-oriented heatmap with different options of distance calculation. Choose ordination method for organization of rows and columns and distance method for cell values."
heatmapOrdinationText = HTML(paste0("Types of ordination methods:<br>",
                                    "<b>NMDS:</b>  Non-metric MultiDimenstional Scaling <br>",
                                    "<b>MDS/PCoA:</b>  principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) <br>",
                                    "<b>DPCoA:</b>  Double Principle Coordinate Analysis <br>",
                                    "<b>DCA:</b>  detrended correspondence analysis  <br>",
                                    "<b>CCA:</b>  correspondence analysis, or optionally, constrained correspondence analysis <br>",
                                    "<b>RDA:</b>  redundancy analysis <br>"))

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

#### themetagenomics-texts ####

themetagenomicsTextTitle = HTML(paste0("<h5>Explore clustering by functional topics in your dataset! Choose covariate of interest to measure its relationship with the samples over topics distribution from the STM. </h5> <i>(for detailed explanation of tool scroll to bottom of page)</i>"))

themetagenomicsTextTopics = HTML(paste0("Pick the number of functional clusters you want to split your OTUs into"))

themetagenomicsTextSigma = HTML(paste0("This sets the strength of regularization towards a diagonalized covariance matrix. Setting the value above 0 can be useful if topics are becoming too highly correlated. Default is 0"))

themetagenomicsText2= HTML(paste0('Below shows topic-to-topic correlations from the samples over topics distribution. The edges represent positive',
                                  ' correlation between two topics, with the size of the edge reflecting to the magnitude of the correlation.',
                                  ' The size of the nodes are consistent with the ordination figure, reflecting the marginal topic frequencies.'))


themetagenomicsText3 = HTML(paste0("Integrates the samples over topics p(s|k) and topics over taxa p(k|t) distributions from the STM, the topic correlations from the p(s|k) component, the covariate effects from the p(s|k) component, and their relationship with the raw taxonomic abundances. The covariate effects for each topic are shown as a scatterplot of posterior weights with error bars corresponding the global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression weights with uncertainty intervals are shown. This is analogous to the mean difference between factor levels in the posterior predictive distribution. For continuous covariates, the points again represent the mean regression weights (i.e., the posterior slope estimate of the covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of the standard deviation of the predicted topic probabilities from the posterior predictive distribution. Colors indicate whether a given point was positive (red) or negative (blue) and did not enclose 0 at a user defined uncertainty interval.",
                                   "The ordination figure maintains the color coding just decribed. The ordination is performed on p(k|t) via either PCoA (using either Jensen-Shannon, Euclidean, Hellinger, Bray-Curtis, Jaccard, or Chi-squared distance) or t-SNE. The latter iterates through decreasing perplexity values (starting at 30) until the algorithm succeeds. The top 2 or 3 axes can be shown. The radius of the topic points corresponds to the topic frequencies marginalized over taxa.",
                                   "The bar plot behaves in accordance with LDAvis. When no topics are chosen, the overall taxa frequencies are shown. These frequencies do not equal the abundances found in the initial abundance table. Instead, they show p(k|t) multiplied by the marginal topic distribution (in counts). To determine the initial order in which taxa are shown, these two distributions are compared via Kullback-Liebler divergence and then weighted by the overall taxa frequency. The coloration of the bars indiciates the taxonomic group the inidividual taxa belong to. The groups shown are determined based on the abundance of that group in the raw abundance table. When a topic is selected, the relative frequency of a given taxa in that topic is shown in red.",
                                   "λ controls relevance of taxa within a topic, which in turn is used to adjust the order in which the taxa are shown when a topic is selected. Relevence is essentially a weighted sum between the probability of taxa in a given topic and the probability of taxa in a given topic relative to the overall frequency of that taxa. Adjusting λ influences the relative weighting such that",
                                   "r = λ x log p(t|k) + λ x log p(t|k)/p(x)",
                                   "The correlation graph shows the topic correlations from p(s|k) ~ MVN(mu,sigma). Again, the coloration described above is conserved. The size of the nodes reflects the magnitude of the covariate posterior regression weight, whereas the width of the edges represents the value of the positive correlation between the connected nodes. By default, the graph estimates are determined using the the huge package, which first performs a nonparanormal transformation of p(s|k), followed by a Meinhuasen and Buhlman procedure. Alternatively, by choosing the simple method, the correlations are simply a thresholded MAP estimate of p(s|k)."))


