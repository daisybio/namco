###### Error messages ######

fileNotFoundError = "The specified file could not be found; please check the file-path and try again."
otuFileNotFoundError = "The specified OTU-file could not be found; please check the file-path and try again."
metaFileNotFoundError = "The specified meta-file could not be found; please check the file-path and try again."
treeFileNotFoundError = "The specified phylogenetic tree file could not be found; please check the file-path and try again."
taxaFileNotFoundError = "The specified taxonomic classification file could not be found; please check the file-path and try again."
biomFileNotFounrError = "The specified biom file could not be found; please check the file-path and try again."
noTaxaInOtuError = "Did not find taxonomy column in the provided OTU file"
differentSamplesInOtuAndMetaError = "There are different samples present in the meta file and the OTU table."
didNotFindSampleColumnError = "The column that was specified as sample-column could not be found in the provided meta file"
didNotFindSpikeColumnError = "Did not find \"amount_spike\"-column in meta-file"
didNotFindWeightColumnError = "Did not find \"total_weight_in_g\"-column in meta-file"
rmSpikesNoMetaError = "Meta file with specific columns is required to remove spikes (check Info&Settings tab)"
wrongTaxaColumnError = "There is a mistake in your taxonomy column; please check that you have no white-spaces inside and not more or less than 6 ; in each row. The wrong row(s) are: "
otuNoMatchTaxaError = "The names of the OTUs in the provided OTU-file and taxonomy file did not correspond; please adapt your files."
treeFileLoadError = "Could not load tree file; check format."
fileEmptyError = "The specified file was empty; please choose different file and try again."
otuOrMetaMissingError = "Missing path to OTU-file."
biomMissingError = "Missing path to biom file."
unequalSamplesError = "OTU-file and meta-file do not contain the same samples; please check for consistent spelling."
duplicateSessionNameError = "The chosen project name is already used; please choose a different name for your project."
noFileError = "No file was chosen."
changeFileEncodingError = "The uploaded file could not be read due to a unknown file encoding."
noTaxaRemainingAfterFilterError = "No taxa are remaining after filtering; please adapt your filtering parameters."
filteringHadNoEffectError = "The applied filtering did not remove any OTUs."
errorDuringDecompression = "There was an issue with your uploaded file(s): if you used a compressed file, check for correct compression method (.tar, .tar.gz, .zip). Did you upload an even number of fastq-files?"
siamcatNotEnoughSamplesError = "In the chosen label, the selected case-value appears only 5 or less times. These are not enough samples to proceed."
noNumericVariablesError = "There are no numeric variables in your dataset; correlations can not be calculated."
timeAndSampleGroupEqualError = "You cannot select the same group variable as a time-point and sample-group value! Please re-think your chosen values."
timeSeriesEqualVariablesError = "You cannot select the same group variable as time-point, sample-group or mean line. Pleasre-think your chosen values."
picrustDifferentialGroupNotFoundError = "The selected sample group for differential analysis could not be found. Differential Analysis will not be performed."
picrustFilesMissingError = "Something went wrong with picrust2, not all files were created. Did your fastq-files have the correct OTU/ASV-names? 
  If you feel like you did nothing wrong, please contanct the author of namco."

inconsistentColumnsForest = "Could not find all variables which were used to build model in the columns of new sample file. Please check for consistent spelling."

noEqualFastqPairsError = "Did not find foreward and reverse fastq file for each sample. Please check your input files again!"
differentSampleNamesFastqError = "The names of the samples in your meta file & the sample names of your fastq-files are not the same or one has more/less then the other! Please check again"

###### Logging Texts #######

log_startText = paste0("#############", Sys.time(), "#############")


###### Info Texts #######

phyloseqSourceText = HTML(paste0("<b>Phyloseq: </b> P. McMurdie, S. Holmes. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. 2013. <a href=\' https://doi.org/10.1371/journal.pone.0061217/\'>  https://doi.org/10.1371/journal.pone.0061217 </a>."))
rheaSourceText = HTML(paste0("<b>RHEA</b>: Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 <a href=\'https://doi.org/10.7717/peerj.2836\'>https://doi.org/10.7717/peerj.2836</a>"))

welcomeText = HTML(paste0("<h2> Welcome to <i>namco</i>, the free Microbiome Explorer</h2>",
                          "<img src=\"Logo.png\" alt=\"Logo\" width=400 height=400>"))

startHereText = HTML(paste0("<img src=\"left-arrow.png\" alt=\"Logo\" width=75 height=75>",
                            "<h3>Start here!</h3>"))

authorsText = HTML(paste0("<b>Authors of this tool:</b>",
                          "Alexander Dietrich (<b>Contact</b> for Issues: alex.dietrich@tum.de),",
                          "  [Benjamin Ölke, ",
                          "Maximilian Zwiebel] <br> ",
                          "Supervisor: Monica Matchado, Dr. Markus List, Prof. Dr. Jan Baumbach <br>",
                          "Chair of Experimental Bioinformatics, TU München <br>"))

contactText = HTML(paste0("<b>Contact for Issues:</b><br>",
                          "<a href = \"mailto: alex.dietrich@tum.de\">alex.dietrich@tum.de</a><br>",
                          "[<b>NOT:</b> alexander.dietrich@tum.de!, this is not me..]"))

inputDataText = HTML(paste0("<p>Namco has 2 options to upload microbiome-data:</p>
                <p><span style='text-decoration: underline;'><b>1) Option 1: OTU-Table and Meta-File:</b></span>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>1.1) OTU-Table:</span> tab-seperated table, where rows represent OTUs amd columns represent samples. Additionally one column in the file can include the <strong>taxonomic information</strong> for the corresponding OTU of that row. This column must be named <b> taxonomy </b>. <br> The order of taxonomies is: <em>Kingdom, Phylum, Class, Order, Family, Genus and Species</em>, seperated by semicolon. If information for any level is missing the space has to be kept empty, but still, the semicolon has to be present. For an OTU with only taxonomic information of the kingdom the entry would look like this: <code>Bacteria;;;;;;</code></p><p>Namco expects un-normalized input data and allows for normalization during file upload; this can be switched off in the upload window if the user has already normalized data.</p></div>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>1.2) Taxonomic Classification File:</span><b>[only needed if no taxonomy column in OTU-file!]</b> tab-seperated file, with the taxonomic classification for each OTU. One column corresponds to a taxonomic level, the name of each column has to be as described in 1). The first column-name can be empty or something like \'taxa\'; the entries in this column are the OTU names, they have to correspond to those in the OTU-table and the phylogenetic tree file (if provided). In total this file has to have 8 columns.</p></div>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>1.3) Meta-file:</span> tab-seperated table with meta information for each sample, mapping at least one categorical or numerical variable to those. One column has to contain the sample-names; enter the name of this column below the meta-file upload field. The entries of this column have to be <b>identical</b> with the column names of the OTU file.</p></div>
                <p><span style='text-decoration: underline;'>1.4) Phylogenetic tree:</span><b>[optional]</b> To access the full functionalities provided by namco in addition to the OTU table and the mapping file, a (rooted) phylogenetic tree of representative sequences from each OTU <strong>in Newick format</strong> is required. Providing an OTU tree is optional, however required to generate specific plots and analysis (indicated by this logo: ", fontawesome::fa("tree", fill="red", height="1.5em"),").</p>
                <p><span style='text-decoration: underline;'><b>2) Option 2: raw fastq-files:</b></span></p>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>2.1) Multiple fastq files:</span> Select multiple .fastq or .fastq.gz files; for each selected file, another paired file has to be selected (so only a even amount of files can be uploaded). For each file, the Ilumina naming convention is expected: <a href=https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm>Ilumnia Naming convention</a>. <b>Important:</b> The sample names of these files will have to be identical with the sample names of the provided meta-file (if provided)!</p></div>
                <div style='text-indent:25px;'><p><span style='text-decoration: underline;'>2.2) Meta-file:</span><b>[optional]</b> tab-seperated table with meta information for each sample, mapping at least one categorical or numerical variable to those. One column has to contain the sample-names; enter the name of this column below the meta-file upload field. Each entry of this column has to be <b>identical</b> with the sample name of one pair of fastq files!.</p></div>"))

testdataText =  HTML(paste0("<p>The testdata was taken from the following experiment: https://onlinelibrary.wiley.com/doi/abs/10.1002/mnfr.201500775. It investigates intestinal barrier integrity in diet induced obese mice. </p>"))


#####basic texts#####

taxBinningText = HTML(paste0("Taxonomic Binning shows you the distribution of taxonomies in your samples/groups for one specific level (Genus, etc.)</p>",
                             "<u>When to use it?</u> Always, it will provide you with an overview, which taxa are highly present and might be of interest for your studies."))


alphaDivText = HTML(paste("Alpha-diversity refers to the diversity <i>within</i> a sample. Different measures exist on how to calculate this value.</p>",
                          "<u>When to use it?</u> Always, as it gives you a good first overview, which samples/sample groups behave different than others.<br><br>",
                          "<u>Additional Information:</u><br>",
                          "<b>Shannon-Index:</b> assumes all species are represented in a sample and that they are randomly sampled</br>",
                          "<b>Simspon-Index:</b> gives more weight to common or dominant species.  In this case, a few rare species with only a few representatives will not affect the diversity.</br>",
                          "These indices are not linear, meaning a shannon index of x is not twice as diverse as a shannon index of 2x. To account for that, the <b>effective</b> index can be calculated, which correspond to the number of equally abundant species that would yield the same index value. (<a href=https://esajournals.onlinelibrary.wiley.com/doi/10.1890/06-1736.1> Jost 2007</a>) </br>",
                          "<b>Richness:</b> simply the summed up occurrence of species per sample (which strongly depends on sequencing depth of (unnormalized) samples)</p>",
                          "See the detailed formulas of calculation at the bottom of the page.<br> <br>",
                          "Additionally to the alpha-diversity, you can perform a pairwise wilcoxon-test (if the sample group is selected), where you can choose which sub-groups you want to compare to another.<br>",
                          "You can choose to display the p-values of the pairwise test or the corresponding <a href=https://www.statology.org/significance-codes-in-r/>significance codes<a>."))

alphaDivFormulas = withMathJax(paste0("For sample j:\
                                      $$Richness_j = \\sum_{i \\in OTU} I(x_{ij}>0)$$ \
                                      $$Shannon-Index_j = \\sum_{i \\in OTU} p_{ij} \\cdot \\ln p_{ij}$$ \
                                      $$effective \\; Shannon-Index_j = \\exp(Shannon-Index_j)$$ \
                                      $$Simpson-Index_j=\\sum_{i \\in OTU} p_{ij}^2 $$ \
                                      $$effective \\; Simpson-Index_j=\\frac{1}{Simpson-Index_j}$$ \
                                      with p_{ij} is the relative abundance of OTU i in sample j"))

betaDivText = HTML(paste0("Beta-diversity refers to the distance between samples. Different measures exist on how to calculate this distance, some of which need a phylogenetic tree ",fontawesome::fa("tree", fill="red", height="1em")," to be calculated. </p>",
                          "<u>When to use it?</u> Always, it can show you how sample-groups might be clustered together. <br><br>
                          <u>Additional Information:</u><br>
                          <b>Bray-Curtis Dissimilarity:</b> dissimilarity based on counts in each sample </br>
                          <b>Unweighted Uni-Frac Distance:</b> uses phylogenetic distance between samples ",fontawesome::fa("tree", fill="red", height="1em"),"</br>
                          <b>Weighted Uni-Frac Distance:</b> also consideres relative abundance of OTUs ",fontawesome::fa("tree", fill="red", height="1em"),"</br>
                          <b>Generalized Uni-Frac Distance:</b> a balanced version of unweighted and weighted Uni-Frac distance, which avoids sensitivity to rare and highly abundant lineages ",fontawesome::fa("tree", fill="red", height="1em"),"</br>
                          <b>Variance adjusted Uni-Frac Distance:</b> Weighted Uni-Frac Distance, which takes variation of weight into account ",fontawesome::fa("tree", fill="red", height="1em"),"</br><br>",
                          "Additionally, Permutational Multivariate Analysis of Variance (PERMANOVA) is performed with the chosen distance matrix on the specified sample groups; the p-values are indicated in the bottom of the ordination plots. A low p-value (<0.05) indicates that the null hypothesis can be rejected, which would state that the sample groups have the same centroids and are equal."))

dimReductionInfoText = HTML(paste0("Here you can use three different methods to see the difference in variance between your samples.</p>",
                                   "<u>When to use it?</u> If you had unexpected results in the beta-diversity analysis or as a second confirmation on your beta-diversity results. <br><br>",
                                   "<u>Additional Information:</u><br>",
                                   "Here are three methods to reduce the high dimensionality of your OTU table into a low-dimensional space (2D or 3D in this case). The provided methods are: PCA, tSNE and UMAP. </br>
                                   <b>PCA</b>(Principal Component Analysis): Maximizing the variance of the OTU abundace values between samples in lower dimension. This results in clustering similar samples next to each other</br>
                                   <b>tSNE</b>(T-distributed Stochastic Neighbor Embedding): Modelling of a high dimensional datapoint into a lower dimension, such that similar points are modelled by nearby points with high probability</br>
                                   <b>UMAP</b>(Uniform manifold approximation and projection): Similar modelling to tSNE, but assuming that data is uniformly distributed on a locally connected Riemannian manifold</p>
                                   With each method you have the option to color samples by meta groups.</br>"))


rarefactionInfoText = HTML(paste0("The rarefaction curve allows you to find samples, where sequencing depth might be much lower than in the rest of your samples. Usually the more common species in a sample are found at first and then - after more sampling - the rare ones are found. This is indicated by a curve which grows rapidly at first and then plateaus. Curves which stop before a plateau is reached can show you samples, where there might be many unrecognized species. You have the option to remove such samples .</p>",
                                  "<u>When to use it?</u> Always, it shows you samples which are not deeply sequenced."))

heatmapText = HTML(paste0("Here you can find a heatmap of the OTU abundance values. There are different methods of ordering the rows and columns, which are explained below. </p>",
                          "<u>When to use it?</u> As an overview of your abundance values. It can help you to find OTUs / sample groups which are very high/low abundant. <br><br>",
                          "<u>Additional Information:</u><br>",
                          "Generate a ecology-oriented heatmap with different options of distance calculation. Instead of the traditionally used hierarchical clustering, here you can use ordination methods like NMDS to organize rows and columns. <br>",
                          "<b>Orientation/Ordination Method</b>: The ordering of the rows (taxa) and columns (samples) can be changed with the orientation (ordination) method. The different methods are named on the bottom of this page.<br>",
                          "<b>Distance Method</b> : There are multiple options to calculate the distance between samples (for a explaination check the beta-diversity tab). Additionally the Jensen–Shannon divergence (jsd) can be used, which involves approximating the data with a continuous probability density function.<br>",
                          "<b>Axis labeling</b>: You can label the columns (samples) by the different sample-groups which were provided in the meta-file.<br>",
                          "<b>Order Samples</b>: If you choose 'yes', the selected orientation method is ignored for the samples; they will rather be ordered by the levels of the selected group (<i>Axis labeling</i>)<br>",
                          "Note: <b>[needs phylogenetic tree file to work",fontawesome::fa("tree", fill="red", height="1.5em"),"]</b>"))
heatmapText2 = HTML(paste0(phyloseqSourceText,"<br> ---- For a detailed explaination of the phyloseq heatmap approach see: <a href=\'https://joey711.github.io/phyloseq/plot_heatmap-examples.html\'> Phyloseq-heatmaps </a>"))
heatmapOrdinationText = HTML(paste0("<u>Types of ordination methods:</u><br>",
                                    "<b>NMDS:</b>  Non-metric MultiDimenstional Scaling <br>",
                                    "<b>MDS/PCoA:</b>  principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) <br>",
                                    "<b>DPCoA:</b>  Double Principle Coordinate Analysis <br>",
                                    "<b>DCA:</b>  detrended correspondence analysis  <br>",
                                    "<b>CCA:</b>  correspondence analysis, or optionally, constrained correspondence analysis <br>",
                                    "<b>RDA:</b>  redundancy analysis <br>"))
neatmapSourceText = HTML(paste0("<b>NeatMap</b>: NeatMap - non-clustering heat map alternatives in R ",
                                "Satwik Rajaram & Yoshi Oono, <b>2010</b>",
                                "<a href=https://doi.org/10.1186/1471-2105-11-45> https://doi.org/10.1186/1471-2105-11-45</a>"))


#####differential texts#####

associationsText = HTML(paste0("Here you can find out which OTUs or taxonomies differ significantly between two sample groups. Note that sometimes there are no results, because nothing was considered significant; you can increase your significance cutoff but then the results might not be as meaningful of course.</p>",
                               "<u>When to use it?</u> If you want to find out which taxa might be differentially abundant between groups such as healthy vs. sick. It can lead you to possible targets, which might explain such a disease.<br><br>",
                               "<u>Additional Information:</u><br>",
                               "First you choose, in which taxonomic level you want to compare the sample groups, by selecting either <i>OTU</i> or a different level in the first select-box.<br>",
                               "You can then choose in the options, which sample groups you want to compare against another",
                               "(eg. meta-label=<i>GENDER</i> and case-level=<i>FEMALE</i>; this would compare the <i>FEMALE</i> samples against all other samples (eg. <i>MALE</i>).) The chosen case-level will be colored red, the rest blue. If only two sample groups exist for that label, the chosen case-level is red, the other label blue. <br>",
                               "A non-parametric Wilcoxon test is then performed between those 2 groups in order to find significantly different OTUs; you can change the cutoff, below which an OTU is considered significant with the p-value slider. <br>",
                               "Additionally some effect size values can be displayed, like the fold-change (difference in abundance), AUC (range from 0-1: ) and prevalence shift.<br><br>",
                               "If no plot loads, no significant results were found and you might have to adapt your significance niveau."))

associationsSourceText = HTML(paste0("<b>SIAMCAT</b>:",
                                     "Wirbel J, Zych K, Essex M, Karcher N, Kartal E, Salazar G, Bork P, Sunagawa S, Zeller G <b>2020</b><br>",
                                     "<a href=https://doi.org/10.1186/s13059-021-02306-1> SIAMCAT: user-friendly and versatile machine learning workflows for statistically rigorous microbiome analyses</a>"))

corrText = HTML(paste0("Here you can find out which OTUs correlate with each other as well as with <u>numeric</u> meta-variables (such as age, weight, dosage of drug,...). Pearson correlation is measured from -1 to 1, with 1 being a perfect positive correlation (<i>high</i> age infers <i>high</i> abundance of taxa).</p>",
                       "<u>When to use it?</u> If you want to detect which meta variable might have an effect on abundance of OTUs. It can also show which OTUs might always appear together/not together (clusters of high/low correlation in OTU heatmap).<br><br>",
                       "<u>Additional Information:</u><br>", 
                       "Additionally, a <i>significance cutoff</i> can be chosen (default = 0.05): All OTUs which have not a single correlation with a p-value lower than this cutoff are not displayed in the final plot. This means, each pairwise correlation that is displayed has at least one significant correlation somewhere. The non-significant correlations of the displayed OTUs can either be set to <i>blank</i>, <i>highlighted</i> or just <i>do nothing</i>.<br>",
                       "We highly advise you to use some sort of normalization for this analysis step (check the 4 normalization methods on the sidebar at the left).<br>",
                       "Note: entries with no values (NA) or a value of infinity will be replaced with the mean of the corresponding variable (eg. mean age of all samples).<br><br>",
                        rheaSourceText))

themetagenomicsSourceText = HTML(paste0("<b>themetagenomics</b>: Stephen Woloszynek, Joshua Chang Mell, Gideon Simpson, and Gail Rosen. Exploring thematic structure in 16S rRNA marker gene surveys. 2017 <a href=\'https://doi.org/10.1371/journal.pone.0219235\'>https://doi.org/10.1371/journal.pone.0219235</a>, <a href=\'https://cran.rstudio.com/web/packages/themetagenomics/vignettes/thematic_structure.html\'>function-topic interactions blog entry</a>"))

topicText = HTML(paste0("Topic Modelling is an idea from language modelling, where words with similar meaning are being clustered together. The same principle holds for this approach: OTUs which share a similar expression will be clustered together in one “topic”. This is a different approach to the standard clustering algorithms and can therefore bring new insights in OTU similarities and differences.</p>",
                        "<u>When to use it?</u> If you want to find similar partners to specific OTUs or if you just want to explore, which OTUs behave similar. <br><br>",
                        themetagenomicsSourceText))

themetagenomicsTextTitle = HTML(paste0("<h5>Explore clustering by functional topics in your dataset! Choose covariate of interest to measure its relationship with the samples over topics distribution from the STM. </h5> ",
                                       "<br>", themetagenomicsSourceText))

themetagenomicsTextTopics = HTML(paste0("Pick the number of functional clusters you want to split your OTUs into"))

themetagenomicsTextSigma = HTML(paste0("This sets the strength of regularization towards a diagonalized covariance matrix. Setting the value above 0 can be useful if topics are becoming too highly correlated. Default is 0"))

themetagenomicsText2= HTML(paste0('Below shows topic-to-topic correlations from the samples over topics distribution. The edges represent positive',
                                  ' correlation between two topics, with the size of the edge reflecting to the magnitude of the correlation.',
                                  ' The size of the nodes are consistent with the ordination figure, reflecting the marginal topic frequencies.<br>',
                                  "<i>If no plot appears, no positive correlation between topics could be found and therefore no edges exist!</i>"))

#####functional texts#####

picrust2Text = HTML(paste0("Here you can apply the <b>PICRUSt2</b> (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) algorithm on your dataset. It aims to predict the functional annotation of each of the provided OTUs as well as performing a copy-number normalization. Function can be predicted in three classes: enzyme classification (EC), KEGG orthology (KO) and the molecular pathway in which the OTU is present (PW). Additionally, a differential functional analysis is run (for which you have to select a sample group), which helps you to find OTUs which have significantly different functions in the sub-groups.</p>",
                           "<u>When to use it?</u> Always, it is a big additional analysis step, which moves into the field of functional analysis of your taxa. You might find interesting pathways or enzymes in which your dataset has an impact on."))

picrust2SourceText = HTML(paste0("<b>PICRUSt2</b>: ",
                                 "Gavin M. Douglas, Vincent J. Maffei, Jesse Zaneveld, Svetlana N. Yurgel, James R. Brown, Christopher M. Taylor, Curtis Huttenhower, Morgan G. I. Langille, <b>2020</b> <br>",
                                 "<a href=https://www.biorxiv.org/content/10.1101/672295v2> PICRUSt2: An improved and customizable approach for metagenome inference </a>"))

aldexSourceText = HTML(paste0("<b>ALDEx2</b>: ",
                              "Fernandes, AD, Macklaim, JM, Linn, TG, Reid, G, Gloor, GB <b>2013</b><br>",
                              "<a href=http://dx.doi.org/10.1371%2/journal.pone.0067019>NOVA-Like Differential Gene Expression Analysis of Single-Organism and Meta-RNA-Seq</a><br>",
                              "Fernandes, D. A, Reid, J., Macklaim, M. J, McMurrough, T.A, Edgell, D.R., Gloor, B. G <b>2014</b><br>",
                              "<a href=http://doi:10.1186/2049-2618-2-15>Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis.</a><br>",
                              "Gloor GB, Macklaim JM, Fernandes AD <b>2016</b><br>",
                              "<a href=http://dx.doi.org/10.1080/10618600.2015.1131161>Displaying Variation in Large Datasets: a Visual Summary of Effect Sizes</a>"))

picrust_pval_info_text = HTML(paste0("<b>Differential analysis using the ALDEx2 package</b><br>",
                                     "Picrust2 can predict the functions of 16S rRNA data. It can predict 3 different classification types of function:<br>",
                                     "1) <b>EC</b> (=Enzyme Classification): classification number of enzymatic function (<a href=https://www.enzyme-database.org/class.php> Enzyme database</a>)<br>",
                                     "2) <b>KO</b> (=KEGG orthology): molecular functions represented in terms of functional orthologs (<a href=https://www.genome.jp/kegg/ko.html> KEGG </a>)<br>",
                                     "3) <b>PW</b> (=pathway): molecular pathway <br>",
                                     "ALDEx2 then performs several analysis steps ((a) generates Monte Carlo samples of the Dirichlet distribution for each sample, (b) converts each instance using a log-ratio transform, then (c) returns test results for two sample Welch t-test).
                                     The first row of plots helps you to find functions, which are significantly different between the two sample-groups you are comparing. The sigificance value can be chosen manually, default is a BH p-value of 0.05.<br>",
                                     "P-values are displayed raw and Benjamini-Hochberg (BH) transformed, which corrects p-values by avoiding false positives (significance is selected by looking at the BH pvalue).<br>",
                                     "The second row plots displays more detailed information on the functions considered significant."))

#####phylogenetic texts#####


phyloTreeText = HTML(paste0("This tab provides a display of the provided (or calculated) phylogenetic tree. Additionally you can label the OTUs by a chosen taxonomic level or a sample group (for the sample group the labeling will be how often a specific OTU appears in a sub-group).</p>",
                            "<u>When to use it?</u> Always, this step belongs to the “standard” analysis steps to get an overview of the relations between OTUs.<br><br>",
                            "<u>Additional Information:</u><br>",
                            "You can also add heatmaps of taxonomic levels or/and meta-groups. The meta-group heatmaps are colored by occurence, meaning how often a specific OTU appears in the corresponding meta-group. <br>",
                            "Most likely you will have to adapt the positioning of the heatmaps in the <i>advanced options</i>, depending on your chosen groups. <br>",
                            "[The aligning of tips has sometimes a strange behaviour, where the offsets will be too large with multiple groups; we are working on fixing this in the future.]"))


#####network texts#####

coOccurrenceInfoText = HTML(paste0("This is a very basic solution on calculating networks with OTU abundance data. It simply connects those pairs of OTUs which both appear in a specific sub-group (they “co-occur”) with green edges, and those who do not appear together in red. You can set the abundance threshold at the top to choose, when an OTU is considered to be not present.</p>",
                                   "<u>When to use it?</u> To find out which OTUs might be more abundant in specific sample groups.<br><br>",
                                   "<u>Additional Information:</u><br>",
                                   "There are 3 steps to building this network, each can be influenced by parameters:<br>",
                                   "<b>Binarization</b>: This step will create a binary matrix out of the abundance OTU matrix. All abundance values, which fall below a certain cutoff, will bes set to 0. This can be interpreted as an OTU in one sample being <i>not present</i>. All remaning OTUs will be considered <i>present</i>.<br>",
                                   "<b>Count calculation</b>: In this step, the number of pairs of <i>present</i> OTUs will be counted (= co-occuring OTUs). This is done seperatly for two groups of samples (eg. case and control), which can be chosen manually. Then for each pair of OTUs either the difference or log2 fold-change between the two sample groups is calculated.<br>",
                                   "<b>Network</b>: The network is then a representation of OTUs (nodes) and count-values (edges) of co-occurring OTU-pairs. Only the pairs with most extreme count-values will be displayed. OTUs can also be colored by their taxonomic group."))

coOcurrenceDistrText = HTML(paste0("This shows the logarithmic distribution in the normalized OTU table (black line is currently selected cutoff)"))

coOcurrenceHeatmapText = HTML(paste0("Heatmap of cutoff-effect: dark fields are being set to 0 in co-occurrence calculation. Y are samples, X are OTUs"))

coOcurrenceCutoffTitleText = HTML(paste0("<h4><b> I) Cutoff Assignment:</b></h4>"))

coOcurrenceCutoffText = HTML(paste0("All values in the normalized OTU-table, which fall below the chosen cutoff are being set to 0. Those OTUs are then counted as <i> not present </i>. "))

coOcurrenceCountsTitleText = HTML(paste0("<h4><b> II) Count Parameters:</b></h4>"))

coOcurrenceNetworkTitleText = HTML(paste0("<h4><b> III) Network:</b></h4>"))

coOcurrenceCountsText = HTML(paste0("Two ways of calculating the counts are possible: <br>",
                                    "Both methods start by counting the co-occurrences of OTUs in all possible sample pairings. Those co-occurrences are then added up to generate a count table for each OTU-pair.",
                                    "By chosing a group from the meta file, this process is executed seperatly for all samples in the group corresponing to a unique covariate. The two tables are then compared: <br>",
                                    "<b> difference</b>: For each OTU pair x and y, calculate: counts(x) - counts(y), where x is the first occuring covariate. <br>",
                                    "<b> log <sub>2</sub> fold-change</b>: For each OTU pair x and y calculate: log<sub>2</sub>(counts(x)+0.001 / counts(y)+0.001), where x is the first occuring covariate."))

compNetworkInfoText = HTML(paste0("Here you have the option to choose from multiple statistical methods in order to calculate a network from the OTU abundances. Each method has a different approach on calculating edge-weights, therefore each can result in a different network representation.</p>",
                                  "<u>When to use it?</u> Here you can try around with different methods to find groups of OTUs, which form hubs in the network or which are far apart and therefore have very different behaviour. <br><br>",
                                  "<u>Additional Information:</u><br>",
                                  "You can choose the method for network calculation as well as how to deal with 0s, normalization and cluster (groups of OTUs) detection.<br>",
                                  "We advise you to apply some appropriate <b>filtering steps</b> before creating a network, since they can become quite overloaded otherwise."))

diffNetworkInfoText = HTML(paste0("The methods to calculate the network are the same as in the “Network inference” tab, but here you can create two separate networks for sample groups with two subgroups (such as healthy vs sick). You can then compare these networks to find OTUs which have different interactions in the two sub-groups.</p>",
                                  "<u>When to use it?</u> This is another tool to find similarities and differences between two sample groups, but with a maybe more clear visual representation. <br><br>"))

taxNetworkInfoText = HTML(paste0("Same as the “Network inference” tab, only that it is not on OTU level, but rather on a chosen taxonomic level (such as Phylum). This can make the networks clearer to look at and can provide insights to higher level interactions.</p>",
                                 "<u>Additional Information:</u><br>",
                                 "It may be that multiple nodes appear without a classification in your chosen level; the reason for this is, that they have differnt values in some other, higher rank. Such cases are then labeled as follows:<br>",
                                 "<i>the first letter of the taxonomic rank and a number, indicating the number of the node without a rank value</i>[<i>the value of the next higher rank which has a label</i>]"))


networkNodeText = HTML(paste("The <i>nodes</i> in the network will be colored by the determined clusters (using the chosen clustering method).<br>",
                             "The <i>edges</i> are colored by the implied association/correlation: green for positive, red for negative."))


diffNetworkSource = HTML(paste0("<b>NetCoMi</b>: ",
                                "Stefanie Peschel, Christian L Müller, Erika von Mutius, Anne-Laure Boulesteix, Martin Depner <b>2020</b><br>",
                                "<a href=https://doi.org/10.1093/bib/bbaa290> NetCoMi: network construction and comparison for microbiome data in R</a>"))


diffTaxNetworkParameterText = HTML(paste0("Parameters:<br><ul>",
                                          "<li>network methods: (all methods use the abundance values of your experiment, <b>normalized</b> by the currently selected method!)<ol>",
                                          "<li><u>spring:</u> neighborhood selection methodolfy outlines in <a href=https://projecteuclid.org/download/pdfview_1/euclid.aos/1152540754>Meinshausen and Buhlmann (2006)</a> method</li>",
                                          "<li><u>pearson:</u> calculates correlation between nodes using the Pearson correlation</li>",
                                          "<li><u>spearman:</u> calculates correlation between nodes using the Spearman correlation</li>",
                                          "<li><u>spieceasi:</u> use <a href=\'https://doi.org/10.1371/journal.pcbi.1004226\'>  SPIEC-EASIs</a> inverse covariance selection method (glasso) to calculate node connections</li>",
                                          "<li><u>bicor:</u> calculates biweight midcorrelation, a median based measure to measure similarity between groups. Can be an alternative to pearson correlation, since it is more robust to outliers</li>",
                                          "<li><u>SparCC:</u> use the <a href= https://doi.org/10.1371/journal.pcbi.1002687>SparCC</a> algorithm to infer correlation of nodes</li>",
                                          "<li><u>euclidian:</u> calculate the euclidian distance between abundance of nodes</li>",
                                          "<li><u>bray:</u> calculate the Bray-Curtis Dissimilarity between nodes</li>",
                                          "<li><u>jsd:</u> calculate the Jensen–Shannon divergence, which involves approximating the data with a continuous probability density function</li>",
                                          "</ol></li>",
                                          "<li>Cluster methods:<ol>",
                                          "<li><u>cluster_fast_greedy:</u> this method tries to find dense subgraphs (communities) via directly optimizing a modularity score. See details in the paper of <a href=10.1103/PhysRevE.70.066111> Clauset et al (2004)</a>.</li>",
                                          "<li><u>hierarchical:</u> use hierarchical clustering based on dissimilarity values</li>",
                                          "</ol></li>",
                                          "<li>Differential associations methods (only important for differential network):<ol>",
                                          "<li><u>permute:</u> permutational tests accoriding to Gill et al, 2010<ol>",
                                          "<li><u>discordant:</u> discordiant algorithm based on Gaussian mixture model (Siska and Kechris, 2016)<ol>",
                                          "<li><u>fisherTest:</u> Fisher's z-test (Fisher, 1992)<ol>",
                                          "<li>Normalization methods: (<b>Important: only use normalization method on un-normalized projects! Otherwise you will normalize your data twice!</b>)<ol>",
                                          "<li><u>mclr:</u> Modified version of the clr transformation without the need for pseudocounts, by calculating the mean of a sample only with non-zero values. See details of approach in the <a href=https://doi.org/10.3389/fgene.2019.00516>paper of Yoon et al (2019)</a>.</li>",
                                          "<li><u>clr:</u> Centered log-ration transformation. Divide each value in a sample by the mean abundance of that sample; to each value a pseudocount of 1 is applied before applying the standard logarithm to the data. The pseudocount is needed to not get negative values.</li>",
                                          "<li><u>rarefy:</u> rarefy the abundance values to an equal senquencing depth</li>",
                                          "<li><u>TSS:</u> Apply total sum scaling: This method removes technical bias related to different sequencing depth, by dividing each feature count by the total library size.</li>",
                                          "</ol></li>",
                                          "<li>zero replacement methonds:<ol>",
                                          "<li><u>pseudocount:</u> To each entry, where the abundance would be 0, a pseudocount of 1 is added. Therefore no more 0s are present in the data to construct the network</li>",
                                          "<li><u>multiplicative replacement:</u> Similar to pseudocounts, but after adding the pseuodcounts, the 'surrounding' values are adjusted to keep the ratios between them, i.e. not modify their relative relationships</li>",
                                          "</ol></li>",
                                          "</ul>"))


#####confounding texts#####
confoundingInfoText = HTML(paste0("Confounders can be defined as variables whose presence affects the variables being studied so that the results do not reflect the actual relationship. This tab helps you find such confounding factors and to be aware of them for the following analysis steps.</p>",
                                  "<u>When to use it?</u> Always, this is a must-do step when analysing an experiment with multiple variables.",
                                  "<u>Additional Information:</u><br>",
                                  "Simply choose a variable of interest and check the result table. <b>[needs phylogenetic tree file to work",fontawesome::fa("tree", fill="red", height="1.5em"),"]</b><br>
                                  The table tells you which of the other variables is considered a confounding factor and if that result is significant (p-value < 0.05). </p>
                                  In the lower half you find the <b>explained variation</b> of each meta variable, meaning which meta variables account for most of the measured variance."))

randomForestText = HTML(paste0("Here you can apply a machine learning method called “random forest” on your dataset. Its goal is to predict the sub-group (or class) of samples, only using the OTU-abundances and/or other group values (such as medication). It also displays which OTUs are most important to build such a model.</p>",
                               "<u>When to use it?</u> To find OTUs, which are decisive for specific sample-groups."))

###### others #####


cutadaptSourceText = HTML(paste0("<b>cutadapt</b>:",
                                 "Marcel Martin,",
                                 "<a href=https://doi.org/10.14806/ej.17.1.200> Cutadapt removes adapter sequences from high-throughput sequencing reads </a>"))

dada2_filter_info = HTML(paste0("<b>Remove spikes</b>: additional script is used to remove spike-in sequences from your input-files; needs meta-file with a <i>total_weight_in_g</i> column, indicating the weight of spike-ins per sample. <br>",
                                "<b>Trim primers</b>: Removes primer sequences from all reads; <i>V3/V4</i> removes the following primers: <code>CCTACGGGNGGCWGCAG</code> (fw-files) and <code>GACTACHVGGGTATCTAATCC</code> (rv-files) <br>",
                                "<b>Truncation</b>: Removes all bases after the specified base-position for fw or rv files<br>",
                                "<b>ASV removal</b>: Remove all ASVs which have abundance values below x% in <i>all</i> samples<br>",
                                "<b>Phylogenetic tree</b>: additional step to build the phylogenetic tree for your ASVs<br>",
                                "<b>Additional filters applied by default</b>:<br>",
                                "- Truncate reads at the first instance of a quality score less than or equal to 2.<br>",
                                "- sequences with more than 0 Ns will be discarded.<br>",
                                "- reads with higher than 2 <i>expected errors</i> will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))<br>",
                                "- discard reads that match against the phiX genome"))

dada2SourceText = HTML(paste0("<b>dada2</b>: ",
                              "Benjamin J Callahan, Paul J McMurdie, Michael J Rosen, Andrew W Han, Amy Jo A Johnson & Susan P Holmes, <b> 2016 </b>,
                              <a href=https://doi.org/10.1038/nmeth.3869> DADA2: High-resolution sample inference from Illumina amplicon data </a>"))



fastqQualityText = HTML(paste0('This plot shows the base quality at each position for one fastq-file. <br>'))


#### combined texts ####


sourcesText = HTML(paste0("This tool was built using source-code from: <br> ",
                          rheaSourceText ,"<br>",
                          themetagenomicsSourceText,"<br>",
                          "<br>",
                          "The following packages were used for big parts of the calculations and/or plots: <br>",
                          phyloseqSourceText, "<br>",
                          picrust2SourceText,"<br>",
                          dada2SourceText, "<br>",
                          aldexSourceText, "<br>",
                          associationsSourceText, "<br>",
                          diffNetworkSource, "<br>",
                          "<br> A full list of used packages will be provieded here soon..."))


