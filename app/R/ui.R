namco_packages <- c(
  "DT", "networkD3", "shiny", "shinyjs", "waiter", "plotly", "shinyBS",
  "fontawesome", "shinyWidgets", "shinydashboard", "shinydashboardPlus"
)

suppressMessages(lapply(namco_packages, require, character.only = T, quietly = T, warn.conflicts = F))
namco_version <- 'v1.1'

source("texts.R")
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = paste0('Microbiome Explorer (', namco_version,')'), 
    titleWidth = 300,
    dropdownMenuOutput("normalizationDropdown")
  ),
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      br(),
      h2("DATA UPLOAD", style = "text-align:center; font-weight:1000"),
      menuItem("Upload OTU/ASV table", tabName = "uploadOTU", icon = icon("file-upload"), badgeLabel = "START", badgeColor = "green"),
      menuItem("Upload raw fastq files", tabName = "uploadFastq", icon = icon("file-upload"), badgeLabel = "START", badgeColor = "green"),
      menuItem("Use MSD data", tabName="uploadMSD", icon=icon("link"), badgeLabel = "START", badgeColor = "green"),
      fluidRow(
        column(12, align = "center", actionBttn("upload_testdata", "Load sample dataset", icon = icon("database"), style = "gradient", color="warning"))
      ),
      hr(),
      h4("Save or restore session:", style = "text-align:center; font-weight:500"),
      fluidRow(
        column(6, align = "center", downloadBttn("saveSession", "Save Session", size = "xs", style = "float")),
        column(6, align = "center", actionBttn("loadSession", "Restore Session", size = "xs", style = "float"))
      ),
      hr(),
      tags$style(HTML("thead {
                    color: #3c8dbc;
                    }")),
      dataTableOutput("datasets"),
      hr(),
      selectInput("normalizationSelect", "Select normalization strategy", c("no Normalization", "by minimum Sampling Depth", "by Rarefaction", "centered log-ratio", "Total Sum Normalization (normalize to 10,000 reads)")),
      actionBttn("normalizationApply", "Apply normalization strategy", style = "pill", color = "primary", size = "sm", block = F),
      hr(),
      menuItem("Welcome!", tabName = "welcome", icon = icon("door-open"), selected = T),
      menuItemOutput("filtering_menu"),
      menuItemOutput("fastq_overview"),
      menuItemOutput("basic_menu"),
      menuItemOutput("differential_menu"),
      menuItemOutput("functional_menu"),
      menuItemOutput("phylo_menu"),
      menuItemOutput("network_menu"),
      menuItemOutput("ml_menu"),
      menuItemOutput("confounding_menu"),
      menuItem("Info & Settings", tabName = "info", icon = icon("info-circle")),
      hr(),
      fluidRow(
        column(12, h4("Choose global color palette", style = "text-align:center; font-weight:500") )
      ),
      selectInput("namco_pallete","Select global color palette for plots:", choices=c("Set1","Set2","Set3","Paired","Dark2","Accent","Spectral"), selected = "Paired")
    ),
    width = 300, minified=F
  ),
  dashboardBody(
    setShadow(class = "dropdown-menu"),
    use_waiter(),
    waiter_show_on_load(html = tagList(spin_rotating_plane(), "Loading necessary packages for NAMCO ...")),
    useShinyjs(),
    tabItems(
      ##### OTU#####
      tabItem(
        tabName = "uploadOTU",
        h2("Upload your OTU/ASV table"),
        HTML("<h5>[For detailed information on how the files have to look, check out the <b>Info & Settings</b> tab on the left!]</h5>"),
        hr(),
        fluidRow(wellPanel(
          fluidRow(
            column(6, wellPanel(fileInput("otuFile", "Select OTU table"), style = "background:#3c8dbc")),
            column(6, wellPanel(fileInput("metaFile", "Select Metadata File"),
                                style = "background:#3c8dbc",
                                textInput("metaSampleColumnOTU", "Name of the sample-column:", value = "SampleID")
            ))
          ),
          fluidRow(
            column(6, wellPanel(
              checkboxInput("taxInOTU", "Click here if the taxonomic classification is stored in a seperate file and not in the OTU-file:", F),
              fileInput("taxFile", "Select Taxonomic classification file")
            )),
            column(6, wellPanel(fileInput("treeFile", "Select Phylogenetic Tree File (optional)", width = "100%"), fontawesome::fa("tree", fill = "red")))
          ),
          hr(),
          fluidRow(
            column(6, textInput("dataName", "Enter a project name:", placeholder = paste0("Namco_project_", Sys.Date()), value = paste0("Namco_project_", Sys.Date()))),
            column(6, actionBttn("upload_otu_ok", "Upload!", size = "lg", color = "success"))
          )
        ))
      ),
      ##### MSD #####
      tabItem(
        tabName="uploadMSD",
        h2("Upload your MSD data directly into Namco"),
        hr(),
        fluidRow(
          wellPanel(
            p("You can either upload a file with mulitple links from MSD and then choose one or paste a single link directly into the corresponding field."),
            fluidRow(
              column(6, 
                     wellPanel(fileInput("msdFile","Upload MSD file with multiple links"),
                               hidden(selectizeInput("msdLinkSelect","Select an experiment from your file", choices=c())),
                               style = "background:#3c8dbc"),
                     wellPanel(textInput("msdLink", "Enter a single MSD link"), style = "background:#3c8dbc")),
              column(6, wellPanel(
                radioGroupButtons("msdOTUType","You can choose if you want to analyse S-OTUs or zOTUs. Both are contained in the MSD data by default.", choices = c("zOTUs", "S-OTUs")),
                p("If you choose to work with zOTUs, the phylogenetic tree has to be built first. This leads to longer upload times.")
              ))
            )
          )
        ),
        fluidRow(
          column(1),
          column(8, 
                 actionBttn("msdStart", "Upload!", size = "lg", color = "success"),
                 column(6, textInput("msdDataName", "Enter a project name:", placeholder = paste0("Namco_project_", Sys.Date()), value = paste0("Namco_project_", Sys.Date()))),
          )
        )
      ),
      ##### FASTQ#####
      tabItem(
        tabName = "uploadFastq",
        h2("Upload fastq sequencing files"),
        HTML("<h5>[For detailed information on file specifications, check out the <b>Info & Settings</b> tab on the left!]</h5>"),
        hr(),
        fluidRow(wellPanel(
          fluidRow(
            column(6, wellPanel(fileInput("fastqFiles", "Select fastq-files or compressed folder", multiple = T, accept = c(".fastq", ".fastq.gz", ".tar", ".tar.gz", ".zip")),
                                style = "background:#3c8dbc",
                                switchInput("fastqIsPaired", "Select type of experiment", onLabel = "paired-end", offLabel = "single-end", value = T, size = "small")
            )),
            column(6, wellPanel(
              fileInput("fastqMetaFile", "Select Metadata File [optional]"),
              textInput("metaSampleColumnFastq", "Name of the sample-column:", value = "SampleID"),
              textInput("sampleNameCutoff", "Sample-name cutoff", value="_L001"),
              p("For details on the sample-name cutoff feature, check out the sample-names explanation in the \'Info & Settings\' tab!"),
              checkboxInput("fastqApplyRelAbundanceFilter","Apply the recommended 0.25% abundance filter", value = T)
            ))
          )
        )),
        fluidRow(
          column(5),
          column(2, 
                 actionBttn("loadFastqc", "Generate read quality profiles", size = "md", color = "warning"),
                 bsTooltip(id = "loadFastqc",placement = 'top', title = "It is highly advised to first check the sequencing quality of your reads in order to set the parameters below correctly.")
        )),
        hr(),
        fluidRow(
          column(1),
          column(10, box(
            title = 'Explore quality profiles of provided fastq files',
            id = 'fastqcBox',
            width = 12, collapsible = T, collapsed = T, status = 'warning',
            selectInput("qualityUploadSelectSample", "Select Sample", choices = c("Press orange button first to generate quality profiles ..."), width = '50%'),
            fluidRow(
              tabBox(
                title = "Sequencing quality of uploaded fastq-files",
                id = "qualityUploadTabBox", width = 10,
                tabPanel(
                  "Forward",
                  fluidRow(column(12, plotOutput("fastq_file_quality_fw_pre")))
                ),
                tabPanel(
                  "Reverse",
                  fluidRow(column(12, plotOutput("fastq_file_quality_rv_pre")))
                )
              )
            )
          )
        ),
        hr(),
        fluidRow(box(
          title='Select one of the following two amplicon sequencing analysis pipelines to process the fastq files:',
          width=12, solidHeader = T, status = 'primary', background = 'gray',
          fluidRow(
            column(6, box(
              title = 'DADA2', width=12, solidHeader = T, status = 'info', 
              dropdownMenu = boxDropdown(
                icon = icon("info-circle"),
                boxDropdownItem('Publication', icon=icon('asterisk'), href='https://doi.org/10.1038/nmeth.3869'),
                boxDropdownItem('Manual', icon=icon('book'), href='https://benjjneb.github.io/dada2/'),
                dropdownDivider(),
                boxDropdownItem('Additional Parameters', icon=icon('table-list'), href='https://docs.google.com/document/d/1A_3oUV7xa7DRmPzZ-J-IIkk5m1b5bPxo59iF9BgBH7I/edit?usp=sharing')
              ),
              div(style = "display: inline-block;vertical-align:top; width: 150px;", numericInput('trim_primers_fw_dada', 'Trim forward primer (by position)', value = 17)),
              bsTooltip(id = "trim_primers_fw_dada",placement = 'top', title = "Remove first x base positions from each forward read; this is where typcally the primers are located and you can insert the length of your forward primers to remove them"),
              div(style = "display: inline-block;vertical-align:top; width: 150px;", numericInput('trim_primers_rv_dada', 'Trim reverse primer (by position)', value = 21)),
              bsTooltip(id = "trim_primers_rv_dada",placement = 'top', title = "Remove first x base positions from each reverse read; this is where typcally the primers are located and you can insert the length of your reverse primers to remove them"),
              p('Primer trimming in DADA2 is position based, so you need to know the length of your used primers. It is also assumed that primer sequences are at the beginning (left) of each read. If primers are already removed from your reads, enter 0.'),
              div(style = "display: inline-block;vertical-align:top; width: 150px;", numericInput("truncFw", "Truncation foreward:", value = 280, min = 1, max = 500, step = 1)),
              bsTooltip(id = "truncFw",placement = 'top', title = "removes all bases after the specified base-position for foreward files; this is used to remove bases with low sequence quality"),
              div(style = "display: inline-block;vertical-align:top; width: 150px;", numericInput("truncRv", "Truncation reverse:", value = 200, min = 1, max = 500, step = 1)),
              bsTooltip(id = "truncRv",placement = 'top', title = "removes all bases after the specified base-position for reverse files; this is used to remove bases with low sequence quality"),
              p('The quality profiles above can guide you to find more fitting cutoff values'),
              radioGroupButtons("buildPhyloTree", "build phylogenetic tree", c("Yes", "No"), direction = "horizontal", selected = "Yes"),
              bsTooltip(id = "buildPhyloTree",placement = 'top', title = "additional step to build the phylogenetic tree for your ASVs [will increase runtime]")
            )),
            column(6, box(
              title = 'LotuS2', width=12, solidHeader = T, status = 'info',
              dropdownMenu = boxDropdown(
                icon = icon("info-circle"),
                boxDropdownItem('Publication (preprint)', icon=icon('asterisk'), href='https://doi.org/10.1101/2021.12.24.474111'),
                boxDropdownItem('Manual', icon=icon('book'), href='http://lotus2.earlham.ac.uk/main.php?site=documentation'),
                dropdownDivider(),
                boxDropdownItem('Additional Parameters', icon=icon('table-list'), href='https://docs.google.com/document/d/1A_3oUV7xa7DRmPzZ-J-IIkk5m1b5bPxo59iF9BgBH7I/edit?usp=sharing')
              ),
              div(style = "display: inline-block;vertical-align:top; width: 250px;", textInput('trim_primers_fw_lotus', 'Trim forward primer (by sequence)', value = 'CCTACGGGNGGCWGCAG')),
              div(style = "display: inline-block;vertical-align:top; width: 250px;", textInput('trim_primers_rv_lotus', 'Trim reverse primer (by sequence)', value = 'GACTACHVGGGTATCTAATCC')),
              bsTooltip(id = 'trim_primers_fw_lotus', placement = 'top', title='Enter sequence of primer to be removed from forward reads'),
              bsTooltip(id = "trim_primers_rv_lotus",placement = 'top', title = "Enter sequence of primer to be removed from forward reads (will be ignored if single-end experiment is selected above)"),
              p('Primer trimming in LotuS2 is sequence based, so you need to know the sequence of your used primers. The tool will search for this sequence in each read and remove it.'),
              selectInput('clustering_lotus','Select sequence clustering algorithm', choices = c('usearch','cdhit', 'swarm', 'uniose', 'dada2')),
              bsTooltip(id = "clustering_lotus",placement = 'top', title = "LotuS2 offers different clustering algorithms from which you can choose. The default is USEARCH"),
              hidden(htmlOutput('dada2_lotus2_warning')),
              textInput('additional_params_lotus', 'Manually enter additional parameters to pipeline',placeholder = c('-id 0.97 -buildPhylo 0')),
              bsTooltip(id = "additional_params_lotus",placement = 'top', title = "LotuS2 has many more parameters that you can check out in their manual (see info box in the right corner). Simply add them to this text box as shown in the example."),
              p('Lotus2 will build a phylogenetic tree by default. If you do not need it, simply insert -buildPhylo 0 in the text field above.')
            ))
          ),
          fluidRow(
            column(2),
            column(2, actionBttn("upload_fastq_dada2", "Start DADA2", size = "lg", color = "success")),
            column(1),
            column(2, textInput("fastqDataName", "Enter a project name:", placeholder = paste0("Namco_project_", Sys.Date()), value = paste0("Namco_project_", Sys.Date()))),
            column(1),
            column(3, actionBttn("upload_fastq_lotus2", "Start LotuS2", size = "lg", color = "success"))
          )
              
        ))
      )),
      ##### MSD
      tabItem(
        tabName="uploadMSD",
        h2("Upload data from MSD"),
        p("Explain MSD + link"),
        hr(),
        fluidRow(wellPanel(
          fluidRow(
            column(6, wellPanel(
              fileInput("uploadMSDTable","Upload the table with MSD links"),
              textInput("uploadMSDLink", "Enter a single MSD link"),
              style = "background:#3c8dbc")),
            column(6, wellPanel(
              checkboxGroupButtons("MSDUseSOTUs","Use S-OTUs or ASVs", choices = c("S-OTUs", "ASVs"))
            ))
          )  
        ))
      ),
      ##### welcome#####
      tabItem(
        tabName = "welcome",
        fluidRow(
          column(12, wellPanel(fluidRow(
            column(2, htmlOutput("logo")),
            column(6, htmlOutput("welcome")),
            column(4, htmlOutput("biomedLogo"))
          ), style="border:4px solid #3c8dbc; margin-bottom: 1px"))
        ),
        fluidRow(
          
          column(3),
          column(6, div(HTML("<center><h3>Start by uploading your data or use our provided sample dataset and try out all the features in <i>NAMCO</i></h3></center>")))
        ),
        hr(),
        fluidRow(
          column(6, fluidRow(box(title="Documentation", htmlOutput("documentation"), solidHeader=T, status="primary",collapsible = T, collapsed = F, width = 12))),
          column(3, 
                 fluidRow(box(title="Issues & Recommendations", htmlOutput("contactText"), solidHeader=T, status="primary",collapsible = T, collapsed = F, width = 12)),
                 fluidRow(box(title="News", htmlOutput("newsText"), solidHeader=T, status="primary", collapsible=T, collapsed=F, width=12))),
          column(3,
                 fluidRow(box(title = "Authors", htmlOutput("authors"), solidHeader = T, status = "primary", collapsible = T, collapsed = T, width = 12)),
                 fluidRow(box(title = "References", htmlOutput("welcome_ref"), solidHeader = T, status = "primary", collapsible = T, collapsed = T, width = 12)))
        ),
        fluidRow(
          column(2),
          column(8, fluidRow(box(title="Overview of NAMCO", htmlOutput("workflow"), solidHeader=T, status="primary",collapsible = T, collapsed = F, width=12)))
        )
      ),
      ##### overview+filter#####
      tabItem(
        tabName = "overview",
        h4("Data Overview & Filtering"),
        fluidRow(
          valueBoxOutput("otus_box1"),
          valueBoxOutput("samples_box1"),
          valueBoxOutput("conditions_box1")
        ),
        h5("Please consider filtering your data!"),
        fluidRow(
          valueBoxOutput("below_abundance_filter_features1"),
          valueBoxOutput("below_prevalence_filter_features1")
        ),
        tabBox(
          id = "filters", width = 12,
          
          tabPanel(
            "Data Overview",
            hr(),
            fluidRow(
              column(12, wellPanel(
                fluidRow(
                  column(3, selectInput("downloadMetaOTUTaxLevel","Combine OTUs/ASVs to this taxonomic level (for download only)", choices=c("OTU/ASV","Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))),
                  column(2, radioButtons('downloadMetaOTUrelativeAbundance',label = 'Download relative or absolute abundance values', choices = c('absolute', 'relative'))),
                  column(7, downloadBttn("downloadMetaOTU", "Download combined (abundance + meta) data", color = "royal",size = "lg"))
                )
              ))
            ),
            h3("Explore the meta-file you uploaded:"),
            fluidRow(
              column(12, wellPanel(
                dataTableOutput("metaTable")
              ))
            )
          ),
          
          tabPanel(
            "Filter Samples & taxonomic levels",
            p("Here you can filter samples and taxonomic levels"),
            hr(),
            fluidRow(
              column(4, wellPanel(
                h4("Filter options for samples"),
                p("Select which sample you want to remove; you can select samples by name or remove whole sample groups."),
                selectInput("filterColumns", "Sample groups", choices = ""),
                selectInput("filterColumnValues", "Group values", choices = "",multiple = T),
                pickerInput("filterSample", "Pick specific samples by name", choices = "", multiple = T, options = list(`actions-box` = T)),
                hr(),
                actionButton("filterApplySamples", "Apply Filter", style = "background-color:blue; color:white; display:inline-block"),
                actionButton("filterResetA", "Restore original dataset", style = "background-color:green; color:white")
              )),
              column(4, wellPanel(
                h4("Filter options for taxa"),
                p("Select which taxonomic features you want to keep (you can select multiple)."),
                selectInput("filterTaxa", "Taxonomic level", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")),
                pickerInput("filterTaxaValues", "Feature name", choices = "", multiple = T, options = list(`actions-box` = T)),
                hr(),
                actionButton("filterApplyTaxa", "Apply Filter", style = "background-color:blue; color:white"),
                actionButton("filterResetB", "Restore original dataset", style = "background-color:green; color:white")
              )),
              column(4, wellPanel(
                h4("Filter History"),
                htmlOutput("filterHistoryText")
              ))
            )
          ),
          tabPanel(
            "Advanced Filtering",
            p("Apply more advanced filterings onto your dataset. Click the checkbox to the left of a filtering function and select a fitting value; multiple combinations of functions can be applied. Click 'Apply Filter' when you are done."),
            p("The plots to the left show you the current distribution of the corresponding value as well as a red line to indicate the currently selected filtering value."),
            p("Note: these functions are applied to the normalized dataset!"),
            hr(),
            fluidRow(
              column(
                3, prettyCheckbox("advFilterMinAbundance", "Filter by minimum abundance", F, status = "primary", shape = "curve", animation = "smooth"),
                p("Remove all OTUs with a summed up abundance value over all samples below ... (left of red line gets removed)")
              ),
              column(3, disabled(numericInput("advFilterMinAbundanceValue", "Value:", 1, 1, 20000, 1))),
              column(6, plotlyOutput("advFilterMinAbundancePlot", height = "200px"))
            ),
            hr(),
            fluidRow(
              column(
                3, prettyCheckbox("advFilterRelAbundance", "Filter by relative abundance", F, status = "primary", shape = "curve", animation = "smooth"),
                p("Remove all OTUs with a relative abundance value in all samples below ... %"),
                p("The red line shows the chosen cutoff in relation to all appearing rel. abundance values.")
              ),
              column(3, disabled(numericInput("advFilterRelAbundanceValue", "Value:", 0.25, 0.001, 100, 0.01))),
              column(6, plotlyOutput("advFilterRelAbundancePlot", height = "200px"))
            ),
            hr(),
            fluidRow(
              column(
                3, prettyCheckbox("advFilterNumSamples", "Filter by occurrence in x samples", F, status = "primary", shape = "curve", animation = "smooth"),
                p("Remove all OTUs which occur at least ... times in all samples (meaning, they have at most ... times an abundance value of 0) (right of red line gets removed)")
              ),
              column(3, disabled(numericInput("advFilterNumSamplesValue", "Value:", 20, 1, 500, 1))),
              column(6, plotlyOutput("advFilterNumSamplesPlot", height = "200px"))
            ),
            hr(),
            fluidRow(
              column(
                3, prettyCheckbox("advFilterMaxVariance", "Filter by highest variance", F, status = "primary", shape = "curve", animation = "smooth"),
                p("Keep only those ... OTUs with the highest variance in abundance (left of red line gets removed)")
              ),
              column(3, disabled(numericInput("advFilterMaxVarianceValue", "Value:", 20, 1, 500, 1))),
              column(6, plotlyOutput("advFilterMaxVariancePlot", height = "200px"))
            ),
            hr(),
            fluidRow(
              column(
                3, prettyCheckbox("advFilterPrevalence", "Filter by prevalence", F, status = "primary", shape = "curve", animation = "smooth"),
                p("Keep only those OTUs whith a prevalence value over ... %")
              ),
              column(3, disabled(numericInput("advFilterPrevalenceValue", "Value:", 10, 0.1, 100, 0.1))),
              column(6, plotlyOutput("advFilterPrevalencePlot", height = "200px"))
            ),
            hr(),
            fluidRow(column(1), column(1, actionButton("filterApplyAdv", "Apply Filter", style = "background-color:blue; color:white")), column(1, actionButton("filterResetC", "Restore original dataset", style = "background-color:green; color:white"))),
          ),
          
          tabPanel(
            'Decontamination',
            p('Remove contaminated OTUs/ASVs using DNA concentration and/or control samples.'),
            hr(),
            fluidRow(
              column(8, box(
                title = span( icon("info"), "Tab-Information"),
                htmlOutput("decontamText"),
                solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
              ))
            ),
            h5("Step 1: Select corresponding data columns"),
            fluidRow(
              column(5, wellPanel(
                p('(a) DNA concentration'),
                selectInput('DNAconcentrationColumn','Select column with DNA concentration info for each sample:', choices = c()),
              )),
              column(5, wellPanel(
                p('(b) Prevalence'),
                selectInput('controlSamplesColumn', 'Select column that indicates control samples from true samples:', choices = c()),
                selectizeInput('controlSamplesName', 'Select, which of the value present in the selected column is the indicator for control samples.', choices = c()),
                box(
                  title='Library size per sample',
                  plotOutput('librarySizePlot'),
                  solidHeader = T, status='info', width=12, collapsible = T, collapsed = T
                ),
                p("The above plot shows the library size of each sample. Control samples are usually located at the lower end of the library size distribution in a dataset."),
                downloadLink("librarySizePlotPDF", "Download as PDF")
              )),
              column(2, wellPanel(
                numericInput('decontamThreshold', 'Probability threshold to reject H0', min = 0, max=1, value = .1, step = 0.01),
                actionBttn("startDecontam", "Find contaminants!", size = "md", color = "warning"),
                p('Scroll down once finished to inspect the detected contaminants...')),
                downloadButton("contamTableDownload", "Download results as table")
              )
            ),
            hr(),
            fixedRow(
              column(5, h5("Step 2: Inspect possible contaminants")),
              column(5),
              column(2, h5('Step 3: select and remove contaminants'))
            ),
            fluidRow(
              column(5, wellPanel(
                selectInput('contamCandidatesSelect', 'Select contaminate candidate to inspect', choices =c()),
                box(
                  title='Feature frequency vs DNA concentration',
                  plotOutput('contamDiagnosticDNA'),
                  solidHeader = T, status='info', width=12, collapsible = T, collapsed = T
                ),
                p('The above plot shows if the selected feature (OTU/ASV) is a robust candidate for a contaminant. 
                  The red line shows the model of a contaminant sequence feature, for which frequency is expected to be inversely proportional to input DNA concentration.
                  This means if the black points follow this red line, this feature is in line with the model and can be considered a robust result.'),
                downloadLink("contamDiagnosticDNA_PDF", "Download as PDF")
              )),
              column(5, wellPanel(
                box(
                  title='Prevalence of features in control vs true samples',
                  plotlyOutput('contamDiagnosticPrev'),
                  solidHeader = T, status='info', width=12, collapsible = T, collapsed = T
                ),
                p('The above plot shows how often the individual OTUs/ASVs were observed in the control and true samples. Contaminants should clearly form a branch that distinguishes them from features that are not contaminated. You can hover over the individual points to find out the name of potential outliers.'),
                downloadLink("contamDiagnosticPrev_PDF", "Download as PDF")
              )),
              column(2, wellPanel(
                pickerInput('contamCandidatesSelectRemove', 'Select all contaminants to remove', choices=c(), multiple = T, options = list(`actions-box` = TRUE)),
                actionBttn("removeContam", "Remove contaminants!", size = "md", color = "warning")
              ))
            )
          ),
          
          tabPanel(
            "Add Meta-Data",
            p("Upload a meta-file, which assigns groups and values to your samples."),
            hr(),
            fluidRow(
              column(6, wellPanel(
                fileInput("metaFileAdditional", "Select Metadata file"),
                textInput("metaAdditionalSampleColumn", "Name of the sample-column:", value = "SampleID"),
                style = "background:#3c8dbc")),
              column(6, actionBttn("upload_meta_ok", "Upload!", size = "lg", color = "success"))
            )
          )
        )
      ),
      ##### fastq overview#####
      tabItem(
        tabName = "fastq_tab",
        h4("fastq Overview"),
        fluidRow(
          valueBoxOutput("otus_box2"),
          valueBoxOutput("samples_box2"),
          valueBoxOutput("conditions_box2")
        ),
        fluidRow(
          tabBox(
            id = "fastq_overview", width = 12,
            tabPanel(
              "Quality and Filtering",
              h3("Analysis of sequence quality for provided fastq files after pipeline"),
              fluidRow(column(
                12,
                selectInput("fastq_file_select_post", label = "Select fastq-pair:", multiple = F, choices = c()),
                fluidRow(
                  column(6, wellPanel(
                    h4("foreward"),
                    div("", plotOutput("fastq_file_quality_fw_post"))
                  )),
                  column(6, wellPanel(
                    h4("reverse"),
                    div("", plotOutput("fastq_file_quality_rv_post"))
                  ))
                )
              )),
              fluidRow(column(10, htmlOutput("fastqQualityText"))),
              hr(),
              fluidRow(
                column(2),
                column(8, box(
                  id='dada2_readloss_box',
                  title='Number of reads after each step in DADA2 pipeline', width=12,
                  solidHeader = T, collapsible = T, collapsed = T, status = 'info',
                  plotlyOutput("fastq_pipeline_readloss")
                ))
              )
            ),
            tabPanel(
              'Pipeline logs',
              h3('Access logging infos from chosen pipeline'),
              p('[Currently only available for LotuS2 pipeline]'),
              hr(),
              fluidRow(
                column(6, box(
                  title='Demultiplexing log', status = 'black', width=12, solidHeader = T,
                  uiOutput('demulti_log_output')
                )),
                column(6, box(
                  title='Complete run log', status = 'black', width=12, solidHeader = T,
                  uiOutput('run_log_output')
                ))
              )
            ),
            tabPanel(
              "Downloads",
              h3("Download options"),
              hr(),
              fluidRow(
                column(
                  12,
                  h3("Download the generated ASV-tables:"), wellPanel(
                    fixedRow(
                      column(4, downloadBttn("download_asv_norm", "Download normalized ASV table", style = "float", size = "sm")),
                      column(4, downloadBttn("download_asv_raw", "Download unnormalized ASV table", style = "float", size = "sm"))
                    )
                  ),
                  h3("Download the ASV sequences:"), wellPanel(
                    fixedRow(column(6, downloadBttn("download_asv_fastq", "Download fasta file of ASV sequences", style = "float", size = "sm")))
                  ),
                  h3("Download the taxonomic classification:"), wellPanel(
                    fixedRow(column(6, downloadBttn("download_taxonomy", "Download taxonomic classification of ASVs", style = "float", size = "sm")))
                  ),
                  h3("Download a phyloseq R-object:"), wellPanel(
                    fixedRow(column(6, downloadBttn("download_phyloseq", "Download phyloseq object", style = "float", size = "sm")))
                  )
                )
              )
            )
          )
        )
      ),
      ##### basic#####
      tabItem(
        tabName = "basics",
        h4("Basic Analysis"),
        fluidRow(
          tabBox(
            id = "basicPlots", width = 12,
            tabPanel(
              "Taxonomic Binning",
              h3("Analyse samples by their taxonomic composition"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("taxBinningText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fixedRow(
                column(9,
                       #shinycssloaders::withSpinner(plotlyOutput("taxaDistribution", height = "auto")),
                       plotlyOutput("taxaDistribution", height = "auto"),
                       downloadLink("taxaPDF", "Download as PDF")
                ),
                column(3, box(
                  width = 12,
                  title = "Options",
                  solidHeader = T, status = "primary",
                  selectInput("taxBinningLevel", "Select taxonomic level to display", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")),
                  selectInput("taxBinningGroup", "Split by sample group", choices = c("None")),
                  selectInput("taxBinningYLabel", "Select label for y-axis", choices = c("None")),
                  switchInput("taxBinningShowY","Show labels on y axis",value=T, onLabel="Show",offLabel="Hide",size="mini", inline = T),
                  switchInput("taxaAbundanceType", "Show relative or absolute abundance", onLabel = "relative", offLabel = "absolute", value = T, size = "mini", inline = T),
                  numericInput("taxBinningTop", "Show top K taxa", value = 10, min = 1, step = 1),
                  checkboxInput('taxBinningRotate',label='Rotate plot', value=F),
                  box(title='Change order of y axis', 
                      checkboxInput('taxBinningOrderManually', label = 'Order manually', value = F),
                      selectizeInput("taxBinningYOrder", 'Manually change order of y-axis by removing values from this list and placing them at a different position (using your cursor)', choices = c(), multiple = T),
                      selectInput("taxBinningOrderReference", "Select taxon by which to order the bars", choices=c("None")),
                      solidHeader = F, status = 'info', width=12, collapsible = T, collapsed = T
                  )
                ))
              )
            ),
            tabPanel(
              "Alpha Diversity",
              h3("Analyse the diversity of species inside of samples"),
              tags$hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("alphaDivText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              fluidRow(
                column(9, wellPanel(
                  plotOutput("alphaPlot", height = "600px"),
                  downloadLink("alphaPDF", "Download as PDF")
                )),
                column(3, wellPanel(
                  selectInput("alphaMethod", "Method:", c("Shannon_Entropy", "effective_Shannon_Entropy", "Simpson_Index", "effective_Simpson_Index", "Richness"), multiple = T, selected = "Richness"),
                  selectInput("alphaGroup", "Group by:", c("-")),
                  radioGroupButtons("alphaScalesFree", "Free y-scale", choices = c("free", "fixed"), selected = "free", direction = "horizontal"),
                  pickerInput("alphaPairs", "Select which sub-group pairs you want to compare with the wilcoxon test", choices = c(), multiple = T, options = list(`actions-box` = T)),
                  radioGroupButtons("alphaSignifView", "Display p-value or significance codes", choices = c("p-value", "codes"), selected = "p-value", direction = "horizontal"),
                  #selectInput("alphaPalette", "Select Color-palette", choices = c("Rainbow (use this if you have many groups)" = "rainbow", "JCO" = "jco", "NPG" = "npg", "AAAS" = "aaas", "NEJM" = "nejm", "Lancet" = "lancet", "JAMA" = "jama", "UCSCGB" = "ucscgb", "Star Trek" = "startrek"), selected = "jco")
                ))
              ),
              br(), br(),
              fluidRow(column(
                12,
                h4("Raw values for alpha diversity scores, including download:"),
                downloadButton("alphaTableDownload", "Download Table")
              )),
              fluidRow(
                column(10, wellPanel(
                  tableOutput("alphaTable")
                ))
              ),
              tags$hr(),
              fluidRow(
                column(1),
                column(10, htmlOutput("alphaDivFormulas"))
              )
            ),
            tabPanel(
              "Beta Diversity",
              h3("Analyse the diversity of species between samples using non-euclidean distance measures"),
              tags$hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("betaDivText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(
                column(9, wellPanel(
                  h5("Hierarchical clustering (Ward's method) of the sample using the chosen distance method"),
                  plotOutput("betaTree", width = "100%", height="600px"),
                  downloadLink("betaTreePDF", "Download as PDF")
                )),
                column(3, wellPanel(
                  selectInput("betaMethod", "Method to calculate distances between samples:", choices = ""),
                  selectInput("betaGroup", "Color samples by the following group:", choices = ""),
                  selectInput("betaGroup2","Add second grouping (by shape) to plot:", choices = ""),
                  selectizeInput("betaLevel", "Display beta-diversitsy of selected group level:", choices = ""),
                  checkboxInput('betaShowLabels','Show sample labels in PCoA and NMDS plots', value = T),
                  downloadLink('betaDownloadDistance', 'Download distance matrix')
                ))
              ),
              hr(),
              h5("Dimensionality reduction based on ecological distances"),
              fluidRow(
                column(12,
                       tabsetPanel(type="tabs",
                                   tabPanel("Non-metric multidimensional scaling (NMDS)",
                                            fluidRow(
                                              column(6, plotOutput("betaDivNMDS", height='600px'), downloadLink('betaDivNMDSPDF', 'Download as PDF')),
                                              column(6, plotOutput("betaDivStress", height='600px'), downloadLink('betaDivStressPDF', 'Download as PDF'))
                                            )),
                                   tabPanel("Principal Coordinate Analysis (PCoA)", 
                                            fluidRow(
                                              column(6, plotOutput("betaDivPcoa", height='600px'), downloadLink('betaDivPocaPDF', 'Download as PDF'))
                                            ))
                       ))
              ),
              hr(),
              fluidRow(
                column(3, p('Change position of text in scatterplot'),)
              ),
              fluidRow(
                column(3, sliderInput('betaDivTextSliderX', 'Move in direction of x-axis', min = -.5, max=.5, step = 0.01, value = 0.1)),
                column(3, sliderInput('betaDivTextSliderY', 'Move in direction of y-axis', min = -.5, max=.5, step = 0.01, value = 0.1))
              )
            ),
            tabPanel(
              "Rarefaction Curve",
              h3("Analysis of species richness with rarefaction curves"),
              tags$hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("rarefactionInfoText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              fluidRow(
                column(8, wellPanel(plotlyOutput("rarefacCurve", width = "100%"), waiter_hide_on_render(id = "rarefacCurve"))),
                column(4, wellPanel(
                  sliderInput("rareToHighlight", "Number of samples with steepest rarefaction slope to be highlighted:", min = 0, value = 1, max = 100, step = 1),
                  br(),
                  verbatimTextOutput("undersampled"),
                  switchInput("excludeSamples", "exclude undersampled samples", value = F)
                )),
                valueBoxOutput("samples_box3")
              )
            ),
            tabPanel(
              "Abundance Heatmaps",
              h3("Generate ecologically-organized heatmaps"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("heatmapText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(
                column(10, wellPanel(
                  plotlyOutput("abundanceHeatmap", height = '500px'),
                  downloadLink("abundanceHeatmapPDF", "Download as PDF")
                )),
                column(
                  2,
                  box(
                    title = "Options",
                    selectInput("heatmapDistance", "Choose distance method", choices = c()),
                    selectInput("heatmapOrdination", "Choose Orientation Method (Ordination)", choices = c("NMDS", "MDS/PCoA", "DPCoA", "DCA", "CCA", "RDA")),
                    selectInput("heatmapSample", "Choose labeling of X-axis", choices = ""),
                    selectInput('heatmapOverlayTaxa', 'Select taxonomic level that is displayed in hover overlay', choices=c("OTU/ASV","Kingdom", "Phylum", "Class", "Order", "Family", "Genus")),
                    checkboxInput("heatmapOrderSamples", "Order samples by selected sample group",value = F),
                    # selectInput("heatmapTrans","Color scale transformation (numerical transformation of observed value)", choices=c("NONE"=NULL, "log(4)[default]"="log_trans(4)", "log(2)"="log2_trans()","log()+1"="log1p()")),
                    solidHeader = T, status = "primary", width = 12
                  )
                )
              ),
              hr(),
              fluidRow(
                column(
                  10,
                  br(),
                  htmlOutput("heatmapOrdinationText"),
                  br(),
                  htmlOutput("heatmapSourceText"),
                  htmlOutput("neatmapSourceText")
                )
              )
            )
          )
        )
      ),
      ##### differential#####
      tabItem(
        tabName = "differential",
        h4("Differential Analysis"),
        fluidRow(
          tabBox(
            id = "differentialPlots", width = 12,
            tabPanel(
              "Associations",
              h3("Explore different measures of association between sample groups"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("associationsText"),
                  br(),
                  htmlOutput("associationsSourceText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(
                column(9, plotOutput("associationsPlot", width = "100%")),
                column(
                  3, selectInput("associations_level", "Choose level of association testing", choices = c("OTU/ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
                  box(
                    width = 12,
                    title = "Options",
                    solidHeader = T, status = "primary",
                    sliderInput("associations_alpha", "Significance level", 0.00001, 0.999, 0.05, 0.001),
                    selectInput("associations_label", "Select meta-label, for which associations are tested", c("")),
                    selectizeInput("associations_case", "Select, which value is considered case (will be compared against all other values in label)", c("")),
                    sliderInput("assiciation_show_numer", "How many significant features do you want to display?", 1, 100, 25, 1),
                    selectInput("associations_sort", "Select how to sort the displayed features", choices = c("p-value", "fold-change", "prevalence shift")),
                    selectInput("associations_panels", "Which additional values do you want to display?", choices = c("fold-change", "AU-ROC", "prevalence"), multiple = T),
                    hr(),
                    actionBttn("associations_start", "Generate Plot...", icon = icon("play"), style = "pill", color = "primary", block = T, size = "md")
                  ),
                  downloadLink("associationsPDF", "Download as PDF"),
                  downloadLink("associationsTable", "Download significant features as table")
                )
              )
            ),
            tabPanel(
              "Correlations",
              h3("Find correlations between OTUs/ASVs and phenotypes"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("corrText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(
                column(9, plotOutput("corrPlot", width = "100%")),
                column(3, box(
                  width = 12,
                  title = "Options",
                  solidHeader = T, status = "primary",
                  selectInput("corrTaxLevel","Select taxonomic level", choices=c("OTU/ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), selected="Phylum"),
                  switchInput("corrIncludeTax", "Include OTUs/ASVs", onLabel = "Yes", offLabel = "No", value = F),
                  pickerInput("corrSelectGroups", "Select meta-variables which you want to include", choices = "", multiple = T, options = list(`actions-box` = T)),
                  numericInput("corrSignifCutoff", "Select significance cutoff", value = 0.05, min = 0.0001, max = 1, step = 0.01),
                  numericInput("corrCorrelationCutoff", "Select absolute correlation cutoff", value = 0.7, min = 0.01, max = 1, step = 0.01),
                  radioGroupButtons("corrPval", "How to display non-significant correlations", choices = c("highlight", "blank", "do nothing"), selected = "do nothing", direction = "horizontal"),
                  sliderInput("corrTextSize","Change size of axis labels", value = 0.4, min=0.1, max=1, step = 0.1)
                ), downloadLink("corrPlotPDF", "Download as PDF"))
              )
            ),
            
            tabPanel(
              "Topic Modeling",
              h3("Explore OTU topics in your dataset"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("topicText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              fluidRow(
                column(
                  3,
                  sliderInput("K", "Pick Number of Topics:", 2, 50, 8, step = 1),
                  htmlOutput("topic_text")
                ),
                # column(3,
                #  sliderInput("sigma_prior","Pick Scalar between 0 and 1:", 0, 1, 0, step=0.01),
                #  htmlOutput("sigma_text")),
                column(
                  3,
                  selectInput("formula", label = "Formula for covariates of interest found in metadata:", choices = "Please provide OTU-table & metadata first!"),
                  selectInput("refs", label = "Binary covariate in formula, indicating the reference level:", choices = "Please provide OTU-table & metadata first!")
                ),
                column(
                  3,
                  # actionButton("themeta","Visualize topics!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                  actionBttn("themeta", "Visualize topics!", icon = icon("play"), style = "pill", size = "md", color = "primary"),
                  br(),
                  downloadButton("downloadThemeta","Download Topic assignment")
                  # ,downloadButton("downloadGeneTable","Download Gene-Table")
                )
              ),
              tags$hr(),
              fluidRow(
                column(3, h4("Input Variables:")),
                column(10, htmlOutput("input_variables")),
                column(1)
              ),
              tags$hr(),
              fixedRow(
                column(1),
                column(10, htmlOutput("text1")),
                column(1)
              ),
              br(),
              fixedRow(
                column(10, plotlyOutput("est", height = "200px"), p("To download as PNG, klick the photo-icon in the plot. PDF downloas is not supported yet."),),
                p("Topics colored red, have a strong association with the chosen reference level;
                  the blue topics on the other hand are associated with the other level within the chosen covariate.
                  (Example: Chosen covariate is Gender and reference level is Female -> Female will be colored red, Male is blue)")
              ),
              br(),
              fixedRow(
                column(1, radioButtons("dim", label = strong("Dim"), choices = list("2D" = "2d", "3D" = "3d"), selected = "2d")),
                column(3, selectInput("dist", label = strong("Method"), choices = list(
                  "Bray Curtis" = "bray",
                  "Jaccard" = "jaccard", "Euclidean" = "euclidean", "Hellinger" = "hellinger", "Chi Squared" = "chi2",
                  "Jensen Shannon" = "jsd"
                ), selected = "jsd")),
                column(1, style = "padding: 25px 0px;", actionButton("reset", "Reset")),
                column(2, numericInput("k_in", label = strong("Topic Number"), value = 0, min = 0, max = 100, step = 1)),
                column(3, sliderInput("lambda", label = strong("Lambda"), min = 0, max = 1, value = 1, step = 0.01)),
                column(2, selectInput("taxon", label = strong("Taxon"), choices = list(
                  "Phylum" = "Phylum",
                  "Class" = "Class", "Order" = "Order", "Family" = "Family", "Genus" = "Genus"
                )))
              ),
              fixedRow(
                column(1, tags$div("Number of components to plot.", class = "capt")),
                column(3, tags$div("Type of distance and method for ordination.", class = "capt")),
                column(1, tags$div("Reset topic selection.", class = "capt")),
                column(2, tags$div("Current selected topic.", class = "capt")),
                column(3, tags$div(paste0(
                  "Relative weighting of selected topic that influences taxa shown in barplot.",
                  " If equal to 1, p(taxa|topic)l if 0, p(taxa|topic)/p(taxa)."
                ), class = "capt")),
                column(2, tags$div("Taxonomic group to dictate bar plot shading", class = "capt"))
              ),
              fixedRow(
                column(6, offset = 0, height = "600px", plotlyOutput("ord")),
                column(6, offset = 0, height = "600px", plotOutput("bar"))
              ),
              fixedRow(
                column(6, p("To download as PNG, klick the photo-icon in the plot. PDF download is not supported yet."), 
                       tags$div(paste0(
                         "Ordination of the samples over topics distribution theta, colored according to",
                         " the weights shown in the scatter plot above. The radius of a given point",
                         " represents the marginal topic frequency. The amount of variation explained is",
                         " annotated on each axis."
                       ), class = "below")),
                column(6, downloadLink("themetaBarPDF", "Download as PDF"), 
                       tags$div(paste0(
                         "Bar plot representing the taxa frequencies. When no topic is selected, the overall",
                         " taxa frequencies are shown, colored based on the selected taxonomy and ordered in",
                         " in terms of saliency. When a topic is chosen, the red bars show the margina taxa",
                         " frequency within the selected topic, ordered in terms of relevency, which in turn",
                         " can be reweighted by adjusting the lambda slider."
                       ), class = "below"))
              ),
              br(),
              fixedRow(
                column(1),
                column(10, htmlOutput("text2")),
                column(1)
              ),
              forceNetworkOutput("corr")
            ),
            
            tabPanel(
              "Time-series analysis",
              h3("Compare time-points and find trends in your dataset"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("timeSeriesText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              p("In order to start your analysis, select the correct meta groups in the options menu & press the 'Perform Analysis' button"),
              fluidRow(
                column(9, tabsetPanel(type="tabs",
                                      tabPanel("Line-Plot", plotOutput("timeSeriesPlot",height = "800px")),
                                      tabPanel("Significant features",
                                               p("Here you can see the significance value for each taxa and meta variable regarding the selected time-points. Since the values are -log10 transformed, the ones at the top behave significantly different over the time-points."),
                                               plotOutput("timeSeriesSignifFeatures", height="800px"), downloadLink("timeSeriesSignifTable","Download as table"))
                )),
                column(3, box(
                  width = 12,
                  title = "Options",
                  solidHeader = T, status = "primary",
                  h4("Fixed Options"),
                  p("If you change a variable here, you need to hit the \'Perform Analysis\' button again."),
                  hidden(selectInput("timeSeriesTaxa","Select taxonomic level you want to analyse", choices=c("OTU/ASV","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))),
                  selectInput("timeSeriesMeasure", "Select which abundance measure you want to compare over the time-points", choices=c("Abundance", "relative Abundance", "Richness","Shannon_Entropy", "effective_Shannon_Entropy", "Simpson_Index", "effective_Simpson_Index")),
                  hidden(numericInput("timeSeriesClusterK", "Select the number of clusters you want to produce (scroll down for help on finding a meaningful value); 0 means no clustering; Note: this will override the selected mean group", min=0, max=100, value=0, step=1)),
                  actionBttn("timeSeriesStart", "Perform analysis", style = "pill", size = "lg", color = "primary"),
                  hr(),
                  h4("Interactive Options"),
                  selectInput("timeSeriesGroup","Select group which represents time-points or something comparable (x-axis)", choices = c()),
                  selectizeInput("timeSeriesTimePointOrder", "Manually change order of time-points (you can delete and add them at the current cursor position)", choices=c(), multiple = T),
                  selectInput("timeSeriesBackground", "Select group which represent the groups over time-points (e.g. patients)", choices=c()),
                  hidden(pickerInput("timeSeriesTaxaSelect", "Select which taxa to show", choices=c(), multiple=T, options=list(`liveSearch` = T,`actions-box` = T))),
                  hidden(selectInput("timeSeriesMeanLine","Select group, for which to display the mean line", choices=c())),
                  sliderInput("timeSeriesLineSize", "Change mean-line size", min = 0, max=5, step=0.1, value=1),
                  selectizeInput("timeSeriesSampleHighlight","Highlight a specific line", choices=c()),
                  selectInput("timeSeriesHighlightColor", "Choose color for highlight", choices=c("red","green","blue","black","orange","purple")),
                  radioGroupButtons("timeSeriesAdjPval","Display default or BH-adjusted p-values", choices=c("default","adjusted"), direction = "horizontal")
                ),
                downloadLink("timeSeriesPlotPDF", "Download as PDF"))
              ),
              hr(),
              h4("If you chose to cluster your samples, you can find additional statistics & information down here:"),
              p("(You first have to select how many clusters you want in the 'Options' menu above)"),
              fluidRow(
                column(6, wellPanel(plotOutput("timeSeriesOptimalClustersPlot"))),
                column(6, wellPanel(plotOutput("timeSeriesClusterSizePlot"), p("(The coloring is controlled by the groups over time-points selection above)")))
              ),
              fluidRow(
                column(6, wellPanel(h4("See, which sample was selected into which cluster:"),
                                    dataTableOutput("timeSeriesClusterContent")))
              )
            ),
            
            tabPanel(
              "Differential statistical Analysis",
              h3("Perform statistical tests on different taxonomic levels between sample groups"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("statTestText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(
                column(12, wellPanel(
                  fluidRow(
                    column(4, 
                           selectInput("statTestMethod", "Select which test you want to perform", choices=c("Wilcoxon test", "Kruskal-Wallis test")),
                           numericInput("statTestCutoff", "Select significance cutoff", value = 0.05, min = 0.01, max = 1, step = 0.01)
                    ),
                    column(4,
                           selectInput("statTestcompLevel", "Select taxonomic level on which to perform test", choices=c("Phylum", "Class", "Order", "Family", "Genus","OTU/ASV"), selected="OTU/ASV"),
                           selectInput("statTestGroup", "Select sample group you want to analyse", choices=c()),
                    ),
                    column(4, actionBttn("statTestStart", "Perform test",icon = icon("play"), style = "pill", size = "md", color = "primary"),
                           p("Select a multiple testing correction method:"),
                           radioGroupButtons("statTestPAdjust",choices = c("Bonferroni"="bonferroni","Benjamini-Hochberg"="BH","None"="none"),direction = "horizontal",individual = T))
                  )
                ))
              ),
              hr(),
              fluidRow(
                column(3, valueBoxOutput("statTestSignifCount")),
                column(3, downloadBttn("statTestDownloadTable", "Download table of significant taxa (p-values)", size = "sm")),
                column(3, downloadBttn("statTestAbundanceDownloadTable", "Download table of significant taxa (abundance)", size = "sm"))
              ),
              hr(),
              fluidRow(
                column(9, plotOutput("statTestPlot")),
                column(3, box(
                  width = 12,
                  title = "Options",
                  solidHeader = T, status = "primary",
                  selectInput("statTestSignifPicker","The significant taxa are listed here; pick one to display the boxplot", choices = c()),
                  hidden(pickerInput("statTestPairPicker", "Select which sub-group pairs you want to display", choices = c(), multiple = T, options = list(`actions-box` = T, `live-search` = TRUE))),
                  sliderInput("statTestLabelSize","Change label size", min = 0, max=100, value=10, step=1)
                ))
              ),
              fluidRow(column(3,downloadLink("statTestPDF", "Download current plot as PDF"))),
              fluidRow(column(3,downloadLink("statTestPDFall","Download all significant features as PDF")))
            )
          )
        )
      ),
      ##### functional#####
      tabItem(
        tabName = "functional",
        h4("Functional Analysis"),
        fluidRow(
          tabBox(
            id = "functionalPlots", width = 12,
            tabPanel(
              "Functional prediction",
              h3("Functional prediction using Picrust2"),
              tags$hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("picrust2Text"),
                  tags$br(),
                  htmlOutput("picrust2SourceText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              fluidRow(
                column(6, wellPanel(
                  h3("Picrust2 parameters & Input Files"),
                  p("Note: Picrust2 will always be applied on the non-normalized dataset automatically."),
                  fileInput("picrustFastaFile", "Upload .fasta file with sequences for your OTUs/ASVs:", accept = c()),
                  checkboxInput("picrust_copy_number_normalization", "Normalize OTU abundances by copy-number", value = T),
                  p("Next to the functional assignment of OTUs, Picrust2 also infers the copy numbers of each 16s-rRNA gene per OTU; you have the option to normalize your abundance values with the copy-numbers by selecting this checkbox.")
                ),
                infoBoxOutput("hasPicrustInfoBox")),
                column(6, 
                       actionBttn("picrust2Start", "Start picrust2 analysis!", icon = icon("play"), style = "pill", color = "primary", block = T, size = "lg"),
                       wellPanel(
                         p("Download zip-archive with raw picrust2 results:"),
                         h4("Please be aware:"),
                         p("This will create a zip archive of all output files, so it might take a few seconds until the download window appears!"),
                         p("This download window will not appear if you use a restored dataset!"),
                         hidden(div(
                           id = "download_picrust_div",
                           downloadButton("download_picrust_raw", "Download picrust2 results as zip archive:")
                         ))
                       ))
              ),
              hr(),
              fluidRow(
                column(6, wellPanel(
                  h3("Differential analysis parameters"),
                  radioGroupButtons("picrustTest", "Select which statistical test you want to perform for the differential analysis", choices = c("Welch's t-test"="t", "Wilcoxon Rank Sum test"="wilcox", "Kruskal Wallis test"="kw"), direction = "horizontal", individual = T),
                  radioGroupButtons("picrustTestNormalization","Select method of normalization for picrust2 results", choices=c("relative abundance"="rel", "centered log-ratio"="clr","none"="none")),
                  selectInput("picrust_test_condition", "Choose condition for which to test differential abundance", choices = c()),
                  selectizeInput("picrust_test_covariate", "Choose covariate against which to compare all other samples", choices = c()),
                  hidden(div(id="aldex2Additional",
                             numericInput("picrust_mc_samples", "Choose number of MC iterations", min = 4, max = 1000, value = 128, step = 4),
                             p("A higher number of MC iterations will increase precision of estimating the sampling error but also increase runtime. For datasets with few samples a higher value can be chosen, with more samples a lower one should be used."),
                  ))
                )),
                column(6, 
                       actionBttn("picrustDiffStart","Start differential analysis",icon = icon("play"), style = "pill", color = "primary", block = T, size = "lg"),
                       wellPanel(
                         p("This button is only activated if you have run the picrust analysis!"),
                         p("You can change the parameters on the left and rerun the analysis & reload the plots by clicking this button."),
                         downloadButton("picrustDiffDownloadEC","Download results of analysis (EC)"),
                         downloadButton("picrustDiffDownloadKO","Download results of analysis (KO)"),
                         downloadButton("picrustDiffDownloadPW","Download results of analysis (PW)")
                       )
                )),
              hr(),
              h3("Differential functional analysis"),
              htmlOutput("aldexSourceText"),
              hr(),
              fluidRow(column(
                12,
                tabBox(
                  title = "Relationships between effect size & p-value",
                  id = "picrust_tabBox", width = 12,
                  tabPanel(
                    "EC",
                    fluidRow(
                      column(6, div("", plotOutput("picrust_ec_effect_plot")), downloadLink("picrust_ec_effectPDF", "Download as PDF")),
                      column(6, div("", plotOutput("picrust_ec_vulcano_plot")), downloadLink("picrust_ec_vulcanoPDF", "Download as PDF"))
                    ),
                    fluidRow(
                      column(
                        8,
                        box(
                          title = "Names of significant functions",
                          verbatimTextOutput("picrust_ec_effect_signif"),
                          status = "info", width = 12, collapsible = T, collapsed = T
                        )
                      ),
                      column(4, valueBoxOutput("picrust_ec_effect_signif_value"))
                    ),
                    h3("Details about significant functions:"),
                    hr(),
                    fluidRow(
                      column(10, plotOutput("picrust_ec_signif_plot"), downloadLink("picrust_ec_signifPDF", "Download as PDF")),
                      column(2, 
                             numericInput("picrust_ec_signif_plot_show", "Maximum number of displayed ECs", 20, min = 1, max = 100, step = 1),
                             pickerInput("picrust_ec_select", "Select specific EC to display", choices=c(), multiple = T, options = list(`liveSearch` = T), width = "fit"),
                             checkboxInput("picrust_show_descripton_ec","Show description of function",value = F),
                             sliderInput("picrust_ylab_size_ec","Change label size on y axis", min=1, max=100, value=10, step=1)
                      )
                    ),
                    p("Here the functions with BH adjusted P-value above the significance threshold are displayed; the boxplot shows the different abundance distributions of a function colored by each sample group. Also the BH-adjusted P-value and effect size is displayed as a barplot.")
                  ),
                  tabPanel(
                    "KO",
                    fluidRow(
                      column(6, div("", plotOutput("picrust_ko_effect_plot")), downloadLink("picrust_ko_effectPDF", "Download as PDF")),
                      column(6, div("", plotOutput("picrust_ko_vulcano_plot")), downloadLink("picrust_ko_vulcanoPDF", "Download as PDF"))
                    ),
                    fluidRow(
                      column(
                        8,
                        box(
                          title = "Names of significant functions",
                          verbatimTextOutput("picrust_ko_effect_signif"),
                          status = "info", width = 12, collapsible = T, collapsed = T
                        )
                      ),
                      column(4, valueBoxOutput("picrust_ko_effect_signif_value"))
                    ),
                    h3("Details about significant functions:"),
                    fluidRow(
                      column(10, plotOutput("picrust_ko_signif_plot"), downloadLink("picrust_ko_signifPDF", "Download as PDF")),
                      column(2, 
                             numericInput("picrust_ko_signif_plot_show", "Maximum number of displayed KOs", 20, min = 1, max = 100, step = 1),
                             pickerInput("picrust_ko_select", "Select specific KO to display", choices=c(), multiple = T, options = list(`liveSearch` = T)),
                             checkboxInput("picrust_show_descripton_ko","Show description of function",value = F),
                             sliderInput("picrust_ylab_size_ko","Change label size on y axis", min=1, max=100, value=10, step=1)
                      )
                    ),
                    p("Here the functions with BH adjusted P-value above the significance threshold are displayed; the boxplot shows the different abundance distributions of a function colored by each sample group. Also the BH-adjusted P-value and effect size is displayed as a barplot.")
                  ),
                  tabPanel(
                    "PW",
                    fluidRow(
                      column(6, div("", plotOutput("picrust_pw_effect_plot")), downloadLink("picrust_pw_effectPDF", "Download as PDF")),
                      column(6, div("", plotOutput("picrust_pw_vulcano_plot")), downloadLink("picrust_pw_vulcanoPDF", "Download as PDF"))
                    ),
                    fluidRow(
                      column(
                        8,
                        box(
                          title = "Names of significant functions",
                          verbatimTextOutput("picrust_pw_effect_signif"),
                          status = "info", width = 12, collapsible = T, collapsed = T
                        )
                      ),
                      column(4, valueBoxOutput("picrust_pw_effect_signif_value"))
                    ),
                    h3("Details about significant functions:"),
                    fluidRow(
                      column(10, plotOutput("picrust_pw_signif_plot"), downloadLink("picrust_pw_signifPDF", "Download as PDF")),
                      column(2, 
                             numericInput("picrust_pw_signif_plot_show", "Set max. number of displayed PWs", 20, min = 1, max = 100, step = 1),
                             pickerInput("picrust_pw_select", "Select specific PW to display", choices=c(), multiple = T, options = list(`liveSearch` = T)),
                             checkboxInput("picrust_show_descripton_pw","Show description of function",value = F),
                             sliderInput("picrust_ylab_size_pw","Change label size on y axis", min=1, max=100, value=10, step=1)
                      )
                    ),
                    p("Here the functions with BH adjusted P-value above the significance threshold are displayed; the boxplot shows the different abundance distributions of a function colored by each sample group. Also the BH-adjusted P-value and effect size is displayed as a barplot.")
                  ),
                  tabPanel(
                    "Information & Options",
                    fluidRow(
                      column(6, wellPanel(
                        h3("Options for Visualization"),
                        numericInput("picrust_signif_lvl", "Change significance level (P-value)", min = 0.01, max = 1, value = 0.05, step = 0.01),
                        p("Here you can set the significance cutoff for the BH adjusted P-value; functions with a p-value below it are considered significant."),
                        numericInput("picrust_signif_lvl_effect", "Change significance level (effect size)", min = -10, max = 10, value = 1, step = 0.01),
                        p("Here you can set the significance cutoff for the effect size; functions with a effect size greater than it it are considered significant. "),
                        sliderInput("picrust_maxoverlaps", "Change number of overlaps for point labels", min = 1, max = 500, value = 50, step = 1),
                        checkboxInput("picrust_signif_label", "Label significant functions in scatterplots", value = F),
                        p("If too many points in close proximity are considered significant, change the number of overlaps, to display more labels.")
                      )),
                      column(6, wellPanel(
                        h3("Information"),
                        htmlOutput("picrust_pval_info_text")
                      ))
                    )
                  ),
                  tabPanel(
                    "Downloads",
                    fluidRow(
                      column(
                        5,
                        wellPanel(downloadBttn("picrust_download_ec", "Download differential analysis of EC numbers (EC)", style = "float", size = "sm")),
                        wellPanel(downloadBttn("picrust_download_ko", "Download differential analysis of KEGG ortholog groups (KO)", style = "float", size = "sm")),
                        wellPanel(downloadBttn("picrust_download_pw", "Download differential analysis of pathways (PW)", style = "float", size = "sm"))
                      )
                    )
                  )
                )
              ))
            )
          )
        )
      ),
      ##### phylogenetic#####
      tabItem(
        tabName = "phylogenetic",
        h4("Phylogenetic Analysis"),
        fluidRow(
          tabBox(
            id = "phylogeneticPlots", width = 12,
            tabPanel(
              "Phylogenetic Tree",
              h3("Phylogenetic Tree of OTU taxa", fontawesome::fa("tree", fill = "red")),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("phyloTreeText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fixedRow(
                column(6, wellPanel(
                  h4("Basic tree visualization options:"),
                  div(
                    id = "phylo_basic",
                    sliderInput("phylo_prune", "Number of OTUs to display (pick the x OTUs with the highest cumulative abundance):", 0, 1, 1, 1),
                    selectInput("phylo_tiplabels", "Label tips (remove OTU labels by choosing \'-\'):", choices = c("taxa_names", "-")),
                    selectInput("phylo_method", "Visualization Method:", choices = c("circular" = "circular", "rectangular" = "rectangular")),
                    selectInput("phylo_edge_length", "Align tips (works better in rectangular mode)", choices = c("Yes", "No")),
                    selectInput("phylo_draw_clado", "Draw Cladogram (dont use calculated branch lengths)", choices = c("Yes" = "none", "No" = "branch.length")),
                    selectInput("phylo_group", "Pick a categorical group from your meta file to add a heatmap", choices = c("NONE")),
                    selectInput("phylo_taxonomy", "Pick a taxonomic group to add a heatmap", choices = c("NONE", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
                  )
                )),
                column(6, wellPanel(
                  h4("Advanced tree visualization options:"),
                  actionButton("phylo_toggle_advanced", "Show/hide advanced options"),
                  hidden(div(
                    id = "phylo_advanced",
                    p("Additional options to manage tree:"),
                    numericInput("phylo_size_tree", "Choose size of tree labels", value = 3, min = 0, step = 1, max = 50),
                    numericInput("phylo_offset", "Choose offset of heatmap to tree", value = 1, min = 0, step = 1, max = 50),
                    numericInput("phylo_width_taxonomy", "Choose width of heatmap boxes (taxonomy heatmap)", value = 1, min = 0.1, max = 1, step = 0.1),
                    numericInput("phylo_width_meta", "Choose width of heatmap boxes (group heatmap); only categorical variables are shown", value = 1, min = 0.1, max = 1, step = 0.1),
                    hr()
                  ))
                ))
              ),
              hr(),
              fixedRow(
                column(12, wellPanel(
                  div(
                    id = "phylo_tree",
                    plotOutput("phyloTree"), style = "height:1200px"
                  )
                ))
              )
            )
          )
        )
      ),
      ##### network#####
      tabItem(
        tabName = "network",
        h4("Network Analysis"),
        fluidRow(
          tabBox(
            id = "netWorkPlots", width = 12,
            tabPanel(
              "Co-occurrence of OTUs/ASVs",
              h3("Co-occurrences network generation"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("basic_info"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(
                column(1),
                column(5, numericInputIcon("binCutoff", "Cutoff for Binarization (OTUs with a smaller value are considered as not present)", min = 0.0000001, max = 100, value = 1, step = 0.01, icon = icon("cut"), width = "400px")),
                column(6, box(
                  title = "Information",
                  htmlOutput("basic_additional"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              h3("You can look at the following two plots to see the effect of your chosen cutoff:"),
              fixedRow(
                column(
                  6,
                  box(
                    title = "Cutoff-Barplot",
                    plotlyOutput("cutoffHist"),
                    htmlOutput("cutoff_text"),
                    solidHeader = T, status = "primary", collapsible = T, collapsed = T, width = 12
                  )
                ),
                column(
                  6,
                  box(
                    title = "Cutoff-Heatmap",
                    plotlyOutput("boolHeat"),
                    htmlOutput("heatmap_text"),
                    solidHeader = T, status = "primary", collapsible = T, collapsed = T, width = 12
                  )
                )
              ),
              tags$hr(),
              htmlOutput("basic_calc_title"),
              fluidRow(
                column(
                  6,
                  wellPanel(
                    radioButtons("useFC", "Calculation of Counts:", c("log2(fold-change)", "difference")),
                    selectInput("groupCol", "Select which sample group is to be compared (minimum of 2 levels in group!):", choices = c("Please Upload OTU & META file first!"), selected = "Please Upload OTU & META file first!"),
                    selectizeInput("groupVar1", "Select variable of group to compare with (reference group)", choices = c("Please Upload OTU & META file first!")),
                    selectInput("groupVar2", "Select variable of group to compare against (choose *all* to compaire against all other variables in group)", choices = c("Please Upload OTU & META file first!"))
                  ),
                  actionBttn("startCalc", " Start count-calculation & (re-)load network!", icon = icon("play"), style = "pill", color = "primary", block = T, size = "lg")
                ),
                column(
                  6,
                  box(
                    title = "Information",
                    htmlOutput("basic_calc_additional"),
                    solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                  )
                )
              ),
              tags$hr(),
              htmlOutput("basic_network_title"),
              fluidRow(
                column(
                  9,
                  forceNetworkOutput("basicNetwork"),
                ),
                column(
                  3,
                  box(
                    title = "Options",
                    sliderInput("networkCutoff", "Number of edges to show (edges are sorted by most extreme values, positive and negative):", 1, 5000, 100, 10),
                    selectInput("netLevel", "Color Taxonomic Level:", choices = c("-", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
                    solidHeader = T, status = "primary", width = 12
                  ),
                  downloadLink("basicNetworkPDF", "Download as HTML"),
                  box(
                    title = "chosen parameters",
                    htmlOutput("chosen_network_params"),
                    solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = F
                  ),
                  p("Green edges: OTU-pair more often occuring in selected reference group."),
                  p("Red edges: OTU-pair more often occuring in other sample group, which is compared against.")
                )
              )
            ),
            
            tabPanel(
              "Network inference",
              h3("Create a single network on your whole dataset"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("compNetworkInfoText"),
                  br(),
                  htmlOutput("networkNodeTextCopyCopy"),
                  br(),
                  htmlOutput("diffNetworkSourceCopyCopy"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              fluidRow(
                column(8, box(
                  title = "Parameter-information",
                  htmlOutput("compNetworkParamsText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(column(
                12,
                wellPanel(
                  fluidRow(
                    column(
                      3,
                      selectInput("compNetworkMeasure", "Choose the measure used for calculation of network", choices = c("spring", "pearson", "spearman", "spieceasi", "bicor", "sparcc")),
                      selectInput("compNetworkColor","How should nodes be colored", choices=c("by detected clusters"="cluster", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
                    ),
                    column(
                      3,
                      selectInput("compNetworkClustMethod", "Choose method how to detect clusters in network", choices = c("cluster_fast_greedy", "hierarchical")),
                      selectInput("compNetworkNormMethod", "Choose normalization method (in order to make counts of different samples comparable)", choices = c("none", "mclr", "clr", "rarefy", "TSS")),
                      selectInput("compNetworkzeroMethod", "Choose method how to replace zeros in data", choices = c("none", "add pseudocount of 1 to data" = "pseudo", "mulitplicative replacement" = "multRepl"))
                    ),
                    box(
                      width = 4,
                      title = "Additional Parameters",
                      solidHeader = T, status = "info", collapsed = T, collapsible = T,
                      hidden(div(
                        id = "compNetworkAdditionalParamsSPRING.EASIDiv",
                        numericInput("compNetworkNlambda", "Number of lambdas", 10, 1, 100, 1),
                        numericInput("compNetworkRepNum", "Number of subsamples for StARS", 20, 1, 100, 1),
                        numericInput("compNetworkLambdaRatio", "Smallest value for lambda", 0.1, 0, 1, 0.01)
                      )),
                      hidden(div(
                        id = "compNetworkAdditionalParamsSPARCCdiv",
                        numericInput("compNetworkIter", "Number of iterations in outer loop", 20, 1, 100, 1),
                        numericInput("compNetworkInnerIter", "Number of iterations in inner loop", 10, 1, 100, 1),
                        numericInput("compNetworkTh", "Threshold for correlations", 0.1, 0, 1, 0.01)
                      ))
                    ),
                    column(
                      2,
                      actionBttn("compNetworkCalculate", "Start Calculation", style = "pill", size = "lg", color = "primary")
                    )
                  )
                )
              )),
              hr(),
              fluidRow(
                column(
                  9,
                  conditionalPanel(
                    condition = "input.compNetworkInteractiveSwitch == true",
                    wellPanel(forceNetworkOutput("compNetworkInteractive", height = "800px"))
                  ),
                  conditionalPanel(
                    condition = "input.compNetworkInteractiveSwitch == false",
                    plotOutput("compNetwork", height = "800px")
                  ),
                  downloadLink("comp_networkPDF", "Download as PDF")
                ),
                box(
                  width = 3,
                  title = "Display options",
                  solidHeader = T, status = "primary",
                  switchInput("compNetworkInteractiveSwitch", label = "Interactive network", value = F, onLabel = "Yes", offLabel = "No"),
                  selectInput("compNetworkLayout", "Layout", choices = c("spring", "circle", "Fruchterman-Reingold" = "layout_with_fr")),
                  hr(),
                  h4("Node options:"),
                  # selectInput("taxNetworkNodeColor","Choose how to color nodes", choices=c("by detected clusters"="cluster", "Kingdom", "Phylum")),
                  selectInput("compNetworkNodeFilterMethod", "Choose method how to filter out nodes (keep top x nodes with ...)", choices = c("none", "highestConnect", "highestDegree", "highestBetween", "highestClose", "highestEigen")),
                  numericInput("compNetworkNodeFilterValue", "Choose x for the node filtering method", value = 100, min = 1, max = 10000, step = 1),
                  #selectInput("compNetworkRmSingles", "How to handle unconnected nodes (all: remove all; inboth: only if unconnected in both networks, none: remove no unconnected nodes)", choices = c("inboth", "none", "all")),
                  selectInput("compNetworkNodeSize", "Choose value which indicates the size of nodes", choices = c("fix", "degree", "betweenness", "closeness", "eigenvector", "counts", "normCounts", "clr", "mclr", "rarefy", "TSS")),
                  hr(),
                  h4("Edge options:"),
                  selectInput("compNetworkEdgeFilterMethod", "Choose method how to filter out edges (threshold: keep edges with weight of at least x; highestWeight: keep first x edges with highest weight)", choices = c("none", "threshold", "highestWeight"), selected = "highestWeight"),
                  numericInput("compNetworkEdgeFilterValue", "Choose x for edge filtering method", value = 300, min = 1, max = 5000, step = 1),
                  hr(),
                  h4("Other options:"),
                  sliderInput("compNetworkTitleSize", "Change title size", min=0, max=10, value=2, step = 0.1)
                )
              ),
              hr(),
              h3("Details about network"),
              fluidRow(column(
                8,
                wellPanel(verbatimTextOutput("compNetworkSummary"))
              ))
            ),
            
            tabPanel(
              "Taxonomic Rank Networks",
              h3("Explore network structures between the discovered taxonomic ranks"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("taxNetworkInfoText"),
                  br(),
                  htmlOutput("networkNodeText"),
                  br(),
                  htmlOutput("diffNetworkSourceCopy"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              fluidRow(
                column(8, box(
                  title = "Parameter-information",
                  htmlOutput("taxNetworkParamsText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(column(
                12,
                wellPanel(
                  fluidRow(
                    column(
                      3,
                      selectInput("taxNetworkRank", "Select taxonomic rank", choices = c("Phylum", "Class", "Order", "Family", "Genus")),
                      selectInput("taxNetworkMeasure", "Choose the measure used for calculation of network", choices = c("spring", "pearson", "spearman", "spieceasi", "bicor", "sparcc"))
                    ),
                    column(
                      3,
                      selectInput("taxNetworkClustMethod", "Choose method how to detect clusters in network", choices = c("cluster_fast_greedy", "hierarchical")),
                      selectInput("taxNetworkNormMethod", "Choose normalization method (in order to make counts of different samples comparable)", choices = c("none", "mclr", "clr", "rarefy", "TSS")),
                      selectInput("taxNetworkzeroMethod", "Choose method how to replace zeros in data", choices = c("none", "add pseudocount of 1 to data" = "pseudo", "mulitplicative replacement" = "multRepl"))
                    ),
                    box(
                      width = 4,
                      title = "Additional Parameters",
                      solidHeader = T, status = "info", collapsed = T, collapsible = T,
                      hidden(div(
                        id = "taxNetworkAdditionalParamsSPRING.EASIDiv",
                        numericInput("taxNetworkNlambda", "Number of lambdas", 10, 1, 100, 1),
                        numericInput("taxNetworkRepNum", "Number of subsamples for StARS", 20, 1, 100, 1),
                        numericInput("taxNetworkLambdaRatio", "Smallest value for lambda", 0.1, 0, 1, 0.01)
                      )),
                      hidden(div(
                        id = "taxNetworkAdditionalParamsSPARCCdiv",
                        numericInput("taxNetworkIter", "Number of iterations in outer loop", 20, 1, 100, 1),
                        numericInput("taxNetworkInnerIter", "Number of iterations in inner loop", 10, 1, 100, 1),
                        numericInput("taxNetworkTh", "Threshold for correlations", 0.1, 0, 1, 0.01)
                      ))
                    ),
                    column(
                      2,
                      actionBttn("taxNetworkCalculate", "Start Calculation", style = "pill", size = "lg", color = "primary")
                    )
                  )
                )
              )),
              hr(),
              fluidRow(
                column(
                  9,
                  conditionalPanel(
                    condition = "input.taxNetworkInteractiveSwitch == true",
                    wellPanel(forceNetworkOutput("taxNetworkInteractive", height = "800px"))
                  ),
                  conditionalPanel(
                    condition = "input.taxNetworkInteractiveSwitch == false",
                    plotOutput("taxNetwork", height = "800px")
                  ),
                  downloadLink("tax_networkPDF", "Download as PDF")
                ),
                box(
                  width = 3,
                  title = "Display options",
                  solidHeader = T, status = "primary",
                  switchInput("taxNetworkInteractiveSwitch", label = "Interactive network", value = F, onLabel = "Yes", offLabel = "No"),
                  selectInput("taxNetworkLayout", "Layout", choices = c("spring", "circle", "Fruchterman-Reingold" = "layout_with_fr")),
                  hr(),
                  h4("Node options:"),
                  # selectInput("taxNetworkNodeColor","Choose how to color nodes", choices=c("by detected clusters"="cluster", "Kingdom", "Phylum")),
                  selectInput("taxNetworkNodeFilterMethod", "Choose method how to filter out nodes (keep top x nodes with ...)", choices = c("none", "highestConnect", "highestDegree", "highestBetween", "highestClose", "highestEigen")),
                  numericInput("taxNetworkNodeFilterValue", "Choose x for the node filtering method", value = 100, min = 1, max = 10000, step = 1),
                  selectInput("taxNetworkRmSingles", "How to handle unconnected nodes (all: remove all; inboth: only if unconnected in both networks, none: remove no unconnected nodes)", choices = c("inboth", "none", "all")),
                  selectInput("taxNetworkNodeSize", "Choose value which indicates the size of nodes", choices = c("fix", "degree", "betweenness", "closeness", "eigenvector", "counts", "normCounts", "clr", "mclr", "rarefy", "TSS")),
                  hr(),
                  h4("Edge options:"),
                  selectInput("taxNetworkEdgeFilterMethod", "Choose method how to filter out edges (threshold: keep edges with weight of at least x; highestWeight: keep first x edges with highest weight)", choices = c("none", "threshold", "highestWeight"), selected = "none"),
                  numericInput("taxNetworkEdgeFilterValue", "Choose x for edge filtering method", value = 300, min = 1, max = 5000, step = 1),
                  hr(),
                  h4("Other options:"),
                  sliderInput("taxNetworkTitleSize", "Change title size", min=0, max=10, value=2, step = 0.1)
                )
              ),
              hr(),
              h3("Details about network"),
              fluidRow(column(
                8,
                wellPanel(verbatimTextOutput("taxNetworkSummary"))
              ))
            ),
            tabPanel(
              "Differential Networks",
              h3("Explore network structures in different sample groups"),
              hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("diffNetworkInfoText"),
                  htmlOutput("networkNodeTextCopy"),
                  br(),
                  htmlOutput("diffNetworkSource"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              fluidRow(
                column(8, box(
                  title = "Parameter-information",
                  htmlOutput("diffNetworkParameterText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              hr(),
              fluidRow(column(
                12,
                wellPanel(
                  fluidRow(
                    column(
                      3,
                      selectInput("diffNetworkSplitVariable", "Choose sample group you want to compare (only groups with 2 levels are shown)", choices = c()),
                      selectInput("diffNetworkMeasure", "Choose the measure used for calculation of network", choices = c("spring", "pearson", "spearman", "spieceasi", "bicor", "sparcc")),
                      selectInput("diffNetworkDiffMethod","Choose method used for calculating differential associations",choices=c("fisherTest","permute","discordant")),
                      hidden(selectInput("diffNetworkColor","How should nodes be colored", choices=c("by detected clusters"="cluster", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")))
                    ),
                    column(
                      3,
                      selectInput("diffNetworkTaxaLevel","Choose taxonomic level",choices=c("OTU/ASV","Kingdom", "Phylum", "Class", "Order", "Family", "Genus")),
                      selectInput("diffNetworkClustMethod", "Choose method how to detect clusters in network", choices = c("cluster_fast_greedy", "hierarchical")),
                      selectInput("diffNetworkNormMethod", "Choose normalization method (in order to make counts of different samples comparable)", choices = c("none", "mclr", "clr", "rarefy", "TSS")),
                      selectInput("diffNetworkzeroMethod", "Choose method how to replace zeros in data", choices = c("none", "add pseudocount of 1 to data" = "pseudo", "mulitplicative replacement" = "multRepl"))
                      # numericInput("diffNetworkSparsMethodParams","A Students t-test is used to select a subset of edges which are connected; choose significance level here",value = 0.05,min = 0.001,max=1,step = 0.001)
                    ),
                    box(
                      width = 4,
                      title = "Additional Parameters",
                      solidHeader = T, status = "info", collapsed = T, collapsible = T,
                      hidden(div(
                        id = "diffNetworkAdditionalParamsSPRING.EASIDiv",
                        numericInput("diffNetworkNlambda", "Number of lambdas", 10, 1, 100, 1),
                        numericInput("diffNetworkRepNum", "Number of subsamples for StARS", 20, 1, 100, 1),
                        numericInput("diffNetworkLambdaRatio", "Smallest value for lambda", 0.1, 0, 1, 0.01)
                      )),
                      hidden(div(
                        id = "diffNetworkAdditionalParamsSPARCCdiv",
                        numericInput("diffNetworkIter", "Number of iterations in outer loop", 20, 1, 100, 1),
                        numericInput("diffNetworkInnerIter", "Number of iterations in inner loop", 10, 1, 100, 1),
                        numericInput("diffNetworkTh", "Threshold for correlations", 0.1, 0, 1, 0.01)
                      ))
                    ),
                    column(
                      2,
                      actionBttn("diffNetworkCalculate", "Start Calculation", style = "pill", size = "lg", color = "primary")
                    )
                  )
                )
              )),
              hr(),
              
              fluidRow(
                column(9,
                       tabsetPanel(id="diffnet_tabs",type="tabs",
                                   tabPanel("2 group network",
                                            conditionalPanel(
                                              condition = "input.diffNetworkInteractiveSwitch == true",
                                              column(6, forceNetworkOutput("groupNetworkInteractive1", height = "800px")),
                                              column(6, forceNetworkOutput("groupNetworkInteractive2", height = "800px"))
                                            ),
                                            conditionalPanel(
                                              condition = "input.diffNetworkInteractiveSwitch == false",
                                              wellPanel(plotOutput("groupNetwork", height = "800px"))
                                            ),
                                            downloadLink("group_networkPDF", "Download as PDF")         
                                   ),
                                   tabPanel("Differential network",
                                            wellPanel(plotOutput("diffNetwork", height="800px")),
                                            downloadLink("diff_networkPDF", "Download as PDF")
                                   )
                       )
                ),
                box(
                  width = 3,
                  title = "Display options",
                  solidHeader = T, status = "primary",
                  hidden(switchInput("diffNetworkInteractiveSwitch", label = "Interactive network", value = F, onLabel = "Yes", offLabel = "No")),
                  selectInput("diffNetworkLayout", "Layout", choices = c("spring", "circle", "Fruchterman-Reingold" = "layout_with_fr")),
                  hr(),
                  h4("Node options:"),
                  hidden(selectInput("diffNetworkNodeFilterMethod", "Choose method how to filter out nodes (keep top x nodes with ...)", choices = c("none", "highestConnect", "highestDegree", "highestBetween", "highestClose", "highestEigen"))),
                  hidden(numericInput("diffNetworkNodeFilterValue", "Choose x for the node filtering method", value = 100, min = 1, max = 10000, step = 1)),
                  hidden(selectInput("diffNetworkRmSingles", "How to handle unconnected nodes (all: remove all; inboth: only if unconnected in both networks, none: remove no unconnected nodes)", choices = c("inboth", "none", "all"))),
                  hidden(selectInput("diffNetworkNodeSize", "Choose value which indicates the size of nodes", choices = c("fix", "degree", "betweenness", "closeness", "eigenvector", "counts", "normCounts", "clr", "mclr", "rarefy", "TSS"))),
                  hr(),
                  h4("Edge options:"),
                  selectInput("diffNetworkEdgeFilterMethod", "Choose method how to filter out edges (threshold: keep edges with weight of at least x; highestWeight: keep first x edges with highest weight)", choices = c("none", "threshold", "highestWeight"), selected = "highestWeight"),
                  numericInput("diffNetworkEdgeFilterValue", "Choose x for edge filtering method", value = 100, min = 1, max = 5000, step = 1),
                  hr(),
                  h4("Other options:"),
                  sliderInput("diffNetworkTitleSize", "Change title size", min=0, max=10, value=2, step = 0.1),
                  sliderInput("diffNetworkLabelSize", "Change node label size", min=0,max=10,value=2,step=0.1)
                )
              ),
              hr(),
              h3("Details about network difference"),
              fluidRow(column(
                8,
                wellPanel(verbatimTextOutput("groupNetworkSummary"))
              ))
            )
          )
        )
      ),
      ##### confounding#####
      tabItem(
        tabName = "confounding",
        h4("Confounding Analysis"),
        fluidRow(
          tabBox(
            id = "confoundingPlots", width = 12,
            tabPanel(
              "Confounding Analysis & Explained Variation",
              h3("Analyse confounding factors", fontawesome::fa("tree", fill = "red")),
              tags$hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("confoundingInfoText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              fluidRow(
                column(1),
                column(2, selectInput("confounding_distance","Select distance method for confounding analysis", choices=c("Bray-Curtis"))),
                column(3, disabled(actionBttn("confounding_start", "Check for confounding factors", icon=icon("play"), style="pill", color="primary",block=T,size="md")))
              ),
              tags$hr(),
              p("This heatmap will show you if there are confounding factors for specific variables. Read the plot from x to y axis like this: tested variable XX has possible confounding factors YY (if the legend says 'yes', YY is a confounder for XX)."),
              fluidRow(
                column(9, plotOutput("confounding_heatmap", height = "800px")),
                column(3,
                       box(
                         width = 12,
                         title = "Options",
                         solidHeader = T, status = "primary",
                         selectInput("confounding_heatmap_type","Select value to show in heatmap", choices = c("is_confounder","direction","pvalue")), 
                         pickerInput("confounding_select_tested_variable", "Select which tested variables to show", choices=c(), multiple = T),
                         sliderInput("confounding_label_size", "Change label text size", min = 0, max=100, value=10, step=0.1),
                         downloadButton("confounding_table_download", "Download Table with results"),
                         downloadLink("confounding_PDF_download","Download as PDF")
                       )
                )
              ),
              tags$hr(),
              h4("Explained Variation:"),
              fluidRow(
                column(1),
                p("The bars represent rsquare (4 digit rounded value is written over bars) and are colored by pvalue. The rsquare value corresponds to the explained variation a variable has"),
                column(7, plotOutput("explainedVariationBar", height = "700px"))
              ),
            )
          )
        )
      ),
      ##### machine learning #####
      tabItem(
        tabName = "machineLearning",
        h4("Machine Learning"),
        fluidRow(
          tabBox(
            id = "machineLearngingPlots", width = 12,
            tabPanel(
              "Random Forests",
              h3("Build a Random Forest machine learning model"),
              tags$hr(),
              fluidRow(
                column(8, box(
                  title = span( icon("info"), "Tab-Information"),
                  htmlOutput("randomForestText"),
                  solidHeader = F, status = "info", width = 12, collapsible = T, collapsed = T
                ))
              ),
              tags$hr(),
              fixedRow(
                column(
                  6, wellPanel(
                    h2("Options for building the model:"),
                    selectInput("forest_variable", "Choose variable of meta file, for which a prediction model will be built:", choices = ""),
                    selectizeInput("forest_covariable", "If chosen variable has more than 2 groups, choose level in variable which will be compared against the other levels:", choices = ""),
                    plotOutput("forest_sample_preview", height = "200px"),
                    hidden(div(
                      id = "forest_continuous_options",
                      radioButtons("forest_continuous_radio", "If a numeric/continuous variable was chosen, select cutoff value to transform variable into 2 distinct groups:", choices = c("Mean", "Median", "Custom (Use slider below)"), inline = T),
                      sliderInput("forest_continuous_slider", "Custom split", 0, 1, 0, .01),
                      plotOutput("forest_continuous_preview", height = "200px")
                    )),
                    hr(),
                    selectInput("forest_type", "Select mode of model calculation", choices = c("random forest"), selected = "randomForest"),
                    selectInput("forest_features", "Select meta-features to build model", choices = "", multiple = T),
                    checkboxInput("forest_otu", "Use OTU relative abundances to predict model", T)
                  ),
                  wellPanel(
                    p("Model parameters:"),
                    verbatimTextOutput("forest_model_parameters")
                  )
                ),
                column(6, wellPanel(
                  actionButton("forest_toggle_advanced", "Show/hide advanced options"),
                  tags$hr(),
                  hidden(
                    div(
                      id = "forest_advanced",
                      div(
                        id = "general_advanced",
                        checkboxInput("forest_clr", "Perform centered log ratio transformation on OTU abundace data", F),
                        sliderInput("forest_partition", "Choose ratio of dataset which will be used for training; the rest will be used for testing the model", 0, 1, .75, .01),
                        selectizeInput("forest_exclude", "Exclude OTUs from model calculation:", choices = "", selected = NULL, multiple = T),
                        p("Resampling options:"),
                        selectInput("forest_resampling_method", "The resampling method", choices = c("boot", "cv", "LOOCV", "LGOCV", "repeatedcv"), selected = "repeatedcv"),
                        numericInput("forest_cv_fold", "Number of folds in K-fold cross validation/Number of resampling iterations for bootstrapping and leave-group-out cross-validation", min = 1, max = 100, step = 1, value = 10),
                        numericInput("forest_cv_repeats", "Number of repeats (Only applied to repeatedcv)", min = 1, max = 100, value = 3, step = 1)
                      ),
                      tags$hr(),
                      div(
                        id = "ranger_advanced",
                        numericInput("forest_ntrees", "Number of decision trees to grow", value = 500, min = 1, max = 10000, step = 1),
                        textInput("forest_mtry", "Number of variables to possibly split at in each node (multiple entries possible, seperate by comma)", "1,2,3"),
                        selectInput("forest_splitrule", "Splitting rule", choices = c("gini", "extratrees", "hellinger"), selected = "gini", multiple = T),
                        textInput("forest_min_node_size", "Minimal node size (multiple entries possible, seperate by comma)", "1,2,3"),
                        selectInput("forest_importance", "Variable importance mode", choices = c("impurity", "impurity_corrected", "permutation"), selected = "impurity")
                      ),
                      # TODO: remove placeholders!!
                      hidden(
                        div(
                          id = "gbm_advanced",
                          textInput("gbm_ntrees", "Number of decision trees to grow (multiple entries possible)", placeholder = "Enter mutiple numbers seperated by comma.."),
                          textInput("gbm_interaction_depth", "The maximum depth of variable interactions. (multiple entries possible)", placeholder = "Enter mutiple numbers seperated by comma.."),
                          textInput("gbm_shrinkage", "The shrinkage parameter applied to each tree in the expansion (multiple entries possible)", placeholder = "Enter mutiple numbers seperated by comma.."),
                          textInput("gbm_n_minobsinoode", "Integer specifying the minimum number of observations in the trees terminal nodes (multiple entries possible)", placeholder = "Enter mutiple numbers seperated by comma..")
                        )
                      )
                    )
                  ),
                  actionBttn("forest_start", "Start model calculation!", icon = icon("play"), style = "pill", color = "primary", block = T, size = "md"),
                  checkboxInput("forest_default", "Use default parameters: (toggle advanced options for more flexibility)", T)
                ))
              ),
              hr(),
              fixedRow(
                column(5),
                column(2, h1("Results")),
                column(5)
              ),
              fixedRow(
                column(4, wellPanel(
                  p("Confusion Matrix for testing-dataset"),
                  plotOutput("forest_con_matrix"),
                  downloadLink("forest_con_matrixPDF", "Download as PDF"),
                  br(),
                  p("Confusion Matrix for full dataset"),
                  plotOutput("forest_con_matrix_full"),
                  downloadLink("forest_con_matrix_fullPDF", "Download as PDF")
                )),
                column(4, wellPanel(
                  p("ROC-Plot: TP-rate vs. FP-rate including AUC for model"),
                  plotOutput("forest_roc"),
                  downloadLink("forest_rocPDF", "Download as PDF"),
                  p("The receiver operating characteristic (ROC) can show you how good the model can distuingish between sample-groups. A perfect ROC-Curve would go from (0,0) to (0,1) to (1,1). This means the model has a perfect measure of seperability. A ROC-Curve that goes diagonally from (0,0) to (1,1) tells you, that the model makes only random predictions."),
                  p("The AUC (area under the curve) is a good measure to compare multiple ROC curves and therefore models. Here a AUC of 1 tells you, that you have a perfect model, AUC of 0.5 is again only random.")
                )),
                column(4, wellPanel(
                  p("Show the top x most important features for building the model. You can change how many features to display by moving the slider."),
                  sliderInput("top_x_features", "Pick x", min = 1, max = 100, value = 20, step = 1),
                  plotOutput("forest_top_features"),
                  downloadLink("forest_top_featuresPDF", "Download as PDF")
                )),
                downloadButton("forest_save_model", "Save model object as RDS file")
              )
            )
          )
        )
      ),
      ##### info#####
      tabItem(
        tabName = "info",
        h4("Information & global settings"),
        fluidRow(
          tabBox(
            id = "info", width = 12,
            tabPanel(
              "Data Input Format",
              tags$hr(),
              fixedRow(
                column(1, ""),
                column(10, htmlOutput("info_inputdata")),
                column(1)
              )
            ),
            tabPanel(
              "Testdata",
              tags$hr(),
              fixedRow(
                column(1, ""),
                column(10, htmlOutput("info_testdata")),
                column(1)
              )
            )
            # tabPanel(
            #   "Settings",
            #   hr(),
            #   fixedRow(
            #     column(1),
            #     column(8, 
            #            selectInput("namco_pallete","Select global color palette for plots:", choices=c("Set1","Set2","Set3","Paired","Dark2","Accent","Spectral"), selected = "Set2"),
            #            p("You can select one of these color palettes to be applied to (almost) all plots. For details about the palettes see this link: https://ggplot2-book.org/scale-colour.html#brewer-scales"))
            #   )
            # )
          )
        )
      )
    )
  )
)
