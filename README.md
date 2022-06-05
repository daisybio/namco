![](app/R/www/Logo.png?raw=true "NamcoLogo")

# Namco
Co-occurrence <b>n</b>etwork <b>a</b>nalysis of <b>m</b>icrobial <b>co</b>mmunities

<b>Namco</b> is a free, open-source R shiny app, that offers multiple analysis tools for studying the microbiome. It includes support for raw fastq sequencing files and pre-processed OTU/ASV abundance tables as input and allows for detailed filtering steps prior to your analysis.

The included analysis tabs start with basic steps like alpha or beta diversity and also include more advanced methods of differential analysis between sample groups, confounder analysis, machine learning, functional prediction and several network implementations to study correlations and connections in your dataset.

You can easily access Namco here: <a href=https://exbio.wzw.tum.de/namco/>https://exbio.wzw.tum.de/namco/</a>

![](app/R/www/namco_workflow_final.png?raw=true "NamcoWorkflow")

### Run Namco locally 

If you are working with sensitive data and you wish to analyse it using Namco on your local machine, we offer a Dockerfile that builds a container including all of Namcos depenencies (due to licensing, you will have to download the USEARCH binary on your own <a href=http://www.drive5.com/usearch/download.html>here</a> and place it in the app/R/data/ directory). 

### Documentation:

If you have questions on how to start your analysis in Namco, check out our Youtube screencast or detailed user-manual! If you encounter issues during the usage of Namco, please post an issue with detailed explanaition in this repository. 

Google document with information on all tabs:
https://docs.google.com/document/d/1A_3oUV7xa7DRmPzZ-J-IIkk5m1b5bPxo59iF9BgBH7I/edit?usp=sharing

Screencast on data-upload:
https://www.youtube.com/embed/dMx2nmXqMfU

### Test dataset

We provide a test dataset, which includes an ASV table, Metadata file, fasta sequences and the phylogentic tree (https://github.com/biomedbigdata/namco/tree/master/app/R/testdata) to test Namco. You can use these input files as template to create your own.

### Citation

Please consider citing us if you used Namco during your analysis:

Namco: A microbiome explorer
Alexander Dietrich, Monica Steffi Matchado, Maximilian Zwiebel, Benjamin Ölke, Michael Lauber, Ilias Lagkouvardos, Jan Baumbach, Dirk Haller, Beate Brandl, Thomas Skurk, Hans Hauner, Sandra Reitmeier, Markus List
bioRxiv 2021.12.15.471754; doi: https://doi.org/10.1101/2021.12.15.471754 
