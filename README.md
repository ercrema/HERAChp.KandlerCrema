# Introduction
This package contains a set of R functions that enables the reproduction of the analysis and associated figures for the book chapter *Analysing cultural frequency data: neutral theory* by Anne Kandler and Enrico Crema for the volume *Handbook of Evolutionary Research in Archaeology*, edited by Anna Prentiss. The package contains two main functions for simulating cultural transmission:

* `transmission()` simulates unbiased and frequency biased cultural transmission and computes turn-over rates and diversity statistic.
* `heteroPopTransmission()` simulates unbiased and frequency biased cultural transmission with an heterogeneous population of social learners.

More details can be found in the help documentation of each function (which can be accessed by typing `?transmission` and `?heteroPopTransmission`) and on the paper. 

# Installation and Vignette

*HERAChp.KandlerCrema* can be installed directly from the github repository using the following command:

`devtools::install_github("ercrema/HERAChp.KandlerCrema")`

Notice that this requires the `devtools` package installed. 

The package contains a vignette which details the workflow required to generate the figures on the book chapter. Some of the figures require a considerable number of simulation runs. In order to reduce running the document is based on 100 simulation runs (rather than the 1,000 used in the book chapter) and some scripts are executed across multiple threads. Nonetheless, rendering the Rmarkdown file and producing the vignette took ca. 40 minutes on an IMac with 8 logical CPUs (3.3 GHz Intel Core i7) and 16 GB Ram. 

There are three ways to obtain the vignette:

1. Create (render) the vignette during the installation. To do so use the following command during the installation:

`devtools::install_github("ercrema/HERAChp.KandlerCrema", build_vignettes=TRUE)`

and access the vignette with the following command after loading the library:

`vignette("ChapterFigures")`

Notice that this option might take a considerable amount of time (see above) and there are no possibility to change the simulation parameters. 


2. Manually download and render the markdown file.

Download the markdown file by clicking [here](https://raw.githubusercontent.com/ercrema/HERAChp.KandlerCrema/master/vignette/ChapterFigures.Rmd),and saving the page ensuring that the file extension is kept to `.Rmd`. Move the file to your current working directory, and type the following command in R:

`rmarkdown::render("ChapterFigures.Rmd")`

This option allows you to change the settings of the R markdown file (e.g. changing the number of simulation runs).

3. Download the rendered vignette file by clicking [here](https://raw.githubusercontent.com/ercrema/HERAChp.KandlerCrema/master/vignette/ChapterFigures.html) and saving the opened page as an html file on your local machine. The file can then be opened by any browser. This option does not require package installation.





