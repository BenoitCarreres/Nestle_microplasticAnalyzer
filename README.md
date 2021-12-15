# Nestle microplastic analyzer

Tool developed by Société des Produits Nestlé S.A. (Nestlé Research) to process Raman spectra (Horiba) data and identify microplastic particles.

---------------------
## Licence

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
     
---------------------

## About

This script was developed as part of a larger application to identify microplastic raman spectra.
This minimalized version was generated to allow reproducibility of our results published in Nature Scientific Reports.
The necessary data and random forest database are available on Zenodo (link bellow).

## Running the script

This software was built and tested on R 4.0.2. HyperSpec 0.99-20171005 is needed to run properly. To install the required dependencies run the following code:
    
    install.packages(c("ranger","splines","igraph"))
    
    install.packages("https://cran.r-project.org/src/contrib/Archive/hyperSpec/hyperSpec_0.99-20171005.tar.gz", repos=NULL, type="source")


The script should be run with all the necessary files in the same folder. It can be called directly via a terminal using Rscript or within Rstudio. Terminal commands at the script location (linux Bash and Windows cmd/PowerShell:

    Linux:    Rscript ramanMpIdentify.R
    Winows:   C:\R\R-4.0.2\bin\Rscript.exe ramanMpIdentify.R

The script will process of the sample files that can be found in the data folders "milk samples" and "WaterSpiked".
Users are required to download the following files from Zenodo (link bellow) and add them in the same folder as the script:
    
    DB_V2.10.RF.rds
    MilkSpiked_thresholds.csv
    WaterSpiked_thresholds.csv

## Script process

Input files are specifically searched using the information provided in the CSV files. The process imports the data, pre-process the spectra to align them to the database settings, baseline, and then the Random Forest model is used to classify. Once classified, each neighboring spectra are clustered into particles, which are then counted by size groups. Once all the files are processes, a summary table is returned with the number of particles per microplastic class. For more details, see the methods of the article given in the link bellow.


Users can also modify this code to perform other analysis of their own HORIBA generated *.spc files.


---------------------------------------------------

Article available in Nature Scientific Reports:
https://doi.org/10.1038/s41598-021-03458-7

Data and random forest model:
https://doi.org/10.5281/zenodo.5607596
