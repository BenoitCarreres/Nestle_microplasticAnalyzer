# Nestle microplastic analyzer
Tool developed by Société des Produits Nestlé S.A. (Nestlé Research) to process Raman spectra (Horiba) data and identify microplastic particles.


This tool was developed as part of a larger application to identify microplastic raman spectra.
This minimalized version was generated to allow reproducibility of our results published in Nature Scientific Reports.
The necessary data and random forest database are available on Zenodo (link bellow). The R script should take an input csv file that will provide the necessary information to analyzed the published dataset.


Users are required to download the following files from Zenodo (link bellow) and add them in the same folder as the script: "DB_V2.10.RF.rds", "MilkSpiked_thresholds.csv", "WaterSpiked_thresholds.csv".

Users can also modify this code to perform other analysis of their own HORIBA generated *.spc files.

This software was built and tested on R 4.0.2

Article available in Nature Scientific Reports:
-link-

Data and random forest model at Zenodo:
10.5281/zenodo.5607597
