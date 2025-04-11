
(1) Citation to the manuscript: This document describes the data files and R-code that accompany ``The effects of resources and population expansion on female-male protein consumption among hunter-gatherers" by Jacob Freeman, Alyssa Ahmann, Robert J. Hard, and Raymond P. Mauldin. The manuscript has been submitted to the Journal of Archaeological Science for re-review on April 11, 2025.  The data used in the paper are contained in 3 data.csv files and the code to analyze the data is contained in SexDifPopResource24.v2.R. This README document describes the csv files and how to get started with the analysis. 

2025 Freeman, J., A. Ahmann, R.J. Hard, R.P. Mauldin. The effects of resources and population expansion on female-male protein consumption among hunter-gatherers, Journal of Archaeological Science, Under Review. 

(2) Description of the files included in the repository and their relationship to the figures and tables in the manuscript
Files: 


(i) File path: data/TexasHGIsotopesMaleFemale.csv 

This file contains information on the stable isotopes of individuals buried in Central Texas and on the Texas Coastal Plain. The files is essential to reproduce the
analysis presented in the main manuscript, including Figures 3, 4, 5, and 6.  

Site--Archaeological site name
Individual ID (paper)--Provenience of individuals at a site	
delta13cap--delta 13C apatitive values	
delta13ccol--delata 13C collagen values	
delta15n--delata 15N collagen values	
C_N--collagen to nitrogen ratio
MidPoint--median date of radiocarbon age on an individual or the median date of associated component for individuals	
PeriodName--Cultural historical time period
PeriodID--Numeric ID for each cultural historical time period.
Age--Estimated age of individual at time of death, as gathered by source material. 
sex--Estimate of male or female for individuals. U indicates unidentified. M=male, F=female.	
Estimated Latitude--Decimal degrees north, obscured
Estimated Longitude--Decimal degrees west, obscured
Body part--Body part (teeth or bone) used for radiocarbon dating and isotope analysis.
Source--source of data	
Setting--Central Texas =CTx, Coastal Plain=TCP
Setting2--Id of inland, riverine and coastal ecosystems


(ii) File path: data/TCPPerCap.csv 

This files contains the time-series of the per capita growth rates of all simulated KDEs, the mean KDE, and per capita growth rate of the mean KDE. The file is essential to reproduce Figures 4 and 6 in the main manuscript and Figures 2 and 3 in the Supporting Material. This file is for those who do not want to reproduce the full analysis of archaeological radiocarbon based on the raw data discussed below.

calBP--30 year time-steps in cal BP.
PeriodID--Cultural historical time period	
PeriodID2--Numeric ID for each cultural historical time period.
v1-v200: The per capita growth rate of each KDE in 30 year intervals
MKDE--The mean value of the 200 KDEs at each 30 year time-step
PerCap--The per capita growth rate of the mean KDE

(iii) File path: data/FinalRCDTexas3.csv 

This file contains all radiocarbon ages associated with archaeological remains collected from Central Texas and the Texas Coastal Plain. This file is essential to reproduce 
TCPPerCap.csv and, ultimately, Figures 4 and 6 in the main manuscript and Figures 2 and 3 in the Supporting Material. These data are provided for researchers interested in reproducing the KDEs used in the study and engaging in their own analysis.

Latitude-decimal degrees north
Longitude-decimal degrees west
Site Name--name of archaeological site	
Trinomial--unique trinomial for each site	
Assay No.--unique id for lab and each radiocarbon sample	
Provenience--provenience of dated material if known
Feature Type--archaeological feature from which material was recovered
Material--material dated	
Taxa dated-species or genera of data material	
Human--yes or no. yes indicates human remains dated	
Age-- radiocarbon age	
Error--standard deviation of raw radiocarbon age	
Corrected/Raw	
Corrected/Normalized 14C age	 
Corrected/Normalized Age	
Delta 13C	
Delta 13C source	
Raw/Measured Age	 
Raw/Measured Age	
AMS or Radiometric	
Comments	
Reference
Region--Central Texas =CTx, Coastal Plain=TCP
Zone--Inland, riverine or coastal ecosystems

(3)The code to analyze the data is contained in SexDifPopResource24.v2.R. The analysis was run in "R version 4.2.2 (2022-10-31 ucrt)" with the Rstudio integrated development environment version ``2023.06.0 Build 421."


(4) Getting started working with the project:

To begin working with the data, (i) open SexDifPopResource24.v2.R., (ii) install and load all of the relevant packages at the beginning of the document, and (iii) follow the comments to begin reproducing the results reported in the manuscript. 
