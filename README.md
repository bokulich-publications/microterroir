# microterroir

All notebooks to analyse data generated within the Lavaux project.


 ## Sensory Analysis 
 `Sensory_Analysis.ipynb`

 Data from tasting with 12 judges of the 2021 microvinification wines according to given sensory characteristics (olfactory and taste). 

 * 2-way ANOVA: primarily significant differences between the judges 
 * PCA: visually check scatter between samples 
 * Box Plots of characteristics
 * Get median sensory values for each plot --> metadata in future analysis! 
 * Spider plots showing distribution of characteristics of each plot 


## GC-MS 
`GC-MS.ipynb`

VOCs analysis with HS-SPME-GC-MS of all 2023 microvinification samples and the end-point fermentation samples (post-MLF) from 2022 and 2021.

* Merge, Normalize and Log2 transform GC-MS raw data 
* Compare MLF samples from different years: 
  * Heatmaps with Z score, euclidian clustering, top variable compounds
  * PCA plot
  * Compare MLF samples with Plot metadata
    * PERMANOVA 
    * Mantel Test with geodisic distances between plots
* Compare 2023 samples 
  * Heatmap
  * PCA plot
* Compare all samples 
  * PCA
  * Diversity Analysis with Richness and Shannon 
  * Statistical differences between groups with ANOVA (box plots)
  
## Correlate GC-MS and Sensory Data
`Correlate_GC_Sensory.ipynb`

* Canonical Correlation Analysis
* PLSR analysis
  * Determine number of components 
  * PLSR with various number of components 
  * Biplot with major loadings 


## Berry Chemistry --> WORK IN PROGRESS 
`Berry_Chemistry.ipynb`

HPLC data of 2021, 2022, 2023 berries 
* Swarmplots of major organic acids and sugars
* ANOVA 
  * Nested (`~ Plot_No * Date`)
  * with Plot Metadata 


## Sensor Data --> WORK IN PROGRESS 
Data from Temperature and Relative Humidity Sensors installed in each plot for 3 years. 

`Sensor_Data_Temperature.ipynb`

* Remove outliers (e.g. when battery was broken or sensor fell to the ground) with Median absolute deviation 
* Compare the median, max and min temperatures of the years
  

`Sensor_Data_RH.ipynb`


## LC-MS 
`LC-MS_pilot.ipynb`

* Venn Diagrams of the LC-MS pilot 


## Soil Analyis 
`Soil_Analysis.ipynb`

Analysis of the soil characteristics of each plot. Data from Soil Conseil (from 10 plots in 2023) and the soil pH measurements. 

* Soil Conseil 
  * Soil texture: stacked bar plots 
  * Swarm plots: Organic Matter, N total 
  * PCA 
* Soil pH 
  * Lineplots over the sampling time 
  * ANOVA between sampling locations and time points 