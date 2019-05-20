# Code Necessary to Recreate Results

### Dependencies
These are the packages all of the programs imported. You will need to have these installed to run various programs:

* matplotlib
* numpy
* pandas
* hdbscan
* cartopy
* scipy
* os, socket, sys
* gzip
* functools
* math
* datetime

&nbsp;
&nbsp;
&nbsp;

#### Processing Raw Data  !!! *This will need to be done before running any other programs* !!! ####

#####"data_processor.py"#####
You will need to define the location of your raw data (*.gz of ENTLN data--> WWLLN data is pulled from internet automatically when merging). You will also need the location of the list of TGFs with dates and times (UsefulTGFs_With_Dates_and_Times2.csv). Every other parameter should be set to give merged data ENTLN and WWLLN data for +/- 10 minutes, and within 500km of Fermi, for the times and dates provided in the TGF list.  

&nbsp;
&nbsp;
&nbsp;

#### Plotting Individual Maps  ####

##### "testing.py" #####
Change "MyENTLN" to the location of the merged and clustered files. This location should be "[LOCATION OF RAW ENTLN DATA]+ENTLN_Processing/Merged/merged_clustering/". Define the TGFs you want to look at and run to produce maps.

&nbsp;
&nbsp;
&nbsp;

#### Creating Histogram Plots ####

##### "create_histogram_plots.py" #####
This should be ready to run. However, you will need to change the data locations of the clustered ENTLN, clustered WWLLN, and clustered MERGED data on lines (94--108) of the "analysis_functions.py" file. 

&nbsp;
&nbsp;
&nbsp;

#### Testing Significance ####

##### "scrambler.py" #####
You will need to define the location of the clustered MERGED data files lines (12--16). Then the script can be run and should display updates of its progress.


##### "non-TGF_clusters.py" #####
Can be used to test significance using non-TGF-Producing storms. Uses data from "NON_TGF_PRODUCING_SAMPLES2.csv", which lists the non-TGF-producing clusters that met predefined criteria. 
