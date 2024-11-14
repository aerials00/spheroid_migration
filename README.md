# spheroid_migration
Image Processing code related to the manuscript: https://www.biorxiv.org/content/10.1101/2024.05.08.593120.abstract

These MATLAB codes are used to plot data in Fig.2, Fig.3 Fig.S7, Fig.S8, Fig.S9, Fig.S19, Fig.S20, Fig.S21 and Fig.S22 in the manuscript. 

MATLAB Code 1: effective-spheroid-radius_peak-height_cell-count
This MATLAB code identifies the spheroid boundary and measures the effective radius over time (i.e. for each time point). 
Additionally, it identifies individual cells that have disseminated from the primary spheroid and when the cell is not connected
to the spheroid or present in a close vicinity. 
Lastly, we perform a cartesian to polar transformation of each image to measure the spheroid protrusions as single peaks and 
measure the average height. 

General instructions:
1. Copy the matlab code (.m file) in the folder where images are stored. 
Make sure images are named sequentially. For eg: S1_t01, S1_t02 and so on. 

2. Run the code in MATLAB. 
Copy the folder director in the field 'datadir ='. 
Copy the correct the file initial with file extension, for eg. S1_*.tif 

This will read all files in your folder for further analysis. 

3. Follow the instructions in the code that will guide you through adjusting the image analysis parameters.
These are listed as comments in front of every line of code (when needed). 

4. Run this code for each spheroid case and save the appropriate .mat file in the right order. 
The .mat file will contain the effective spheroid radius, cell count, peak height and many other parameters used for plotting. 


MATLAB code 2: plotting normalized spheroid radius for all spheroids and calculating onset time 
This MATLAB code compiles the data for each spheroid analysis performed using MATLAB code 1 to plot the normalized increase in effective
spheroid radius (with standard deviation) for one particular set of spheroids (For eg. MV3 - control, n = 5 spheroids). 

From this, the code calculates the time at which the spheroid achieves 10% increase from its initial value (i.e. 1.1) referred to as onset time. 
Onset time for each spheroid is saved in the .mat file that can be used to plot. 

Change the directory and filename to plot together different conditions. For eg. MV3_control, MV3_MMP, MV3_TGF-B for 2.4 mg/mL collagen. 
