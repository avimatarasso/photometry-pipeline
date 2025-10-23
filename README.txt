All code below was adapted by Avi Matarasso from the Bruchas lab pipeline, originally developed by Christian Pedersen.

If you are recording from a TdT machine, you will initiate conversion of TDT files to matlab data and .mat files from the tdtExtract.m script.

Next, you can run your .mat files through the photometryScript.m file to process the .mat files and find individual and averaged traces, as well as heat maps.

tdtExtract - Initial extract script of all tdt files from the session folders. 

photometryScript - Script to actually graph preprocessed tdt files. After initial extraction, use this script.

Other functions and scripts:
alignEvent - lines up each event within the full trace defined in a text file where the first column is the event times and the second column is the event label, separated by a space or a comma. Will output the data in a vector defined by the timeBefore and timeAfter variable, with the timings of the events in timingIdxs.

centAndNormData - centers the data around the baseline, defined here as the timeBefore an event you define in the photometryScript file. Note, this is not normalizing, but you can adjust this code to normalize the data here. 

createPlots - plots individual, across mouse, time-averaged plots and heatmaps of photometry data defined of the eventLabel defined within timeBefore to timeAfter the event. Also includes some full trial plots. 

defineDir - define the files within the working directory that you are planning to work on.

defineSubject - using the photometry file name, region defined, and subject index, this code will define the names of each subject and text files. This code is very individualized. My typical way of naming mice is CAGE#-MOUSE#_Sensor_Region_date_time, or as an example '659-5_GRABDA_CA1-251001-162616'.

finalSave - save all the variables and figures made in the photometryScript.

findBaselineAndNormalize - Not currently used in the photometryScript. the code will find the baseline within one event and normalize the data to that baseline.

fitAndSub - Apply a linear least squared fit of the isobestic to the sensor, and subtract the sensor signal from the isobestic. The function also plots these time series.

headEntryAdjust - If you have a repeated event and you only want the first event after another event (i.e. - head entry after a reward) within 10s (default), include headEntry in the savename string. This script will adjust the events to just look at the one immediately after each reward. 

LLS - Actual calculation of the linear least squared fit of the isobestic to the sensor.

mseb - handy function to graph both the avged time series and the SEM. Credit: Andreas Trier Poulsen, atpo@dtu.dk

plotIndividual - plots individual average time series and heat maps.

tdt2mat - Import data from TDT system recording into Matlab structure

tdtCONVERT - extract all tdt files and data streams from their proprietary structures.

graphing variables and code
brewermap is a nice way to visualize colorschemes and palettes.
cmapVar is a mat file that contains some nice color axis for heatmap plotting (BkBl for black to blue, BkPl for black to purple, and BkRd for black to red.

