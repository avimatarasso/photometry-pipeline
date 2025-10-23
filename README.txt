If you are recording from a TdT machine, you will initiate conversion of TDT files to matlab data and .mat files from the tdtExtract.m script.

Next, you can run your .mat files through the photometryScript file to process the .mat files and find individual and averaged traces, as well as heat maps.

alignEvent - lines up each event within the full trace defined in a text file where the first column is the event times and the second column is the event label, separated by a space or a comma. Will output the data in a vector defined by the timeBefore and timeAfter variable, with the timings of the events in timingIdxs.

centAndNormData - centers the data around the baseline, defined here as the timeBefore you define in the photometryScript file. you can uncomment some lines of code to find what happens when you normalize to the standard deviation of the baseline too.

