## GLINT Data analysis code

Analysis code for GLINT nuller data analysis. Very much work-in-progress.

It contains 3 main scripts (with a bunch of associated functions):

* `glintReadData.m` - Reads the raw acquisition files (binary DAQ streams) and saves appropriately binned data to an intermediate binned data file. Also makes a bunch of useful plots to help choose suitable data ranges, give a preview of the histogram, etc.
* `glintFitData.m` - Fits ASC / NSC models to the binned data files. Note that NSC is only implemented on GPU, so a CUDA-capable GPU is required. Tune `nLoops` and `nSamps` as required. Also does basin-hopping. Saves output as a fitted parameter file.
* `glintReadFittedData.m` - Takes the fitted parameter file and makes a nice plot of (the best) result, complete with fitted parameter values. When multiple fits are in a file (i.e. due to basin-hopping) default behaviour is to plot the best one. Also makes some useful diagnostic plots.