# drc_wrapper
## lab-specific functionality for the drc package
## wrapper function for drc with summary data files as input.  
 
 
### Input args:  

1. **`summaryFile_root_dir`**: root dir with files (character).   

2. **`doseFun`**: LL.4 for 4-parameter log logistic or W1.4/ W2.4 for 4 par Weibull (character).  
3.  **`featureName`**: Measure of interest (character).   
 
4. **`EDs`** Dose response derived measures vector of numerics.    
5.  **`controlTreat`** (character) if provided will divide by this value on plate to plate and cell line basis.  
6.  **`minmaxNorm`** (boolean) if TRUE will perform minmax normalization on plate to plate and cell line basis.  
7.  **`zeroDoseTreatment`** (character) if any what treatment should be plate-cell_line specific zero concentration?  
8.  **`nConc`**: Number of concentrations in your dose responses, dose-response curves with less than *nConc* are removed.  
9. **`minRespL`** minimum absoluut response level defined as max-min for the dose response curve. Note if *minmaxNorm* = TRUE or
10.  **`controLTreat`** is defined the scale of your input values will change accordingly.  
11.  **`startpars`** list containing 4 parameters start values in a vector per slot; index of the slot should correspond to curve. Slots can be NULL  *DEBUG*: only for debugging script, leave as FALSE. The order is b c d e with b = slope, c = min, d = max, e = inflection point  
12. **`calibrated`**: setting to FALSE is usefull to find parameter solutions one by one and adjust start-values. No dose response estimates are calculated unless set to TRUE (logical)  
13.  **`plotProcData`** Sript is slower, set to FALSE only after modifying filter or normalization steps. plot normalized non fitted data, handy for checking normalization/ data and estimate starting pars (logical) 
14.  **`lowerl`** lower bounds for parameters, passed on to ED function, see drc::ED for details 
15.  **`Upperl`** high bounds for parameters, passed on to ED function, see drc::ED for details 
16. **`rmData`**: data.frame with identical headers and entries to remove. Possible headers are 'cell_line', 'treatment', 'dose_uM', 'plateID', 'replID', 'plateWellID', 'variable' Use to remove e.g. treatmant-dose combinations where GFP response declines due to cell death  
17.   **`test1`** (integer or NULL) usefull to start testing on a single curve  
18.   **`maxOfTime`** (logical). In case you have time course data you can choose to select the max value from each individual time curve.  
19.   **`finalAnalysis`**: attempts drm function with identical c & d parameters and separate = TRUE, performs mixed linear modeling with interaction terms. Will likely only work with very high quality data or low amount of curves. See drm manual for details.    

### Default args:  
```R

summaryFile_root_dir = ".", 
                        doseFun = "LL.4", 
                        featureName = NULL,  
                        EDs = c(10, 25, 50, 80),  
                        controlTreat = NULL, # to divide by  
                        minmaxNorm = FALSE,   
                        zeroDoseTreatment = NULL, # to define for each treatment a zero concentration using this control-treatment  
                        nConc = NULL,  
                        minRespL = NULL,  
                        startpars = NULL,  
                        DEBUG = FALSE,  
                        calibrated = FALSE,  
                        plotProcData = TRUE,  
                        lowerl = NULL, upperl = NULL,  
                        rmData = NULL,  
                        test1 = NULL, 
                        maxOfTime = NULL,
                        finalAnalysis = FALSE, ...)  

```     

### Possible workflow:                         

use b c d e parameters for the 4-parameter log logistic functions as described below.  
keep calibrated set to FALSE untill the fits are converging.  

<img src="https://latex.codecogs.com/gif.latex?f%28x%29%20%3D%20c%20&plus;%20%5Cfrac%7Bd-c%7D%7B1&plus;%5Cexp%28b%28%5Clog%28x%29-%5Clog%28e%29%29%29%7D"></img>


with b: Hill's slope of the curve (i.e. this is related to the steepness of the curve at point c).   
with c: min value  
with d: max value  
with e: the point of inflection (i.e. the point on the S shaped curve )  


Weibull W1.4:  


$ f(x)=c + (d - c)exp(-exp(b(log(x) - log(e)))). $

with b: Hill's slope of the curve (i.e. this is related to the steepness of the curve at point c).  
with c: min value  
with d: max value  
with e: the point of inflection (i.e. the point on the S shaped curve )  
 
startpars in order b c d e   
 
#### step 1: setting drc_wrapper arguments 

set root dir of summary data files  
set feature name of interest.  
set controlTreat if you require fold change with respect to defined controlTreat  
set minmaxNorm if you require minmax normalization  
set zeroDoseTreatment so all curves have a zero concentration  
set nConc to number of concentrations in your dose response curves   
set minRespL to minimum response level (recommended as will remove non responsive curves)  
  
test with `test1 = 1` (or other curve number) and plot fit for manual inspection  
  
#### step 2: calibration  

define 1 set of start parameters eg: `startpars[[1]] <- c(1,1,1,2)`
run with `calibrated = FALSE` and `plotProcData = TRUE`  
adjust minRespL startpars or remove data with rmData if numerical optimizations are not converging  
You might want to remove dropping stress response signals at higher concentrations if you know cells are dying.   
check your plot in /output (`plotProcData = TRUE` creates figure of your data) if you change rmData, minRespL or nConc  

#### step 3: perform ED estimations   

run with `calibrated = TRUE`  
  
#### step 4: manual inspection  

repeatedly run the lines with `i= 0: ` 
```R

i=i+1  
plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")  
log(data_result$ED_est[[i]], base=10)  
```
#### step 5: re-iterate based on results.

*optional:*   
sometimes it is possible to find a full numerical optimization solution for the set of dose-responses (identical min and max ll4 parameters).  
run with `finalAnalysis = TRUE ` 
with `plot(data_result$model.out)` and `summary(data_result$model.out)` you can check the results.  

### example startpars & rmData  

#### examples use of startparts:
```R

startpars = list()  

startpars[[1]] <- c(2,0.1,1,0.008)  
startpars[[13]] <- c(10,1,0.7,0.1)  
startpars[[15]] <- c(10,0.7,1,0.1)  
startpars[[23]] <- c(500,1,0.7,2)  
``` 

#### examples use of rmData:
```R 

rmData <- data.frame(cell_line = "SRXN1", treatment = "Acrylamide", dose_uM = 10000)  
 
# LL.4 or W2.4 or W1.4  
source("~/src/drc_wrapper.R")  
  
# run drc_wrapper function.
# Result is a list with 3 slots: model.out, ED_est & CurveName  

data_result <- drc_wrapper(summaryFile_root_dir = "summaryDataFiles",   
                           featureName = "Cytoplasm_Intensity_IntegratedIntensity_Image_Rhodamine",   
                           controlTreat = "DMSO", nConc = 10, minRespL = NULL,
                            minmaxNorm = FALSE, DEBUG = FALSE,   
                            doseFun = "W1.4", zeroDoseTreatment = "DMSO",   
                           type = "relative", reference = "control",  
                           startpars = NULL, calibrated = FALSE, plotProcData = FALSE,  
                           rmData = NULL, lowerl = NULL , upperl = NULL,  
                           finalAnalysis = FALSE, test1 = 2)  
  
# inspect plots and modify  

i=0  
  
i=i+1  
plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")  
data_result$ED_est[[i]]  
  
#to check ED and sd:  

ED(data_result$model.out[[i]], respL = 0.5)  
  
# write data to file (after running with calibrated = TRUE )  
write.table(data_result$result_ED,
            file = 'output/ED_cytoplasm_int_integr.txt',
            sep = '\t', row.names = FALSE)  

write.table(data_result$result_ED,
            file = 'output/ED_SD_cytoplasm_int_integr.txt',
            sep = '\t', row.names = FALSE)  
```

