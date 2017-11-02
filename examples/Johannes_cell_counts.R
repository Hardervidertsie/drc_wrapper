
#rm(list=ls)
#I would like to create the ECs of the following parameters:
#  -Fraction GFP-positive cells (absolute)
#-Fraction PI-positive cells (absolute)
#-Integrated GFP intenstity (relative)
#-Cell number decline compared to medium control  (relative)
#For comparison with the results of our partners I would like to have EC10, 15 and 50.
#In the data upload sheet, we have to fill in BMCs but I will use EC 10 an 15 for that.



# wrapper function for drc with summary data files as input.

# input arg: root dir with files (character). Measure of interest (character). 
# doseFun (character). ED Dose response derived measures vector of integers.
# controlTreat (character) if provided will divide by this value in treatment column on plate to plate and cell line basis
# minmaxNorm (boolean) if TRUE will perform minmax normalization on plate to plate and cell line basis
# zeroDoseTreatment (character) if any what treatment should be plate-cell_line specific zero concentration?
# startpars list containing 4 parameters in a vector per slot; index of the slot should correspond to curve. Slots can be NULL
# minRespL minimum absoluut response level defined as max-min for the dose response curve. Note if minmaxNorm = TRUE or
# controLTreat the scale of your input values will change accordingly.
# plotProcData plot normalized non fitted data, handy for checking normalization/ data and estimate starting pars
# rmData: data.frame with identical headers and entries to remove. Possible headers are 'cell_line', 'treatment', 'dose_uM', 'plateID', 
# 'replID', 'plateWellID', 'variable'
#   Use to remove e.g. treatmant-dose combinations where GFP response declines due to cell death

## use b c d e parameters for the 4-parameter log logistic functions as described below.
# keep calibrated set to FALSE untill the fits are converging.

# f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
# with b: Hill's slope of the curve (i.e. this is related to the steepness of the curve at point c).
# with c: min value
# with d: max value
# with e: the point of inflection (i.e. the point on the S shaped curve )


# startpars in order b c d e 
# the ED function requires the arguments reference (control/upper), type (relative/absolute) (see ?ED for further arguments)

# step 1: setting drc_wrapper arguments
# set root dir of summary data files
# set feature name of interest.
# set controlTreat if you require fold change with respect to defined controlTreat
# set minmaxNorm if you require minmax normalization
# set zeroDoseTreatment so all curves have a zero concentration
# set nConc to number of concentrations in your dose response curves 
# set minRespL to minimum response level (recommended as will remove non responsive curves)

# step 2: calibration
# define 1 set of start parameters eg: startpars[[1]] <- c(1,1,1,2)
# run with calibrated = FALSE and plotProcData = TRUE
# adjust minRespL startpars or remove data with rmData if numerical optimizations are not converging
# You might want to remove dropping stress response signals at higher concentrations if you know cells are dieing. 
# check your plot in /output (plotProcData = TRUE creates figure of your data) if you change rmData, minRespL or nConc

# step 3: perform ED estimations 
# run with calibrated = TRUE

# step 4: manual inspection
# repeatedly run the lines with i= 0:
#i=i+1
#plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")
#log(data_result$ED_est[[i]], base=10)

# step 5: re-iterated based on results.

# optional: 
# sometimes it is possible to find a full numerical optimization solution for the set of dose-responses (identical min and max ll4 parameters).
# run with finalAnalysis = TRUE
# with plot(data_result$model.out) and summary(data_result$model.out) you can check the results.



# example of the startpars.
startpars = list()

startpars[[1]] <- c(1,1,1,5)
startpars[[4]] <- c(1,1,1,1)
startpars[[10]] <- c(0,0,0,10)
startpars[[15]] <- c(1.5,1,1,1)
startpars[[18]] <- c(1.5,1,1,1)
startpars[[29]] <- c( 10    ,     0.7911     ,    1.2396   ,     20.3251  )
startpars[[33]] <- c( 1    ,     0.7911     ,    1.2396   ,     10.3251  )
startpars[[34]] <- c( 10    ,     0.2911     ,    1.2396   ,     0.53251  )
startpars[[37]] <- c( 1    ,     0.2911     ,    1.2396   ,     0.53251  )
startpars[[48]] <- c( 1    ,     1     ,    2   ,     50  )
startpars[[51]] <- c( 50    ,     0.5     ,    1.5   ,     1  )
startpars[[53]] <- c( 1    ,     0.5     ,    1.5   ,     10  )
startpars[[54]] <- c( 0    ,     0.1     ,    1.5   ,     1  )
# examples use of rmData:

rmData <- data.frame(cell_line = "SRXN1", treatment = "Acrylamide", dose_uM = 10000)


source("drc_wrapper.R")

# run drc_wrapper function. Result is a list with 3 slots: model.out, ED_est & CurveName
data_result <- drc_wrapper(summaryFile_root_dir = "summaryDataFiles", 
                           featureName = "Nuclei_Number_Object_Number", 
                           controlTreat = "Medium", nConc = 10, minRespL = NULL, minmaxNorm = FALSE, DEBUG = FALSE, 
                            doseFun = "LL.4", zeroDoseTreatment = "Medium", 
                           type = "relative", reference = "upper",
                           startpars = startpars, calibrated = TRUE, plotProcData = FALSE,
                           rmData = NULL, lowerl = c(NA, 0, 0, NA) , upperl = c(NA, 5, 5, NA),
                           finalAnalysis = FALSE)

# inspect plots and modify
i=0

i=i+1
plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")
log(data_result$ED_est[[i]], base=10)


# write data to file (after running with calibrated = TRUE )
write.table(data_result$result_ED, file = 'output/ED_cell_counts.txt', sep = '\t', row.names = FALSE)
write.table(data_result$result_ED, file = 'output/ED_SD_cell_counts.txt', sep = '\t', row.names = FALSE)



