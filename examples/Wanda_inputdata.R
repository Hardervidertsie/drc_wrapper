
# example of the startpars.
startpars = list()
startpars[[1]] <- c(2, 0.1, 1, 0.01)
startpars[[2]] <- c(2, 0, 1.2, 2)

#startpars[[5]] <- c(2, 0.1, 1, 0.01)

#startpars[[23]] <- c(500,1,0.7,2)

# examples use of rmData:

rmData <- data.frame(cell_line = "SRXN1", treatment = "Acrylamide", dose_uM = 10000)

# LL.4 or W2.4 or W1.4
source("D:/src/drc_wrapper/drc_wrapper.R")

# run drc_wrapper function. Result is a list with 3 slots: model.out, ED_est & CurveName
data_result <- drc_wrapper(summaryFile_root_dir = "summaryDataFiles", 
                           featureName = "Cytoplasm_Intensity_IntegratedIntensity_Image_Rhodamine", 
                           controlTreat = "DMSO", nConc = 9, minRespL = NULL, minmaxNorm = FALSE, DEBUG = FALSE, 
                            doseFun = "W1.4", zeroDoseTreatment = "DMSO",
                           EDs = c(10, 25, 50, 80),
                           type = "relative", reference = "control",
                           startpars = startpars, calibrated = TRUE, plotProcData = FALSE,
                           rmData = NULL, lowerl = NULL, upperl = NULL,
                           finalAnalysis = FALSE, test1 = NULL)

# inspect plots and modify
i=0
i=i+1
plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")
data_result$ED_est[[i]]

# to check ED and sd:
ED(data_result$model.out[[i]], respL = 0.5)

# write data to file (after running with calibrated = TRUE )
write.table(data_result$result_ED, file = 'output/ED_cytoplasm_int_integr.txt', sep = '\t', row.names = FALSE)
write.table(data_result$result_SD, file = 'output/ED_SD_cytoplasm_int_integr.txt', sep = '\t', row.names = FALSE)



