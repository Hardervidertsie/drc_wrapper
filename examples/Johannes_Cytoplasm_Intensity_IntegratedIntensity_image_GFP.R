

# example of the startpars.
startpars = list()

startpars[[1]] <- c(-9,0.01,0.5,1)

startpars[[4]] <- c(-15,0.01,0.25, 500 )
startpars[[5]] <- c(-1,0.01,0.8,1000)
startpars[[7]] <- c(-15,0,0.5,100)




# examples use of rmData:

rmData <- data.frame(cell_line = "SRXN1", treatment = "Acrylamide", dose_uM = 10000)


source("D:/src/drc_wrapper/drc_wrapper.R")

# run drc_wrapper function. Result is a list with 3 slots: model.out, ED_est & CurveName
data_result <- drc_wrapper(summaryFile_root_dir = "summaryDataFiles", 
                           featureName = "Cytoplasm_Intensity_IntegratedIntensity_image_GFP", 
                           EDs = c(10, 25, 50, 80),
                           controlTreat = NULL, nConc = 10, minRespL = 0.2, minmaxNorm = TRUE, DEBUG = FALSE, 
                            doseFun = "LL.4", zeroDoseTreatment = "DMSO[0.2%]", 
                           type = "relative", reference = "control",
                           startpars = startpars, calibrated = FALSE, plotProcData = FALSE,
                           rmData = rmData, lowerl = NULL , upperl = NULL,
                           finalAnalysis = FALSE, test1 = NULL, maxOfTime =FALSE)

# inspect plots and modify
i=0


i=i+1
plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")
log(data_result$ED_est[[i]], base=10)


# write data to file (after running with calibrated = TRUE )
write.table(data_result$result_ED, file = 'output/ED_cytoplasm_int_integr.txt', sep = '\t', row.names = FALSE)
write.table(data_result$result_ED, file = 'output/ED_SD_cytoplasm_int_integr.txt', sep = '\t', row.names = FALSE)



