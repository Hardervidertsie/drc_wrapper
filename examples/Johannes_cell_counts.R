


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
                           startpars = NULL, calibrated = FALSE, plotProcData = FALSE,
                           rmData = NULL, lowerl = c(NA, 0, 0, NA) , upperl = c(NA, 5, 5, NA),
                           test1 = NULL, maxOfTime = FALSE,
                           finalAnalysis = FALSE)

# inspect plots and modify
i=0

i=i+1
plot(data_result$model.out[[i]],  main = data_result$CurveName[[i]], type = "all")
log(data_result$ED_est[[i]], base=10)


# write data to file (after running with calibrated = TRUE )
write.table(data_result$result_ED, file = 'output/ED_cell_counts.txt', sep = '\t', row.names = FALSE)
write.table(data_result$result_ED, file = 'output/ED_SD_cell_counts.txt', sep = '\t', row.names = FALSE)



