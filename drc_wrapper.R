drc_wrapper <- function(summaryFile_root_dir = ".",   
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
                        maxOfTime = FALSE,
                        finalAnalysis = FALSE, ...) {
  
  if(DEBUG){
    summaryFile_root_dir <- "summaryDataFiles"
    doseFun <- "LL.4"
    featureName <- "Cytoplasm_Intensity_IntegratedIntensity_Image_Rhodamine"	
    EDs <- c(10,25,50, 80)
    controlTreat <- "DMSO"
    zeroDoseTreatment <- "DMSO"
    minmaxNorm <- FALSE
    nConc <- 10
    
    plotProcData <- FALSE
    rmData <- NULL
    
  }  
  options(stringsAsFactors = FALSE)
  # load require packages
  require(drc)
  require(reshape2)
  require(magic)
  require(ggplot2)
  if(grep("4", doseFun) != 1) {
    stop("Currently only 4-parameter models implemented in drc_wrapper")
  }
  
  if(!is.null(test1) & calibrated){
    stop("when testing set calibrated to FALSE")
  }
  
  if(!is.null(test1) & !is.numeric(test1)){
    stop("test1 must be integer")
  }
  if(is.numeric(test1)){
    if(round(test1, digits=0) != test1){
      stop("test1 must be integer")
    }
  }
  
  if(!is.null(startpars)){
    if(!is.list(startpars)){
      stop("startpars is not list")
    }
  }
  
  if(
    any( unlist(sapply(startpars, "[[", 2)) > upperl[2]) | any(unlist(sapply(startpars, "[[", 2)) < lowerl[2])){
    stop("start values not within specified lower and upper limits")
  }
  
  if(
    any( unlist(sapply(startpars, "[[", 3)) > upperl[3]) | any(unlist(sapply(startpars, "[[", 3)) < lowerl[3])){
    stop("start values not within specified lower and upper limits")
  }
  
  
  if(!dir.exists('drc_output')){
    dir.create('drc_output')
    message("\`drc_output\` directory created")
  }
  if(is.null(featureName)) {
    stop("please provide featureName")
  }
  if(!is.logical(minmaxNorm)){
    stop("expected logical minmaxNorm")
  }
  
  # check input
  if(!is.null(zeroDoseTreatment) & !is.null(controlTreat)){
    warning("Setting zero dose treatment of fold changes")
  }
  if(minmaxNorm & !is.null(controlTreat)){
    stop("Minmax normalization of fold changes is not allowed, 
         choose either minmaxNorm or control, not both")
  }
  if(is.null(nConc)){
    stop("Please specify number of concentrations")
  }
  
  if(!doseFun %in% unlist(getMeanFunctions(display = FALSE))) {
    stop(paste(doseFun, " not found, please run \"getMeanFunctions()\" for available models, example: \"LL2.4\""))
  }
  
  
  file_list <- dir(summaryFile_root_dir)[grepl(".txt", dir(summaryFile_root_dir))]
  if( length(file_list) == 0 ) {
    stop("no files in provided root dir")
  }
  data_list = alist()
  filepath_list <- paste(summaryFile_root_dir, file_list, sep ='/')
  
  for(i in seq_along( file_list) ) {
    data_list[[i]] <- read.delim(file = filepath_list[i], sep ="\t", header = TRUE)
  }
  
  my_data <- do.call('rbind', data_list)
  
  
  
  
  
  
  # validate input arguments and data
  if(!all(c('cell_line', 'treatment', 'dose_uM', 'plateID', 
            'replID', 'variable', 'value') %in%
          colnames(my_data))
  ) {
    stop("drc_wrapper requires the following headers:\n
         'cell_line', 'treatment', 'dose_uM', 'plateID', 
         'replID', 'variable', 'value' ")
  }
  
  if(!is.logical(maxOfTime)){
    stop("logical maxOfTime expected")
  }
  
  if(!maxOfTime) {
  if( !is.null(my_data$timeID) | !is.null(my_data$timeAfterExposure) ) {
    if(length(unique(my_data$timeID)) > 1 | 
       length(unique(my_data$timeAfterExposure)) > 1 ) {
      warning("Multiple unique time points detected, 
              dose response modeling input is not suited for time course data, 
              please aggregate first.")
      }
    }
  }
  
  # aggregate over time
  if(maxOfTime){
    
    my_data <- aggregate(value ~ cell_line + treatment + dose_uM + plateID + 
                         replID + variable, data = my_data, FUN = max )
  
print(writeLines("aggregated: (value ~ cell_line + treatment + dose_uM + plateID + 
                         replID + variable, FUN = max) \n 
                 consider performing time-removing aggregation yourself")
      )
warning("aggregated over time. Consider performing time-removing aggregation yourself")
    }
  
  
  
  if(!is.null(controlTreat)) {
    if(!controlTreat %in% my_data$treatment) {
      stop("Specified controlTreat not found in treatment entries")
    } 
  }
  
  if(!featureName %in% my_data$variable) {
    stop("Specified featureName not found in variable entries")
  } 
  
  my_data <- my_data[my_data$variable %in% featureName, ]
  
  
  
  if(!is.null(zeroDoseTreatment)) {
    # add dmso as start concentration for all compounds
    if(!zeroDoseTreatment %in% my_data$treatment) {
      stop("zeroDoseTreatment not found in treatment entries")
    }
    
    
    all_treats <- unique(my_data$treatment)
    
    ii=1;i=2
    outputlist = alist()
    outputcomps = alist()
    my_data$tmp_plateID_cell_line <- paste0(my_data$plateID, my_data$cell_line)
    plateCellIDs <- unique(my_data$tmp_plateID_cell_line)
    zerDoseData <- my_data[ my_data$treatment == zeroDoseTreatment, ]
    
    for( ii in seq_along(plateCellIDs)){
      buffer <- zerDoseData[ zerDoseData$tmp_plateID_cell_line == plateCellIDs[ii] ,]
      for(i in seq_along(all_treats)){
        buffer$treatment <- all_treats[i]
        if(i == 1){
          compresult <- buffer
        }  
        if(i != 1){
          compresult <- rbind(compresult,buffer)
        }
      }
      outputlist[[ii]] <- compresult
    }
    zerDoseData <- do.call('rbind' ,outputlist)
    my_data <- rbind(my_data, zerDoseData)
  }
  # if relative to controlTreat TRUE: divide by  
  if(!is.null(controlTreat)){
    
    control_data <- my_data[ my_data$treatment == controlTreat, ]
    
    control_data_aggr <- aggregate( data = control_data, 
                                    value ~ plateID + replID + cell_line,
                                    mean )
    colnames(control_data_aggr)[ colnames(control_data_aggr) == "value"] <- "controlValue"
    
    my_data <- merge(my_data , control_data_aggr, 
                     by = c('plateID', 'replID', 'cell_line' ),
                     all.x = TRUE, all.y=FALSE)
    my_data$value <- my_data$value / my_data$controlValue
    my_data$controlValue <- NULL
    my_data$variable <- paste0(my_data$variable, "vs", controlTreat )
  }
  
  my_data <- aggregate(data = my_data, 
                       value ~  treatment + dose_uM + plateID + replID + variable + cell_line ,
                       mean)
  
  
  if( minmaxNorm ) { 
    
    normFun <- function(dataIn) {
      
      minV <- aggregate(data = dataIn, value ~ plateID + cell_line, 
                        FUN = function(x) min(x, na.rm = TRUE))
      maxV <- aggregate(data = dataIn, value ~ plateID + cell_line, 
                        FUN = function(x) max(x, na.rm = TRUE))
      
      colnames(minV)[colnames(minV) == "value"] <- "minV"
      colnames(maxV)[colnames(maxV) == "value"] <- "maxV"
      
      dataIn <- merge(dataIn, minV, by = c("plateID", "cell_line"), all.x = TRUE)
      dataIn <- merge(dataIn, maxV, by = c("plateID", "cell_line"), all.x = TRUE)
      
      dataIn$norm_value <- (dataIn$value - dataIn$minV) / (dataIn$maxV - dataIn$minV)
      dataIn$minV <- NULL
      dataIn$maxV <- NULL
      dataIn$value <- NULL
      colnames(dataIn)[colnames(dataIn) == "norm_value"] <- "value"
      dataIn$variable <- paste0("mmn_", dataIn$variable)
      return(dataIn)  
    }
    my_data <- normFun(my_data)
    
  }
  
  
  # dose response modelin:
  # check # concentrations per treatment
  
  my_data$CurveName <- paste(my_data$plateID, my_data$cell_line, my_data$treatment)
  
  indrm <- table(my_data$CurveName) < 10
  rmCurves <- names(table(my_data$CurveName))[indrm]
  
  if(sum(indrm)!=0) {
    write.table(rmCurves, file = 'drc_output/removed_curves.txt', sep ='\t', col.names = NA)
    warning("Curves removed becuase n-concentrations < specified number, 
            consult \'drc_output/removed_curves/txt\'")
  }
  
  my_data <-  my_data[ !my_data$CurveName %in% rmCurves, ]
  
  
  
  my_data$CurveName <- paste(my_data$cell_line, my_data$treatment)
  
  my_data_a <- aggregate(value ~treatment + dose_uM + replID + cell_line + CurveName, data = my_data, mean )
  
  all.curves <- unique(my_data_a$CurveName)
  
  if(!is.null(rmData)){
    
    if(!all(colnames(rmData) %in% colnames(my_data_a))  ){
      stop("rmData headers not found in aggregated data")
    }
    
    
    my_data_a$dose_uM <- as.character(my_data_a$dose_uM)
    indToRm <- match( colnames(rmData) ,colnames(my_data_a) )
    rmData$rmNames<- apply(rmData, 1, function(x) {paste(as.character(x),collapse = "_" )})
    my_data_a$rmNames <- apply(my_data_a[, indToRm], 1, function(x) {paste(as.character(x), collapse = "_")})
    my_data_a <- my_data_a[!my_data_a$rmNames %in% rmData$rmNames, ]
    
    all.curves <- unique(my_data_a$CurveName)
    my_data_a$dose_uM <- as.numeric(my_data_a$dose_uM)
  }
  
  
  if(!is.null(minRespL)){
    
    # code to remove non responsive curves here
    
    respdata <- aggregate( value ~ CurveName, data = my_data_a, 
                           FUN = function(x){
                             max(x) - min(x)
                           })  
    rmCurvesResp <- respdata$CurveName[respdata$value < minRespL]
    my_data_a <- my_data_a[!my_data_a$CurveName %in% rmCurvesResp, ]
    all.curves <- unique(my_data_a$CurveName)
  }
  
  
  # plot non fitted result
  if(plotProcData ){
    
    pdfsize <- 4 + round(length(all.curves)/3, digits=0)
    
    pdf(file = paste0('drc_output/', featureName, '_allCurves.pdf'), height = pdfsize, width = pdfsize)
    
    p<-ggplot(data = my_data_a, aes( x = log(dose_uM+0.0001), y = value ))  +
      geom_point(aes(color = factor(replID))) + facet_wrap(~ CurveName) + xlab("log(dose_uM + 0.0001)")
    print(p)
    dev.off()
  }
  
  
  
  if(!is.null(startpars)){
  if(length(startpars) == 1){
    startpars[1:length(all.curves)] <- startpars[1] 
  } else{
    startpars[[length(all.curves)+1]] =list()
    startpars[sapply(startpars, is.null)] <- startpars[1]
    startpars[[length(all.curves)+1]] =NULL
  }
  }
  
  #slope min max inflection
  if(!is.null(test1) & !is.null(startpars)) {
    model.out = alist()
    i <- test1
    sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
    print(paste(i, ": ", all.curves[i]))
    model.out[[i]]<- eval(
      parse( text =
               paste("drm(value ~ dose_uM, type = \"continuous\", start = startpars[[i]],
                     lowerl = lowerl, upperl = upperl, separate = TRUE,
                     data = sel_data, fct =", paste0(doseFun, '()', ")"))))
    print(model.out[[i]])
  } 
  
  if(!is.null(test1) & is.null(startpars)){
    model.out = alist()
    i <- test1
    sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
    print(paste(i, ": ", all.curves[i]))
    model.out[[i]]<- eval(
      parse( text =
               paste("drm(value ~ dose_uM, type = \"continuous\",
                     lowerl = lowerl, upperl = upperl, separate = TRUE,
                     data = sel_data, fct =", paste0(doseFun, '()', ")"))))
    print(model.out[[i]])
  }
  
  if(!calibrated & is.null(test1)){
   if(!is.null(startpars)){ 
    
    model.out = alist()
    
    for(i in seq_along(all.curves)){
      
      sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
      print(paste(i, ": ", all.curves[i]))
      model.out[[i]]<- eval(
        parse( text =
                 paste("drm(value ~ dose_uM, type = \"continuous\", start = startpars[[i]],
                       lowerl = lowerl, upperl = upperl, separate = TRUE,
                       data = sel_data, fct =", paste0(doseFun, '()', ")")
        )
        )
                 )
      print(model.out[[i]])
      
      
    }
   } else{
     model.out = alist()
     
     for(i in seq_along(all.curves)){
       
       sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
       print(paste(i, ": ", all.curves[i]))
       model.out[[i]]<- eval(
         parse( text =
                  paste("drm(value ~ dose_uM, type = \"continuous\", 
                       lowerl = lowerl, upperl = upperl, separate = TRUE,
                       data = sel_data, fct =", paste0(doseFun, '()', ")")
                  )
         )
       )
       print(model.out[[i]])
       
     }
   }
    
  }
  
  #   }
  # )
  # all.curves <- unique(my_data_a$CurveName)
  # 
  # 
  # i=8
  # doseFun <- "LL.4"
  # startpars <- c(0.1, 0.1, 0.1, 2)
  # 
  # plot( sel_data[, "dose_uM"], sel_data[, "value"]  )
  # 
  # for(i in seq_along(all.curves)){
  # print(i)
  #     sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
  #     sel_data$dose_uM <- log(sel_data$dose_uM + 1, base = 10)
  #   eval(
  #     parse( text =
  #              paste("drm(value ~ dose_uM, type = \"continuous\", logDose = 10, start = startpars,
  #                    data = sel_data, separate = TRUE, fct =", paste0(doseFun, '()', ")")
  #                    )
  #     )
  #     )
  # }
  
  
  
  if(calibrated){
    
    
    
    #==============
    ED_est = list()
    
    model.out = alist()
    if(!is.null(startpars)){
      
    for(i in seq_along(all.curves)){
      
      sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
      
      model.out[[i]]<- eval(
        parse( text =
                 paste("drm(value ~ dose_uM, type = \"continuous\", start = startpars[[i]],
                       lowerl = lowerl, upperl = upperl, separate = TRUE,
                       data = sel_data, fct =", paste0(doseFun, '()', ")")
        )
        )
                 )
      ED_est[[i]] <- ED( model.out[[i]], respLev =  EDs,  ...)
      
    }
    
    
  } else{
    for(i in seq_along(all.curves)){
      
      sel_data <- my_data_a[ my_data_a$CurveName == all.curves[i], ]
      
      model.out[[i]]<- eval(
        parse( text =
                 paste("drm(value ~ dose_uM, type = \"continuous\", 
                       lowerl = lowerl, upperl = upperl, separate = TRUE,
                       data = sel_data, fct =", paste0(doseFun, '()', ")")
                 )
        )
      )
      ED_est[[i]] <- ED( model.out[[i]], respLev =  EDs,  ...)
      
    }
  }
    output = list(model.out = model.out, ED_est = ED_est, CurveName = all.curves)
    
  } else {
    output = list(model.out = model.out, CurveName = all.curves)
  }
  
  
  
  if(finalAnalysis & !is.null(startpars)){
    # starter is in the form bb..cdee..
    if(!is.null(startpars)){
    starter <- c(sapply(startpars, "[[", 1),
                 mean(sapply(startpars, "[[", 2)), mean(sapply(startpars, "[[", 3)),
                 sapply(startpars, "[[", 4))
    
    
    model.out <- eval(
      parse( text = 
               paste("drm(value ~ dose_uM, 
                     type = \"continuous\", data = my_data_a, 
                     curveid = CurveName, start = starter,
                     separate = TRUE, 
                     lowerl = lowerl, upperl = upperl,
                     pmodels = data.frame(CurveName, 1, 1, CurveName),
                     fct =", paste0(doseFun, '()', ")")
      )
      )
      )
    
    ED_est <- ED(model.out, respLev = EDs)
    output = list(model.out = model.out, ED_est = ED_est, curveNames = all.curves)
  } else{
  } 
    model.out <- eval(
      parse( text = 
               paste("drm(value ~ dose_uM, 
                     type = \"continuous\", data = my_data_a, 
                     curveid = CurveName,
                     separate = TRUE, 
                     lowerl = lowerl, upperl = upperl,
                     pmodels = data.frame(CurveName, 1, 1, CurveName),
                     fct =", paste0(doseFun, '()', ")")
      )
      )
      )
    
    ED_est <- ED(model.out, respLev = EDs)
    output = list(model.out = model.out, ED_est = ED_est, curveNames = all.curves)
    
  }
  
  
  if(calibrated | finalAnalysis){
  
    result <-do.call('rbind', output$ED_est)
    nd <- nrow(result) / length(output$CurveName) 
    rownames(result) <- paste(rep(output$CurveName, each = nd),
                              rownames(result), sep = '_')
    
    result <- data.frame(result)
    result$respL <- sapply(
      strsplit(x = rownames(result), split = "_"),
      "[[",2
    )
    result$treatment <- sapply(
      strsplit(x = rownames(result), split = "_"),
      "[[",1
    )
    
    
    result_ED <- dcast(treatment ~ respL, data = result, value.var = 'Estimate')
    result_SD <- dcast(treatment ~ respL, data = result, value.var = 'Std..Error')
    
      output$result_ED <- result_ED
      output$result_SD <- result_SD
  }
  
  
  output
  
  } 