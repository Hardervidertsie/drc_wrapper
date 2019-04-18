### steven explanation splines

################################################################
#install.packages("minpack.lm")
require(minpack.lm)

## get data in similar format

indir <- "J:/Workgroups/FWN/LACDR/TOX/Data Wanda/eutox/case study 4 mitochondria/14_paper/figure_MMP/"
outdir <- "J:/Workgroups/FWN/LACDR/TOX/Data Wanda/eutox/case study 4 mitochondria/14_paper/figure_MMP"
# swetox
swetox <- read.delim(paste0(indir, 'swetox_MMP_longformat.txt'))
swetox_sub <- subset(swetox, select=c("treatment","log_M","replID","cell_line","percentage_of_control"))


##### todo = perform line fitting
# per treatment

head(swetox_sub)

##############HOE MAAK IK HIER EEN LOOP VAN ZODAT DE FUNCTIE VOOR ALLE TREATMENTS WORDT UITGEVOERD EN OPGESLAGEN (zoals bij de drc-wrapper)


ncurves <- length(unique(swetox_sub$treatment))

parlist = NULL
parlist[[ncurves+1]] =list()
unique(swetox_sub$treatment)
skipfit <- c("Capsaicin")

# pas eventueel startwaarden aan als geen fit mogelijk

parlist[[9]] <- list(Bottom = 50, Top = 150,
                     LogEC50= -5.5, HillSlope= -11.5)
parlist[[22]] <- list(Bottom = 50, Top = 140,
                     LogEC50= -6, HillSlope= -5)

curve.nls2= alist()
curve.names = alist()
for(i in seq_along(unique(swetox_sub$treatment))) {

  tstdata = na.omit(subset(swetox_sub, treatment %in% unique(swetox_sub$treatment)[i] ))
  if(is.null(parlist[[i]])){
    parlist[[i]] <- list(Bottom = min(tstdata$percentage_of_control), Top = max(tstdata$percentage_of_control),
                       LogEC50=median(tstdata$log_M), HillSlope= -1)
  }
   


if(!any(skipfit %in% unique(swetox_sub$treatment)[i])  ) {
print(paste0( "attempting to fit: ", unique(swetox_sub$treatment)[i], ", curve number ", i ))
curve.nls2[[i]] <- nlsLM(formula = percentage_of_control ~ Bottom + (Top - Bottom)/(1+10^((LogEC50-log_M)*HillSlope)),
                    data = tstdata, start = parlist[[i]], control = nls.lm.control(maxiter = 200))

curve.names[[i]] <- unique(swetox_sub$treatment)[i]
} else{
  curve.nls2[[i]] <- NA
  curve.names[[i]] <- unique(swetox_sub$treatment)[i]
}
}

# check curve manually : 

plot(percentage_of_control ~ log_M, data = tstdata,main = unique(tstdata$treatment))

# maak tabel met resultaten



summary(curve.nls2[[1]])$parameters

names(curve.nls2)  <- unlist(curve.names)

extract_pars <- lapply(curve.nls2, function(x) summary(x)['parameters'])
extract_pars <- sapply(extract_pars, "[[", 1)

results <- data.frame( Estimate = as.numeric(sapply(extract_pars, "[", 1)),
                       Std_error = as.numeric(sapply(extract_pars, "[", 2)),
                       t_value = as.numeric(sapply(extract_pars, "[", 3)),
                       Pr = as.numeric(sapply(extract_pars, "[", 4)),
                       did_converge = unlist( sapply(sapply(curve.nls2, "[", 'convInfo'), "[", "isConv")))

rownames(results) <- unlist(curve.names)

write.table(results,   file = paste0(outdir, "/resultstable.txt"), sep = "\t", row.names = TRUE, col.names = NA)


###########################
# plot for check
pdf(paste0(outdir, "/plot_fitting.pdf"), width = 24, height = 24)
par(mfrow = c(6,4))
op <- par('mar')
for(i in seq_along(curve.nls2)){
  
  tstdata = na.omit(subset(swetox_sub, treatment %in% unique(swetox_sub$treatment)[i] ))
    par(mar = rep(2, 4))
  
print(plot(percentage_of_control ~ log_M, data = tstdata,main = unique(tstdata$treatment)))

if(!is.na(curve.nls2[[i]])){
print(lines(seq(from = min(tstdata$log_M), to = max(tstdata$log_M), length.out=100), 
      predict(curve.nls2[[i]], 
              newdata = data.frame(log_M = seq(from =  min(tstdata$log_M), to = max(tstdata$log_M), length.out=100))), col = 'red')
)
  }
}


dev.off()
par(mar = op)

# Error in nlsModel(formula, mf, start, wts) :singular gradient matrix at initial parameter estimates
# 3, 9, 22
# Warning message:In nls.lm(par = start, fn = FCT, jac = jac, control = control, lower = lower,  :lmdif: info = -1. Number of iterations has reached `maxiter' == 50.
# 7, 8, 10, 13, 18, 21 
# Error in `[[.default`(tstdata$treatment, i) : subscript out of bounds
# 18, na 18 wordt alles in 1 graph geplot



##### todo = get top value per compound for normalisation
################## create dataframe with all coef values per treatment


## perform normalisation

## plot all curves in one plot


###################
# VU
VU <- read.delim(paste0(indir, 'VU_MMP_longformat.txt'))
VU_mean <- aggregate(percentage_of_control ~ treatment+conc_M+conc_uM+log_M+replID+cell_line, VU, mean)
VU_sub <- subset(VU_mean, select=c("treatment","log_M","replID","cell_line","percentage_of_control"))

# LU
indir1 <- 'I:/Data/Objective 2 use test platform/WS035-37-39 EuToxRisk/'
LU <- read.delim(paste0(indir1, 'WS035-37-39_Rhodamine_normtodmso_average_23-8-2017.txt'))
LU<-subset(LU,grepl("^24$", LU$timeID))
LU$treatment <- as.character(LU$treatment)
LU$normal_value <- LU$normal_value*100
LU$conc_M <- LU$dose_uM/1000000
LU$log_M <- log10(LU$conc_M)
names(LU)[names(LU)=="normal_value"] <- "percentage_of_control"
names(LU)[names(LU)=="dose_uM"] <- "conc_uM"
LU$cell_line <- "HepG2"
LU_sub <- subset(LU, select=c("treatment","log_M","replID","cell_line","percentage_of_control"))



