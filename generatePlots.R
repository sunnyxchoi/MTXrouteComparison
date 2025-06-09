
########################## LOAD PACKAGES ##########################

library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(RMINC)
library(MRIcrotome)
library(reshape2)
library(glue)
library(tidyverse)

########################## LOAD DATA ##########################

# Use this statement below in all scripts to load everything in this environment
load("/projects/schoi/MTXcomparison/analysis/aug2024analysis/linearModelSex.RData")


########################FDR CORRECTION#########################
###############################################################


### Example of q-value across multiple outcomes
### Need to set selcols and eliminate outer loop OR
### Set colpat patterns as needed for your model

q <- 0.1
# dose + route figures FDR
for (colpat in c("^MTX.._IT..._Pr(>|t|)","^MTX.._IT._Pr(>|t|)","^MTX.._IV_Pr(>|t|)","^age_flag..:route_flagIV_Pr(>|t|)", "^age_flag..:Dose._Pr(>|t|)", "^age_flag..:Dose..._Pr(>|t|)")){
  #Satterthwaite approach
  selcols<-colnames(Mres)[grep(colpat,colnames(Mres))]
  allpvals<-c(Mres[,selcols])
  allqvals<-p.adjust(allpvals,method="fdr")
  qcols<-sub("Pr...t..","qvalue",selcols)
  for (k in 1:length(selcols)){
    newcol<-allqvals[( (k-1)*nrow(Mres)+1 ):( k*nrow(Mres) )]
    if (qcols[k] %in% colnames(Mres)){
      Mres[,qcols[k]]<-newcol
    } else {
      Mres<-cbind(Mres,newcol)
      colnames(Mres)[ncol(Mres)]<-qcols[k]
      Mres_FDR <- Mres
    }
  }
}

Mqval <- NULL  # Initialize Mqval as NULL
for (colpat in c("^MTX.._IT..._qvalue","^MTX.._IT._qvalue","^MTX.._IV_qvalue","^age_flag..:route_flagIV_qvalue")){
  selcols<-colnames(Mres_FDR)[grep(colpat,colnames(Mres_FDR))]
  for (k in selcols){
    Mqval<-cbind(Mqval, Mres_FDR[, k, drop=FALSE])
    Msig<-Mqval<q
  }
}

# Save the matrix Mres_FDR as a CSV file
# write.csv(Mres_FDR, "/projects/schoi/MTXcomparison/analysis/finalScripts/Mres_FDR.csv", row.names = TRUE)
# write.csv(Msig, "/projects/schoi/MTXcomparison/analysis/finalScripts/Msig.csv", row.names = TRUE)

########### PLOTTING ############

### Effect calculations ###

eff_MTX24_IT0.5 <- (Mres_FDR[, "MTX24_IT0.5_Estimate"])/(Mres_FDR[, "age_flag24:Dose0.5_Estimate"])*(Msig[, "MTX24_IT0.5_qvalue"])
eff_MTX42_IT0.5 <- (Mres_FDR[, "MTX42_IT0.5_Estimate"])/(Mres_FDR[, "age_flag42:Dose0.5_Estimate"])*(Msig[, "MTX42_IT0.5_qvalue"])
eff_MTX63_IT0.5 <- (Mres_FDR[, "MTX63_IT0.5_Estimate"])/(Mres_FDR[, "age_flag63:Dose0.5_Estimate"])*(Msig[, "MTX63_IT0.5_qvalue"])
eff_MTX24_IT1 <- (Mres_FDR[, "MTX24_IT1_Estimate"])/(Mres_FDR[, "age_flag24:Dose1_Estimate"])*(Msig[, "MTX24_IT1_qvalue"])
eff_MTX42_IT1 <- (Mres_FDR[, "MTX42_IT1_Estimate"])/(Mres_FDR[, "age_flag42:Dose1_Estimate"])*(Msig[, "MTX42_IT1_qvalue"])
eff_MTX63_IT1 <- (Mres_FDR[, "MTX63_IT1_Estimate"])/(Mres_FDR[, "age_flag63:Dose1_Estimate"])*(Msig[, "MTX63_IT1_qvalue"])
eff_MTX24_IT2.5 <- (Mres_FDR[, "MTX24_IT2.5_Estimate"])/(Mres_FDR[, "age_flag24:Dose2.5_Estimate"])*(Msig[, "MTX24_IT2.5_qvalue"])
eff_MTX42_IT2.5 <- (Mres_FDR[, "MTX42_IT2.5_Estimate"])/(Mres_FDR[, "age_flag42:Dose2.5_Estimate"])*(Msig[, "MTX42_IT2.5_qvalue"])
eff_MTX63_IT2.5 <- (Mres_FDR[, "MTX63_IT2.5_Estimate"])/(Mres_FDR[, "age_flag63:Dose2.5_Estimate"])*(Msig[, "MTX63_IT2.5_qvalue"])
eff_MTX24_IT5 <- (Mres_FDR[, "MTX24_IT5_Estimate"])/(Mres_FDR[, "age_flag24:Dose5_Estimate"])*(Msig[, "MTX24_IT5_qvalue"])
eff_MTX42_IT5 <- (Mres_FDR[, "MTX42_IT5_Estimate"])/(Mres_FDR[, "age_flag42:Dose5_Estimate"])*(Msig[, "MTX42_IT5_qvalue"])
eff_MTX63_IT5 <- (Mres_FDR[, "MTX63_IT5_Estimate"])/(Mres_FDR[, "age_flag63:Dose5_Estimate"])*(Msig[, "MTX63_IT5_qvalue"])

eff_MTX24_IV <- (Mres_FDR[, "MTX24_IV_Estimate"])/(Mres_FDR[, "age_flag24:Dose5_Estimate"] + Mres_FDR[, "age_flag24:route_flagIV_Estimate"])*(Msig[, "MTX24_IV_qvalue"])
eff_MTX42_IV <- (Mres_FDR[, "MTX42_IV_Estimate"])/(Mres_FDR[, "age_flag42:Dose5_Estimate"] + Mres_FDR[, "age_flag42:route_flagIV_Estimate"])*(Msig[, "MTX42_IV_qvalue"])
eff_MTX63_IV <- (Mres_FDR[, "MTX63_IV_Estimate"])/(Mres_FDR[, "age_flag63:Dose5_Estimate"] + Mres_FDR[, "age_flag63:route_flagIV_Estimate"])*(Msig[, "MTX63_IV_qvalue"])


# Looking at smallest structures
# Sort the vector in ascending order
sort_eff_MTX24_IT5 <- sort(eff_MTX24_IT5)
sort_eff_MTX24_IV <- sort(eff_MTX24_IV)


### Plots ###

source('/projects/dfernandes/immune_strains/final_analysis_nov18/anatomy_pruning/tutorial/label_functions.R')

mask <- '/projects/schoi/MTXcomparison/pipelines/hpftransfer/MTX17022024/atlas/nlin_labels/nlin_new.mnc'
mappings <- '/projects/schoi/DSURQEE_60micron_invivo/DSURQEE_R_mapping.csv'
dt.struc.run <- read.csv(mappings)
structnames <- dt.struc.run$Structure

anatVol <- mincArray(mincGetVolume('/projects/schoi/MTXcomparison/pipelines/hpftransfer/MTX17022024/MTX17022024_nlin/MTX17022024-nlin-3.mnc'))


plot_MTX24_IT0.5 <- struct_wise_to_vxl_wise(eff_MTX24_IT0.5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX42_IT0.5 <- struct_wise_to_vxl_wise(eff_MTX42_IT0.5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX63_IT0.5 <- struct_wise_to_vxl_wise(eff_MTX63_IT0.5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX24_IT1 <- struct_wise_to_vxl_wise(eff_MTX24_IT1, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX42_IT1 <- struct_wise_to_vxl_wise(eff_MTX42_IT1, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX63_IT1 <- struct_wise_to_vxl_wise(eff_MTX63_IT1, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX24_IT2.5 <- struct_wise_to_vxl_wise(eff_MTX24_IT2.5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX42_IT2.5 <- struct_wise_to_vxl_wise(eff_MTX42_IT2.5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX63_IT2.5 <- struct_wise_to_vxl_wise(eff_MTX63_IT2.5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX24_IT5 <- struct_wise_to_vxl_wise(eff_MTX24_IT5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX42_IT5 <- struct_wise_to_vxl_wise(eff_MTX42_IT5, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX63_IT5 <- struct_wise_to_vxl_wise(eff_MTX63_IT5, mincGetVolume(mask), structure.names = structnames, defs = mappings)

plot_MTX24_IV <- struct_wise_to_vxl_wise(eff_MTX24_IV, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX42_IV <- struct_wise_to_vxl_wise(eff_MTX42_IV, mincGetVolume(mask), structure.names = structnames, defs = mappings)
plot_MTX63_IV <- struct_wise_to_vxl_wise(eff_MTX63_IV, mincGetVolume(mask), structure.names = structnames, defs = mappings)


########################## SAVE THE ENVIRONMENT ##########################

# Saving the environment
save.image("/projects/schoi/MTXcomparison/analysis/finalScripts/generatePlots.RData")
# Use this statement below in all scripts to load everything in this environment
#load("/projects/schoi/MTXcomparison/analysis/finalScripts/generatePlots.RData")

