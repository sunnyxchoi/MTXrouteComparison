library(tidyverse)
library(multcomp)
library(lme4)
library(lmerTest)
library(ggplot2)
library(RMINC)
library(grid)
library(gridExtra)
library(MRIcrotome)

########################## LOAD DATA ##########################

# Use this statement below in all scripts to load everything in this environment
load("/projects/schoi/MTXcomparison/analysis/finalScripts/generatePlots.RData")



#######################################################################
# utility functions                                                   #
#######################################################################

clamp<-function(x,lower=-Inf,upper=Inf){
  ifelse(x>upper,upper,ifelse(x<lower,lower,x))
}

anatResToMincArray<-function(anatVecToOutput,labelArray,atlasDefs,sigflag=NULL){
  dfAtlasLabels<-read.csv(atlasDefs,stringsAsFactors=FALSE)
  if (any( !(names(anatVecToOutput)%in%dfAtlasLabels$Structure) )){
    print("WARNING (anatResToMincArray): Guessing structure matches based on regex.")
    #ditch all the spaces and punctuation and try matching
    srchpatt<-gsub("ZZZZZ","..*",gsub("\\.\\.","[^a-z]*",gsub("\\.$","",
                                                              gsub("[[:punct:]]",".",gsub("'","ZZZZZ",gsub(" ",".",dfAtlasLabels$Structure))))))
    for (j in 1:length(srchpatt)){
      mtchind<-grep(paste0("^",srchpatt[j],"$"),names(anatVecToOutput))
      dfAtlasLabels$Structure[j]<-names(anatVecToOutput)[mtchind]
    }
  }    
  if (typeof(anatVecToOutput)=="character"){
    Aout<-array("#00000000",dim=dim(labelArray))
  } else {
    Aout<-array(0,dim=dim(labelArray)) #rep(0,length(labelimg))
  }
  if (is.null(sigflag)){
    sigflag<-rep(TRUE,length(anatVecToOutput))
  }
  for (j in 1:length(anatVecToOutput)){
    print(paste0("We entered the for loop! Stage:", j))
    cstruct<-names(anatVecToOutput)[j]
    clabels<-dfAtlasLabels[dfAtlasLabels$Structure==cstruct,c("right.label","left.label")]
    if (sigflag[j]){
      outval<-anatVecToOutput[cstruct]
      Aout[abs(labelArray-as.numeric(clabels[1]))<0.5]<-outval
      Aout[abs(labelArray-as.numeric(clabels[2]))<0.5]<-outval 
    }
  }
  print("Returning Aout")
  return(Aout)
}


genOverlayArray<-function(structurenames,outvals,man_colVol,manRcolVol,statLowVol,statHighVol,transparencyvec,labelVol,labelmappings){
  dftmp<-data.frame(structure=structurenames,outval=outvals)
  dftmp$opaquecolour<-ifelse(dftmp$outval>0.0, 
                             man_colVol[as.integer(clamp( length(man_colVol)*(dftmp$outval-statLowVol)/(statHighVol-statLowVol) ,lower=1,upper=length(man_colVol)))],
                             man_RcolVol[as.integer(clamp( length(man_colVol)*(-dftmp$outval-statLowVol)/(statHighVol-statLowVol) ,lower=1,upper=length(man_colVol)))]
  )
  colvec<-colorspace::adjust_transparency(dftmp$opaquecolour,transparencyvec)
  names(colvec)<-dftmp$structure
  overlayVol<-anatResToMincArray(colvec,labelVol,labelmappings,NULL)
  return(overlayVol)
}


genHorizLegend <- function(clrscale,txtlabel="Value (a.u.)",loval=0,hival=1,
                           textsize=4,topM=0.25,bottomM=-0.7,limtextvposadj=0,limtexthposadj=0,lbltextvposadj=0){
  dummydf<-data.frame(id=rep(1,length(clrscale)),val=1:length(clrscale))
  nchar<-max(length(as.character(loval)),length(as.character(hival)))
  val=seq(1,length(clrscale))
  xpos=val/length(clrscale)
  hlgndplot<-ggplot(dummydf) +
    geom_tile(aes(x = xpos, y=0, height=0.5,fill = val)) +
    scale_y_continuous(expand=c(0,0.0),limits=c(bottomM,topM),breaks=1)+
    scale_x_continuous(limits=c(-limtexthposadj,1.0+limtexthposadj),expand = c(0,0)) +
    scale_fill_gradientn(colors=clrscale) +
    theme_minimal()+
    theme(legend.position="none",panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.border=element_rect(colour="transparent",fill="transparent"),panel.background=element_rect(fill="transparent",colour=NA),
          axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),
          plot.margin=margin(0,2*nchar*textsize,3*textsize,2*nchar*textsize,"mm"))+
    annotate("text", x = c(-limtexthposadj), y = c(-0.25-limtextvposadj),hjust=0,vjust=1,label = c(as.character(loval)),color="black",size=textsize,angle=0)+
    annotate("text", x = c(1.0+limtexthposadj), y=c(-0.25-limtextvposadj),hjust=1,vjust=1,label = c(as.character(hival)),color="black",size=textsize,angle=0)+
    annotate("text", x = c(0.5), y=c(-0.5-lbltextvposadj),hjust=0.5,vjust=1,label=txtlabel,color="black",size=textsize,angle=0)
  return(hlgndplot)
}


source('/projects/dfernandes/immune_strains/final_analysis_nov18/anatomy_pruning/tutorial/label_functions.R')

crop_func = function(x) {
  nsx = names(attributes(x)) 
  sx = attributes(x)[c('class','filename','likeVolume')]
  x = x[,40:320,]
  attributes(x) = c(sx,attributes(x),sizes=list(rev(dim(x))))
  x}

#############################################################
# Inputs for neuroanatomical maps
#############################################################

#full-width layout of image panels
nlin3avg<-'/projects/schoi/MTXcomparison/pipelines/hpftransfer/MTX17022024/MTX17022024_nlin/MTX17022024-nlin-3.mnc'
nlin3labels<-'/projects/schoi/MTXcomparison/pipelines/hpftransfer/MTX17022024/atlas/nlin_labels/nlin_new.mnc'
labelmappings<-'/projects/schoi/DSURQEE_60micron_invivo/DSURQEE_R_mapping.csv'

anatVol <- mincArray(mincGetVolume(nlin3avg))
labelVol <- mincArray(mincGetVolume(nlin3labels))

statLowVol <- 1.0
statHighVol <- 15.0

man_colVol<-colorRampPalette(c("#881802","#c05004","#d08008","#ffa060"))(128)
man_RcolVol<-colorRampPalette(c("#023858","#0470a0","#08a0b0","#60c0c0"))(128)


########################## IMPORT ATLAS STRUCTURE ##########################

dt.struc <- read.csv('/projects/schoi/DSURQEE_60micron_invivo/DSURQEE_R_mapping.csv')
#dt.struc$Structure <- str_replace_all(dt.struc$Structure, "[^[:alnum:]]", " ") %>% str_squish() %>% str_replace_all(., ' ', '_')
#dt.struc$Structure <- make.unique(dt.struc$Structure)
structure_list <- dt.struc$Structure



#image mappings output
for (croute in c("IT")){
  Lgrobs=list()
  for (cdose in c("0.5", "5")){
    for (ctpt in c("24", "42", "63")){
      ccolref<-paste0("MTX",ctpt,"_", croute, cdose)
      if (!paste0(ccolref,"_Pr(>|t|)")%in%colnames(Mres)){ 
        ccolref<-paste0("MTX",ctpt,"_", croute, cdose)
      }
      Valpha<-ifelse( Mres[,paste0(ccolref,"_Pr(>|t|)")]>0.05 ,0.0,
                      ifelse( (Mres[,paste0(ccolref,"_Pr(>|t|)")]<0.05)&(Mres[,paste0(ccolref,"_qvalue")]>=0.1) ,0.3,1.0))[structure_list]
      #names(Valpha)<-names(structure_list)
      overlayVol<-genOverlayArray(names(Valpha),
                                  (100*Mres[,paste0(ccolref,"_Estimate")]/Mres[,paste0("age_flag",ctpt,":Dose", cdose, "_Estimate")])[structure_list],
                                  man_colVol,manRcolVol,statLowVol,statHighVol,Valpha,labelVol,labelmappings)
      print("We passed overlay!")
      lGrpST <- sliceSeries(begin=105, end = 105, ncol=1, nrow=1, dimension=3) %>% 
        anatomy(crop_func(anatVol), low=700, high=1400) %>% 
        overlay(crop_func(overlayVol)) %>% 
        grobify()
      Lgrobs[[ccolref]]<-lGrpST
    }
  }
  Nwht<-round(2*(statLowVol/statHighVol)*length(man_RcolVol)/(1.0-statLowVol/statHighVol))
  lgndtxt<-c("0.5"="0.5 mg/kg (%)",
             "5.0"="5.0 mg/kg (%)")[cdose]
  llgnd<-genHorizLegend( c(rev(man_RcolVol),rep("#FFFFFF",Nwht),man_colVol),lgndtxt,-statHighVol,statHighVol,
                         textsize=3.7,topM=0.6,bottomM=-0.8,limtextvposadj=0.1,lbltextvposadj=-0.1)
  colHead1<-textGrob("IT MTX 0.5 mg/kg",gp=gpar(fontsize=10))
  colHead2<-textGrob("IT MTX 5.0 mg/kg",gp=gpar(fontsize=10))
  png(file = paste0("/projects/schoi/MTXcomparison/figures/", "fig5", "_",croute, ".png"), width = 2.5, height = 5, units = "in", res = 300)
  grid.arrange(grobs=list(colHead1,colHead2,
                          Lgrobs[[grep("MTX24_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT5",names(Lgrobs))]],
                          Lgrobs[[grep("MTX42_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX42_IT5",names(Lgrobs))]],
                          Lgrobs[[grep("MTX63_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX63_IT5",names(Lgrobs))]],
                          llgnd),
               widths=c(1,1),
               heights=c(0.2,1,1,1,0.05,0.8),
               ncol=2,nrow=6,
               layout_matrix = rbind(
                 c(1, 2),  # Column headers
                 c(3, 4),  # Plots for P24
                 c(5, 6), # Plots for P42
                 c(7, 8), # Plots for P63
                 c(NA, NA),  # Empty row for spacing (optional, can be removed)
                 c(9, 9))   # Legend spanning all columns
  )
  dev.off()
}


# Run this to display IT neuroanatomical maps on Rstudio plots
grid.arrange(grobs=list(colHead1,colHead2,
                        Lgrobs[[grep("MTX24_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT5",names(Lgrobs))]],
                        Lgrobs[[grep("MTX42_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX42_IT5",names(Lgrobs))]],
                        Lgrobs[[grep("MTX63_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX63_IT5",names(Lgrobs))]],
                        llgnd),
             widths=c(1,1),
             heights=c(0.2,1,1,1,0.05,0.8),
             ncol=2,nrow=6,
             layout_matrix = rbind(
               c(1, 2),  # Column headers
               c(3, 4),  # Plots for P24
               c(5, 6), # Plots for P42
               c(7, 8), # Plots for P63
               c(NA, NA),  # Empty row for spacing (optional, can be removed)
               c(9, 9))   # Legend spanning all columns
)


#image mappings output
for (croute in c("IV")){
  Lgrobs=list()
  for (ctpt in c("24", "42", "63")){
    ccolref<-paste0("MTX",ctpt,"_", croute)
    if (!paste0(ccolref,"_Pr(>|t|)")%in%colnames(Mres)){ 
      ccolref<-paste0("MTX",ctpt,"_", croute)
    }
    Valpha<-ifelse( Mres[,paste0(ccolref,"_Pr(>|t|)")]>0.05 ,0.0,
                    ifelse( (Mres[,paste0(ccolref,"_Pr(>|t|)")]<0.05)&(Mres[,paste0(ccolref,"_qvalue")]>=0.1) ,0.3,1.0))[structure_list]
    #names(Valpha)<-names(structure_list)
    overlayVol<-genOverlayArray(names(Valpha),
                                (100*Mres[,paste0(ccolref,"_Estimate")]/Mres[,paste0("age_flag",ctpt,":Dose", cdose, "_Estimate")])[structure_list],
                                man_colVol,manRcolVol,statLowVol,statHighVol,Valpha,labelVol,labelmappings)
    print("We passed overlay!")
    lGrpST <- sliceSeries(begin=105, end = 105, ncol=1, nrow=1, dimension=3) %>% 
      anatomy(crop_func(anatVol), low=700, high=1400) %>% 
      overlay(crop_func(overlayVol)) %>% 
      grobify()
    Lgrobs[[ccolref]]<-lGrpST
  }
  Nwht<-round(2*(statLowVol/statHighVol)*length(man_RcolVol)/(1.0-statLowVol/statHighVol))
  lgndtxt<-c("IV"="Volume Change (%)")[croute]
  llgnd<-genHorizLegend( c(rev(man_RcolVol),rep("#FFFFFF",Nwht),man_colVol),lgndtxt,-statHighVol,statHighVol,
                         textsize=3.7,topM=0.6,bottomM=-0.8,limtextvposadj=0.1,lbltextvposadj=-0.1)
  colHead1<-textGrob("IV MTX 5.0 mg/kg",gp=gpar(fontsize=10))
  png(file = paste0("/projects/schoi/MTXcomparison/figures/", "fig5", "_",croute, ".png"), width = 1.25, height = 5, units = "in", res = 300)
  grid.arrange(grobs=list(colHead1,
                          Lgrobs[[grep("MTX24_IV",names(Lgrobs))]],
                          Lgrobs[[grep("MTX42_IV",names(Lgrobs))]],
                          Lgrobs[[grep("MTX63_IV",names(Lgrobs))]],
                          llgnd),
               widths=c(1),
               heights = c(0.2, 1, 1, 1, 0.05, 0.8),
               # heights=c(0.2,1,1,1,0.05,0.8),
               ncol=1,nrow=6,
               layout_matrix = rbind(
                 c(1),  # Column headers
                 c(2),  # Plots for P24
                 c(3), # Plots for P42
                 c(4), # Plots for P63
                 c(NA),  # Empty row for spacing (optional, can be removed)
                 c(5))   # Legend spanning all columns
  )
  dev.off()
}


# Run this to display IV neuroanatomical maps on Rstudio plots
grid.arrange(grobs=list(colHead1,
                        Lgrobs[[grep("MTX24_IV",names(Lgrobs))]],
                        Lgrobs[[grep("MTX42_IV",names(Lgrobs))]],
                        Lgrobs[[grep("MTX63_IV",names(Lgrobs))]],
                        llgnd),
             widths=c(1),
             heights = c(0.2, 1, 1, 1, 0.05, 0.8),
             # heights=c(0.2,1,1,1,0.05,0.8),
             ncol=1,nrow=6,
             layout_matrix = rbind(
               c(1),  # Column headers
               c(2),  # Plots for P24
               c(3), # Plots for P42
               c(4), # Plots for P63
               c(NA),  # Empty row for spacing (optional, can be removed)
               c(5))   # Legend spanning all columns
)


########################## FIG 5B: DENSITY PLOTS ##########################

#P24
change_MTX24_IT0.5 <- (Mres[, "MTX24_IT0.5_Estimate"])/(Mres[, "age_flag24:Dose0.5_Estimate"]) *100
change_MTX24_IT5 <- (Mres[, "MTX24_IT5_Estimate"])/(Mres[, "age_flag24:Dose5_Estimate"]) *100
change_MTX24_IV <- (Mres[, "MTX24_IV_Estimate"])/(Mres_FDR[, "age_flag24:Dose5_Estimate"] + Mres_FDR[, "age_flag24:route_flagIV_Estimate"]) *100

#P42
change_MTX42_IT0.5 <- (Mres[, "MTX42_IT0.5_Estimate"])/(Mres[, "age_flag42:Dose0.5_Estimate"]) *100
change_MTX42_IT5 <- (Mres[, "MTX42_IT5_Estimate"])/(Mres[, "age_flag42:Dose5_Estimate"]) *100
change_MTX42_IV <- (Mres[, "MTX42_IV_Estimate"])/(Mres_FDR[, "age_flag42:Dose5_Estimate"] + Mres_FDR[, "age_flag42:route_flagIV_Estimate"]) *100

#P63
change_MTX63_IT0.5 <- (Mres[, "MTX63_IT0.5_Estimate"])/(Mres[, "age_flag63:Dose0.5_Estimate"]) *100
change_MTX63_IT5 <- (Mres[, "MTX63_IT5_Estimate"])/(Mres[, "age_flag63:Dose5_Estimate"]) *100
change_MTX63_IV <- (Mres[, "MTX63_IV_Estimate"])/(Mres_FDR[, "age_flag63:Dose5_Estimate"] + Mres_FDR[, "age_flag63:route_flagIV_Estimate"]) *100


# Density lines with all doses on same plot

# P24

# Set the margins
par(mar = c(5, 5, 4, 4))

# Set the background to white (mimicking theme_bw)
par(bg = "white")

# Plot the first density line
plot(density(change_MTX24_IT0.5), xlim = c(-25, 15), ylim = c(0, 0.24), 
     main = "P24", 
     xlab = "% Volume Change", ylab = "Density", col = "blue", lwd = 5, 
     cex.axis = 1.2, cex.lab = 1.3, font.lab = 2, cex.main = 1.6)  # Modify axis and label sizes

# Add the other density lines
lines(density(change_MTX24_IT5), col = "green", lwd = 5)
lines(density(change_MTX24_IV), col = "orange", lwd = 5)

# Optional: Adjust gridlines to make the plot cleaner
grid(lwd = 0.5, col = "gray60")

# Legend appears on lines so just copy and paste onto figure at end
graphics::legend("right", legend = c("IT 0.5 mg/kg", "IT 5.0 mg/kg", "IV 5.0 mg/kg"), 
                 col = c("blue", "green", "orange"), lwd =5, xpd=FALSE)


# P42
# Set the margins
par(mar = c(5, 5, 4, 4))

# Set the background to white (mimicking theme_bw)
par(bg = "white")

# Plot the first density line
plot(density(change_MTX42_IT0.5), xlim = c(-25, 15), ylim = c(0, 0.24), 
     main = "P42", 
     xlab = "% Volume Change", ylab = "Density", col = "blue", lwd = 5, 
     cex.axis = 1.2, cex.lab = 1.3, font.lab = 2, cex.main = 1.6)  # Modify axis and label sizes

# Add the other density lines
lines(density(change_MTX42_IT5), col = "green", lwd = 5)
lines(density(change_MTX42_IV), col = "orange", lwd = 5)

# Optional: Adjust gridlines to make the plot cleaner
grid(lwd = 0.5, col = "gray60")

# Legend appears on lines so just copy and paste onto figure at end
graphics::legend("right", legend = c("IT 0.5 mg/kg", "IT 5.0 mg/kg", "IV 5.0 mg/kg"), 
                 col = c("blue", "green", "orange"), lwd =5, xpd=FALSE)



# P63
# Set the margins
par(mar = c(5, 5, 4, 4))

# Set the background to white (mimicking theme_bw)
par(bg = "white")

# Plot the first density line
plot(density(change_MTX63_IT0.5), xlim = c(-25, 15), ylim = c(0, 0.24), 
     main = "P63", 
     xlab = "% Volume Change", ylab = "Density", col = "blue", lwd = 5, 
     cex.axis = 1.2, cex.lab = 1.3, font.lab = 2, cex.main = 1.6)  # Modify axis and label sizes

# Add the other density lines
lines(density(change_MTX63_IT5), col = "green", lwd = 5)
lines(density(change_MTX63_IV), col = "orange", lwd = 5)

# Optional: Adjust gridlines to make the plot cleaner
grid(lwd = 0.5, col = "gray60")

# Legend appears on lines so just copy and paste onto figure at end
graphics::legend("right", legend = c("IT 0.5 mg/kg", "IT 5.0 mg/kg", "IV 5.0 mg/kg"), 
                 col = c("blue", "green", "orange"), lwd =5, xpd=FALSE)




### Perform two-sided Mann-Whitney U tests for each comparison at P24 ###

# IV 5.0 mg/kg vs IT 5.0 mg/kg at P24
mw_IV5_vs_IT5_two_sided_P24 <- wilcox.test(change_MTX24_IV, change_MTX24_IT5, alternative = "two.sided")

# IV 5.0 mg/kg vs IT 0.5 mg/kg at P24
mw_IV5_vs_IT0.5_two_sided_P24 <- wilcox.test(change_MTX24_IV, change_MTX24_IT0.5, alternative = "two.sided")

# Print the p-values for the comparisons at P24
print(mw_IV5_vs_IT5_two_sided_P24$p.value)
print(mw_IV5_vs_IT0.5_two_sided_P24$p.value)


# Perform two-sided Mann-Whitney U tests for each comparison at P42

# IV 5.0 mg/kg vs IT 5.0 mg/kg at P42
mw_IV5_vs_IT5_two_sided_P42 <- wilcox.test(change_MTX42_IV, change_MTX42_IT5, alternative = "two.sided")

# IV 5.0 mg/kg vs IT 0.5 mg/kg at P42
mw_IV5_vs_IT0.5_two_sided_P42 <- wilcox.test(change_MTX42_IV, change_MTX42_IT0.5, alternative = "two.sided")

# Print the p-values for the comparisons at P42
print(mw_IV5_vs_IT5_two_sided_P42$p.value)
print(mw_IV5_vs_IT0.5_two_sided_P42$p.value)


# Perform two-sided Mann-Whitney U tests for each comparison at P63

# IV 5.0 mg/kg vs IT 5.0 mg/kg at P63
mw_IV5_vs_IT5_two_sided_P63 <- wilcox.test(change_MTX63_IV, change_MTX63_IT5, alternative = "two.sided")

# IV 5.0 mg/kg vs IT 0.5 mg/kg at P63
mw_IV5_vs_IT0.5_two_sided_P63 <- wilcox.test(change_MTX63_IV, change_MTX63_IT0.5, alternative = "two.sided")

# Print the p-values for the comparisons at P63
print(mw_IV5_vs_IT5_two_sided_P63$p.value)
print(mw_IV5_vs_IT0.5_two_sided_P63$p.value)



########################## FIG 5C: SEX EFFECTS NEUROANATOMICAL MAPS ##########################

# Use this statement below in all scripts to load everything in this environment
load("/projects/schoi/MTXcomparison/analysis/finalScripts/generateSexPlots.RData")


statLowVol <- 1.0
statHighVol <- 15.0

# Pink gradient
man_colVol <- colorRampPalette(c("#C71585", "#FF1493", "#FF69B4", "#FFB6C1"))(128)
# Green gradient
man_RcolVol <- colorRampPalette(c("#006400", "#32CD32", "#7FFF00", "#ADFF2F"))(128)


#image mappings output
for (croute in c("IT")){
  Lgrobs=list()
  for (cdose in c("0.5", "5")){
    for (ctpt in c("24", "42", "63")){
      ccolref<-paste0("MTX",ctpt,"_", croute, cdose, "_sex")
      if (!paste0(ccolref,"_Pr(>|t|)")%in%colnames(Mres_sexFDR)){ 
        ccolref<-paste0("MTX",ctpt,"_", croute, cdose, "_sex")
      }
      Valpha<-ifelse( Mres_sexFDR[,paste0(ccolref,"_Pr(>|t|)")]>0.05 ,0.0,
                      ifelse( (Mres_sexFDR[,paste0(ccolref,"_Pr(>|t|)")]<0.05)&(Mres_sexFDR[,paste0(ccolref,"_qvalue")]>=0.1) ,0.3,1.0))[structure_list]
      #names(Valpha)<-names(structure_list)
      overlayVol<-genOverlayArray(names(Valpha),
                                  (100*Mres_sexFDR[,paste0(ccolref,"_Estimate")]/Mres_sexFDR[,paste0("age_flag",ctpt,":Dose", cdose, "_Estimate")])[structure_list],
                                  man_colVol,manRcolVol,statLowVol,statHighVol,Valpha,labelVol,labelmappings)
      print("We passed overlay!")
      lGrpST <- sliceSeries(begin=105, end = 105, ncol=1, nrow=1, dimension=3) %>% 
        anatomy(crop_func(anatVol), low=700, high=1400) %>% 
        overlay(crop_func(overlayVol)) %>% 
        grobify()
      Lgrobs[[ccolref]]<-lGrpST
    }
  }
  Nwht<-round(2*(statLowVol/statHighVol)*length(man_RcolVol)/(1.0-statLowVol/statHighVol))
  lgndtxt<-c("0.5"="0.5 mg/kg (%)",
             "5.0"="5.0 mg/kg (%)")[cdose]
  llgnd<-genHorizLegend( c(rev(man_RcolVol),rep("#FFFFFF",Nwht),man_colVol),lgndtxt,-statHighVol,statHighVol,
                         textsize=3.7,topM=0.6,bottomM=-0.8,limtextvposadj=0.1,lbltextvposadj=-0.1)
  colHead1<-textGrob("IT MTX 0.5 mg/kg",gp=gpar(fontsize=10))
  colHead2<-textGrob("IT MTX 5.0 mg/kg",gp=gpar(fontsize=10))
  png(file = paste0("/projects/schoi/MTXcomparison/figures/", "fig5", "_",croute, "_sex_", ".png"), width = 2.5, height = 5, units = "in", res = 300)
  grid.arrange(grobs=list(colHead1,colHead2,
                          Lgrobs[[grep("MTX24_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT5",names(Lgrobs))]],
                          Lgrobs[[grep("MTX42_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX42_IT5",names(Lgrobs))]],
                          Lgrobs[[grep("MTX63_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX63_IT5",names(Lgrobs))]],
                          llgnd),
               widths=c(1,1),
               heights=c(0.2,1,1,1,0.05,0.8),
               ncol=2,nrow=6,
               layout_matrix = rbind(
                 c(1, 2),  # Column headers
                 c(3, 4),  # Plots for P24
                 c(5, 6), # Plots for P42
                 c(7, 8), # Plots for P63
                 c(NA, NA),  # Empty row for spacing (optional, can be removed)
                 c(9, 9))   # Legend spanning all columns
  )
  dev.off()
}


# Run this for IT sex effects neuroanatomical maps in Rstudio plots
grid.arrange(grobs=list(colHead1,colHead2,
                        Lgrobs[[grep("MTX24_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT5",names(Lgrobs))]],
                        Lgrobs[[grep("MTX42_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX42_IT5",names(Lgrobs))]],
                        Lgrobs[[grep("MTX63_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX63_IT5",names(Lgrobs))]],
                        llgnd),
             widths=c(1,1),
             heights=c(0.2,1,1,1,0.05,0.8),
             ncol=2,nrow=6,
             layout_matrix = rbind(
               c(1, 2),  # Column headers
               c(3, 4),  # Plots for P24
               c(5, 6), # Plots for P42
               c(7, 8), # Plots for P63
               c(NA, NA),  # Empty row for spacing (optional, can be removed)
               c(9, 9))   # Legend spanning all columns
)


#image mappings output
for (croute in c("IV")){
  Lgrobs=list()
  for (ctpt in c("24", "42", "63")){
    ccolref<-paste0("MTX",ctpt,"_", croute, "_sex")
    if (!paste0(ccolref,"_Pr(>|t|)")%in%colnames(Mres_sexFDR)){ 
      ccolref<-paste0("MTX",ctpt,"_", croute, "_sex")
    }
    Valpha<-ifelse( Mres_sexFDR[,paste0(ccolref,"_Pr(>|t|)")]>0.05 ,0.0,
                    ifelse( (Mres_sexFDR[,paste0(ccolref,"_Pr(>|t|)")]<0.05)&(Mres_sexFDR[,paste0(ccolref,"_qvalue")]>=0.1) ,0.3,1.0))[structure_list]
    #names(Valpha)<-names(structure_list)
    overlayVol<-genOverlayArray(names(Valpha),
                                (100*Mres_sexFDR[,paste0(ccolref,"_Estimate")]/Mres_sexFDR[,paste0("age_flag",ctpt,":Dose", cdose, "_Estimate")])[structure_list],
                                man_colVol,manRcolVol,statLowVol,statHighVol,Valpha,labelVol,labelmappings)
    print("We passed overlay!")
    lGrpST <- sliceSeries(begin=105, end = 105, ncol=1, nrow=1, dimension=3) %>% 
      anatomy(crop_func(anatVol), low=700, high=1400) %>% 
      overlay(crop_func(overlayVol)) %>% 
      grobify()
    Lgrobs[[ccolref]]<-lGrpST
  }
  Nwht<-round(2*(statLowVol/statHighVol)*length(man_RcolVol)/(1.0-statLowVol/statHighVol))
  lgndtxt<-c("IV"="Volume Change (%)")[croute]
  llgnd<-genHorizLegend( c(rev(man_RcolVol),rep("#FFFFFF",Nwht),man_colVol),lgndtxt,-statHighVol,statHighVol,
                         textsize=3.7,topM=0.6,bottomM=-0.8,limtextvposadj=0.1,lbltextvposadj=-0.1)
  colHead1<-textGrob("IV MTX 5.0 mg/kg",gp=gpar(fontsize=10))
  png(file = paste0("/projects/schoi/MTXcomparison/figures/", "fig5", "_",croute, "_sex_", ".png"), width = 1.25, height = 5, units = "in", res = 300)
  grid.arrange(grobs=list(colHead1,
                          Lgrobs[[grep("MTX24_IV",names(Lgrobs))]],
                          Lgrobs[[grep("MTX42_IV",names(Lgrobs))]],
                          Lgrobs[[grep("MTX63_IV",names(Lgrobs))]],
                          llgnd),
               widths=c(1),
               heights = c(0.2, 1, 1, 1, 0.05, 0.8),
               # heights=c(0.2,1,1,1,0.05,0.8),
               ncol=1,nrow=6,
               layout_matrix = rbind(
                 c(1),  # Column headers
                 c(2),  # Plots for P24
                 c(3), # Plots for P42
                 c(4), # Plots for P63
                 c(NA),  # Empty row for spacing (optional, can be removed)
                 c(5))   # Legend spanning all columns
  )
  dev.off()
}


# Run this for IV sex effects neuroanatomical maps in Rstudio plots
grid.arrange(grobs=list(colHead1,
                        Lgrobs[[grep("MTX24_IV",names(Lgrobs))]],
                        Lgrobs[[grep("MTX42_IV",names(Lgrobs))]],
                        Lgrobs[[grep("MTX63_IV",names(Lgrobs))]],
                        llgnd),
             widths=c(1),
             heights = c(0.2, 1, 1, 1, 0.05, 0.8),
             # heights=c(0.2,1,1,1,0.05,0.8),
             ncol=1,nrow=6,
             layout_matrix = rbind(
               c(1),  # Column headers
               c(2),  # Plots for P24
               c(3), # Plots for P42
               c(4), # Plots for P63
               c(NA),  # Empty row for spacing (optional, can be removed)
               c(5))   # Legend spanning all columns
)


# Use below for sex gradient legend

sex_fig.24routehorizontal<- sliceSeries(nrow = 1, ncol = 1, begin=105, end =105,dimension = 3)%>%
  anatomy(anatVol, low = 700, high = 1400) %>% 
  overlay(mincArray(sex_plot_MTX24_IV*100), low=1, high=15, symmetric =T, col = man_colVol, man_RcolVol) %>%
  addtitle('IV (5.0 mg/kg)') %>% 
  MRIcrotome::legend('% Volume Change') %>% 
  contourSliceIndicator(anatVol, c(700, 1400)) %>% draw()



# Figure 5D

########################## LOAD PACKAGES ##########################

library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(RMINC)
library(MRIcrotome)
library(reshape2)
library(glue)
library(grid)
library(gridExtra)

########################## LOAD DATA ##########################

# Use this statement below in all scripts to load everything in this environment
load("/projects/schoi/MTXcomparison/analysis/MTX17022024/med_MTX.RData")

df_dose <- med_df_MTX %>% filter(Dose == 5)
df_dose$Group <- paste(df_dose$Treatment, df_dose$Sex, sep = "_")


MTX24_IT0.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "M")
MTX24_IT0.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "F")
MTX24_IT1_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "M")
MTX24_IT1_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "F")
MTX24_IT2.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "M")
MTX24_IT2.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "F")
MTX24_IT5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "M")
MTX24_IT5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "F")

MTX42_IT0.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "M")
MTX42_IT0.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "F")
MTX42_IT1_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "M")
MTX42_IT1_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "F")
MTX42_IT2.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "M")
MTX42_IT2.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "F")
MTX42_IT5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "M")
MTX42_IT5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "F")

MTX63_IT0.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "M")
MTX63_IT0.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "F")
MTX63_IT1_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "M")
MTX63_IT1_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "F")
MTX63_IT2.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "M")
MTX63_IT2.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "F")
MTX63_IT5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "M")
MTX63_IT5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "F")


########################## SETTING UP IV COEFFICIENTS ##########################
################################################################################

MTX24_IV_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "M")
MTX24_IV_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "F")
MTX42_IV_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "M")
MTX42_IV_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "F")
MTX63_IV_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "M")
MTX63_IV_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "F")


########################## GENERAL FLAGS ##########################

age_flag <- as.character(med_df_MTX$Age)
sex_flag <- med_df_MTX$Sex
dose_flag <- as.character(med_df_MTX$Dose)
study_id <- med_df_MTX$Study_ID
pregnancy_cage <- med_df_MTX$Pregnancy_cage
med_df_MTX$Group <- paste(med_df_MTX$Route, med_df_MTX$Dose, sep = "_")
group <- med_df_MTX$Group
route_flag <- med_df_MTX$Route



med_df_MTX_ranef <- med_df_MTX
for (j in 1:length(structure_list)){
  cstruct <- structure_list[j]
  print(paste(as.character(j), cstruct, sep=": "))
  first_form <- glue("`{cstruct}` ~ -1 + age_flag + age_flag:route_flag + age_flag:sex_flag +
                      MTX24_IT0.5_M + MTX24_IT0.5_F + 
                      MTX24_IT1_M + MTX24_IT1_F + 
                      MTX24_IT2.5_M + MTX24_IT2.5_F + 
                      MTX24_IT5_M + MTX24_IT5_F + 
                      MTX42_IT0.5_M + MTX42_IT0.5_F + 
                      MTX42_IT1_M + MTX42_IT1_F + 
                      MTX42_IT2.5_M + MTX42_IT2.5_F + 
                      MTX42_IT5_M + MTX42_IT5_F + 
                      MTX63_IT0.5_M + MTX63_IT0.5_F + 
                      MTX63_IT1_M + MTX63_IT1_F + 
                      MTX63_IT2.5_M + MTX63_IT2.5_F + 
                      MTX63_IT5_M + MTX63_IT5_F + 
                      MTX24_IV_M + MTX24_IV_F + 
                      MTX42_IV_M + MTX42_IV_F + 
                      MTX63_IV_M + MTX63_IV_F +
                      (1|study_id)")
  crlm <- lmer(as.formula(first_form), data = med_df_MTX_ranef, REML = FALSE)
  # Add the random effect of litter to the data frame
  med_df_MTX_ranef$ID_random_effect <- ranef(crlm)$study_id[as.character(med_df_MTX_ranef$Study_ID), ]
  # Subtract the random effect from the structures volumes
  med_df_MTX_ranef[cstruct] <- med_df_MTX_ranef[cstruct] - med_df_MTX_ranef$ID_random_effect
}



mean_ci <- function(indata){ mnse <- mean_se(indata); mnse$ymin <- mnse$y-1.96*(mnse$y-mnse$ymin); mnse$ymax <- mnse$y+1.96*(mnse$ymax-mnse$y); return(mnse) }
#med_df_MTX_ranef$Treatment <- fct_rev(med_df_MTX_ranef$Treatment)
med_df_MTX_ranef$Sex <- fct_rev(med_df_MTX_ranef$Sex)
df_dose_ranef <- med_df_MTX_ranef %>% filter(Dose == 5, Route == "IV")

########################## FIG 5D - LOC ##########################
# Lateral_orbital_cortex

plot_LOC_sex <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"Lateral_orbital_cortex"],colour=Sex))+
  geom_point(alpha=0.4, aes(shape=Sex, colour=Sex), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Sex, colour=Sex), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Treatment) +
  scale_color_manual(values = c("M" = "#036E00", "F" = "#FF00DC")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Lateral Orbital Cortex") +
  theme(text=element_text(size=16,  family="Arial"), 
        axis.title = element_text(size = 16),  # Increase axis title size
        axis.text = element_text(size = 14)) + # Increase axis labels size)
  theme(legend.position = "none", axis.title.x=element_blank())

########################## FIG 5D - SC ##########################
# colliculus_superior

plot_CS_sex <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"colliculus_superior"],colour=Sex))+
  geom_point(alpha=0.4, aes(shape=Sex, colour=Sex), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Sex, colour=Sex), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Treatment) +
  scale_color_manual(values = c("M" = "#036E00", "F" = "#FF00DC")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Superior Colliculus") +
  theme(text=element_text(size=16,  family="Arial"), 
        axis.title = element_text(size = 16),  # Increase axis title size
        axis.text = element_text(size = 14)) + # Increase axis labels size)
  theme(legend.position = "none", axis.title.x=element_blank())

########################## FIG 5D - CA1 ##########################
# CA1Rad

plot_CA1Rad_sex <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"CA1Rad"],colour=Sex))+
  geom_point(alpha=0.4, aes(shape=Sex, colour=Sex), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Sex, colour=Sex), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Treatment) +
  scale_color_manual(values = c("M" = "#036E00", "F" = "#FF00DC")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="CA1 Stratum Radiatum") +
  theme(text=element_text(size=16,  family="Arial"), 
        axis.title = element_text(size = 16),  # Increase axis title size
        axis.text = element_text(size = 14)) + # Increase axis labels size)
  theme(legend.position = "none")


########################## PUT ALL PANELS TOGETHER ##########################

png(filename="/projects/schoi/MTXcomparison/figures/fig5_sex.png", width = 4, height = 8, units = 'in', res = 600)
grid.arrange(grobs = list(plot_LOC_sex, plot_CS_sex, plot_CA1Rad_sex),
             ncol=1,nrow=3)
dev.off()


########################## RUN TO OUTPUT PLOTS IN RSTUDIO ##########################

plot_LOC_sex
plot_CS_sex
plot_CA1Rad_sex

########################## TAKE LEGEND FROM HERE ##########################

png(filename="/projects/schoi/MTXcomparison/figures/fig5_sex_legend.png", width = 12, height = 12, units = 'in', res = 600)
ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"Lateral_orbital_cortex"],colour=Sex))+
  geom_point(alpha=0.4, aes(shape=Sex, colour=Sex), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Sex, colour=Sex), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Treatment) +
  scale_color_manual(values = c("M" = "#036E00", "F" = "#FF00DC")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Lateral Orbital Cortex") +
  theme(text=element_text(size=16,  family="Arial"))
dev.off()

