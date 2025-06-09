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
    print(paste0("WE ENTERED THE FOR LOOP YAY! Stage:", j))
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


#############################################################
# Generation of image panel displays for all figures
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



########################## FIG 4A: NEUROANATOMICAL MAPS ########################## 

#image mappings output
for (croute in c("IT")){
  Lgrobs=list()
  for (cdose in c("0.5", "1", "2.5", "5")){
    for (ctpt in c("24")){
      ccolref<-paste0("MTX",ctpt,"_", croute, cdose)
      if (!paste0(ccolref,"_Pr(>|t|)")%in%colnames(Mres)){ 
        ccolref<-paste0("MTX",ctpt,"_", croute, cdose)
      }
      Valpha<-ifelse( Mres[,paste0(ccolref,"_Pr(>|t|)")]>0.05 ,0.0,
                      ifelse( (Mres[,paste0(ccolref,"_Pr(>|t|)")]<0.05)&(Mres[,paste0(ccolref,"_qvalue")]>=0.1) ,0.5,1.0))[structure_list]
      #names(Valpha)<-names(structure_list)
      overlayVol<-genOverlayArray(names(Valpha),
                                  (100*Mres[,paste0(ccolref,"_Estimate")]/Mres[,paste0("age_flag",ctpt,":Dose", cdose, "_Estimate")])[structure_list],
                                  man_colVol,manRcolVol,statLowVol,statHighVol,Valpha,labelVol,labelmappings)
      print("We passed overlay!")
      lGrpST <- sliceSeries(begin=220, end = 90, ncol=2, nrow=6) %>% 
        anatomy(anatVol, low=700, high=1400) %>% 
        overlay(overlayVol) %>% 
        grobify()
      Lgrobs[[ccolref]]<-lGrpST
    }
  }
  Nwht<-round(2*(statLowVol/statHighVol)*length(man_RcolVol)/(1.0-statLowVol/statHighVol))
  lgndtxt<-c("0.5"="0.5 mg/kg (%)",
             "1.0"="1.0 mg/kg (%)",
             "2.5"="2.5 mg/kg (%)",
             "5.0"="5.0 mg/kg (%)")[cdose]
  llgnd<-genHorizLegend( c(rev(man_RcolVol),rep("#FFFFFF",Nwht),man_colVol),lgndtxt,-statHighVol,statHighVol,
                         textsize=5,topM=0.6,bottomM=-0.8,limtextvposadj=0.1,lbltextvposadj=-0.1)
  colHead1<-textGrob("IT MTX 0.5 mg/kg",gp=gpar(fontsize=18))
  colHead2<-textGrob("IT MTX 1.0 mg/kg",gp=gpar(fontsize=18))
  colHead3<-textGrob("IT MTX 2.5 mg/kg",gp=gpar(fontsize=18))
  colHead4<-textGrob("IT MTX 5.0 mg/kg",gp=gpar(fontsize=18))
  # Save the output to a PDF
  #pdf(file = paste0("/projects/schoi/MTXcomparison/figures/", "fig3", croute, "v2.pdf"), width = 12, height = 8)
  # Save the output to a PNG with resolution specification
  png(file = paste0("/projects/schoi/MTXcomparison/figures/", "fig3", croute, "v2.png"), width = 12, height = 8, units = "in", res = 300)
  grid.arrange(grobs=list(colHead1,colHead2,colHead3,colHead4,
                          Lgrobs[[grep("MTX24_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT1",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT2.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT5",names(Lgrobs))]],
                          llgnd),
               widths=c(2,2,2,2),
               heights=c(0.2,1.5,0.05,0.3),
               ncol=4,nrow=4,
               layout_matrix = rbind(
                 c(1, 2, 3, 4),  # Column headers
                 c(5, 6, 7, 8),  # Plots for each dosage
                 c(NA, NA, NA, NA),  # Empty row for spacing (optional, can be removed)
                 c(9, 9, 9, 9))   # Legend spanning all columns
  )
  dev.off()
}


# Run this to display in R plots
grid.arrange(grobs=list(colHead1,colHead2,colHead3,colHead4,
                        Lgrobs[[grep("MTX24_IT0.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT1",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT2.5",names(Lgrobs))]],Lgrobs[[grep("MTX24_IT5",names(Lgrobs))]],
                        llgnd),
             widths=c(2,2,2,2),
             heights=c(0.2,1.5,0.05,0.3),
             ncol=4,nrow=4,
             layout_matrix = rbind(
               c(1, 2, 3, 4),  # Column headers
               c(5, 6, 7, 8),  # Plots for each dosage
               c(NA, NA, NA, NA),  # Empty row for spacing (optional, can be removed)
               c(9, 9, 9, 9))   # Legend spanning all columns
)




# For legend generation

fig.ITdosecoronal<- sliceSeries(nrow = 4, ncol = 1, begin=220, end = 90,dimension = 2)%>%
  anatomy(anatVol, low = 700, high = 1400) %>% 
  overlay(mincArray(plot_MTX24_IT5*100), low=1, high=15, symmetric =T, col = man_colVol, man_RcolVol) %>%
  addtitle('5.0 mg/kg') %>% 
  MRIcrotome::legend('% Volume Change') %>% 
  contourSliceIndicator(anatVol, c(700, 1400)) %>% draw()


########################## FIG 4B: DENSITY PLOTS ##########################

#P24
### Fig 3, dose analysis - plot % change for each dose (normalize to controls)
change_MTX24_IT0.5 <- (Mres[, "MTX24_IT0.5_Estimate"])/(Mres[, "age_flag24:Dose0.5_Estimate"]) *100
change_MTX24_IT1 <- (Mres[, "MTX24_IT1_Estimate"])/(Mres[, "age_flag24:Dose1_Estimate"]) *100
change_MTX24_IT2.5 <- (Mres[, "MTX24_IT2.5_Estimate"])/(Mres[, "age_flag24:Dose2.5_Estimate"]) *100
change_MTX24_IT5 <- (Mres[, "MTX24_IT5_Estimate"])/(Mres[, "age_flag24:Dose5_Estimate"]) *100


# Density lines instead of histogram with all doses on same plot

# Set the background to white (mimicking theme_bw)
par(bg = "white")

# Adjust margins to make room for the legend
par(mar = c(5, 4, 4, 4) + 0.1)  # Increase the right margin

# Plot the density line
plot(density(change_MTX24_IT0.5), xlim = c(-25, 15), ylim = c(0, 0.24), main = "P24", 
     xlab = "% Volume Change", ylab = "Proportion of Structures", col = "blue", lwd = 5,
     cex.main = 1.5,  # Increase title size
     cex.lab = 1.4,   # Increase axis label size
     cex.axis = 1.2)   # Increase tick label size

# Add the second density line
lines(density(change_MTX24_IT1), col = "red", lwd = 5)

# Add the third density line
lines(density(change_MTX24_IT2.5), col = "purple", lwd = 5)

# Add the fourth density line
lines(density(change_MTX24_IT5), col = "green", lwd = 5)

# Optional: Adjust gridlines to make the plot cleaner
grid(lwd = 0.5, col = "gray60")

### Legend appears on lines so just copy and paste onto figure at end
# Optional: Add legend to distinguish between the lines
graphics::legend("right", legend = c("0.5 mg/kg", "1.0 mg/kg", "2.5 mg/kg", "5.0 mg/kg"), 
                 col = c("blue", "red", "purple", "green"), lwd = 5, xpd=FALSE)


#P63
### Fig 3, dose analysis - plot % change for each dose (normalize to controls)
change_MTX63_IT0.5 <- (Mres[, "MTX63_IT0.5_Estimate"])/(Mres[, "age_flag63:Dose0.5_Estimate"]) *100
change_MTX63_IT1 <- (Mres[, "MTX63_IT1_Estimate"])/(Mres[, "age_flag63:Dose1_Estimate"]) *100
change_MTX63_IT2.5 <- (Mres[, "MTX63_IT2.5_Estimate"])/(Mres[, "age_flag63:Dose2.5_Estimate"]) *100
change_MTX63_IT5 <- (Mres[, "MTX63_IT5_Estimate"])/(Mres[, "age_flag63:Dose5_Estimate"]) *100


# Density lines instead of histogram with all doses on same plot

# Set the background to white (mimicking theme_bw)
par(bg = "white")

# Adjust margins to make room for the legend
par(mar = c(5, 4, 4, 4) + 0.1)  # Increase the right margin

# Plot the density line
plot(density(change_MTX63_IT0.5), xlim = c(-25, 15), ylim = c(0, 0.24), main = "P63", 
     xlab = "% Volume Change", ylab = "Proportion of Structures", col = "blue", lwd = 5,
     cex.main = 1.5,  # Increase title size
     cex.lab = 1.4,   # Increase axis label size
     cex.axis = 1.2)   # Increase tick label size

# Add the second density line
lines(density(change_MTX63_IT1), col = "red", lwd = 5)

# Add the third density line
lines(density(change_MTX63_IT2.5), col = "purple", lwd = 5)

# Add the fourth density line
lines(density(change_MTX63_IT5), col = "green", lwd = 5)

# Optional: Adjust gridlines to make the plot cleaner
grid(lwd = 0.5, col = "gray60")

### Legend appears on lines so just copy and paste onto figure at end
# Optional: Add legend to distinguish between the lines
graphics::legend("right", legend = c("0.5 mg/kg", "1.0 mg/kg", "2.5 mg/kg", "5.0 mg/kg"), 
                 col = c("blue", "red", "purple", "green"), lwd = 5, xpd=FALSE)




### CALCULATING PERCENT VOLUME CHANGE FOR FIGURE 4A

change_MTX24_IT0.5 <- (Mres[, "MTX24_IT0.5_Estimate"])/(Mres[, "age_flag24:Dose0.5_Estimate"]) *100
change_MTX24_IT1 <- (Mres[, "MTX24_IT1_Estimate"])/(Mres[, "age_flag24:Dose1_Estimate"]) *100
change_MTX24_IT2.5 <- (Mres[, "MTX24_IT2.5_Estimate"])/(Mres[, "age_flag24:Dose2.5_Estimate"]) *100
change_MTX24_IT5 <- (Mres[, "MTX24_IT5_Estimate"])/(Mres[, "age_flag24:Dose5_Estimate"]) *100

change_MTX42_IT0.5 <- (Mres[, "MTX42_IT0.5_Estimate"])/(Mres[, "age_flag42:Dose0.5_Estimate"]) *100
change_MTX42_IT1 <- (Mres[, "MTX42_IT1_Estimate"])/(Mres[, "age_flag42:Dose1_Estimate"]) *100
change_MTX42_IT2.5 <- (Mres[, "MTX42_IT2.5_Estimate"])/(Mres[, "age_flag42:Dose2.5_Estimate"]) *100
change_MTX42_IT5 <- (Mres[, "MTX42_IT5_Estimate"])/(Mres[, "age_flag42:Dose5_Estimate"]) *100

change_MTX63_IT0.5 <- (Mres[, "MTX63_IT0.5_Estimate"])/(Mres[, "age_flag63:Dose0.5_Estimate"]) *100
change_MTX63_IT1 <- (Mres[, "MTX63_IT1_Estimate"])/(Mres[, "age_flag63:Dose1_Estimate"]) *100
change_MTX63_IT2.5 <- (Mres[, "MTX63_IT2.5_Estimate"])/(Mres[, "age_flag63:Dose2.5_Estimate"]) *100
change_MTX63_IT5 <- (Mres[, "MTX63_IT5_Estimate"])/(Mres[, "age_flag63:Dose5_Estimate"]) *100




# Perform Mann-Whitney U tests (Wilcoxon rank-sum test) for each comparison for Figure 4B

# P24
mw_test_5_vs_2.5_P24 <- wilcox.test(change_MTX24_IT5, change_MTX24_IT2.5, alternative = "less")
mw_test_2.5_vs_1_P24 <- wilcox.test(change_MTX24_IT2.5, change_MTX24_IT1, alternative = "less")
mw_test_1_vs_0.5_P24 <- wilcox.test(change_MTX24_IT1, change_MTX24_IT0.5, alternative = "less")

# Extract p-values
p_values_24 <- c(mw_test_5_vs_2.5_P24$p.value, 
              mw_test_2.5_vs_1_P24$p.value, 
              mw_test_1_vs_0.5_P24$p.value)

# Adjust p-values for multiple comparisons using Bonferroni correction
adjusted_p_values_24 <- p.adjust(p_values_24, method = "bonferroni")

# Print adjusted p-values
print(adjusted_p_values_24)



# P63
mw_test_5_vs_2.5_P63 <- wilcox.test(change_MTX63_IT5, change_MTX63_IT2.5, alternative = "less")
mw_test_2.5_vs_1_P63 <- wilcox.test(change_MTX63_IT2.5, change_MTX63_IT1, alternative = "less")
mw_test_1_vs_0.5_P63 <- wilcox.test(change_MTX63_IT1, change_MTX63_IT0.5, alternative = "less")

# Extract p-values
p_values_P63 <- c(mw_test_5_vs_2.5_P63$p.value, 
              mw_test_2.5_vs_1_P63$p.value, 
              mw_test_1_vs_0.5_P63$p.value)

# Adjust p-values for multiple comparisons using Bonferroni correction
adjusted_p_values_63 <- p.adjust(p_values_P63, method = "bonferroni")

# Print adjusted p-values
print(adjusted_p_values_63)

