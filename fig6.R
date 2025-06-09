### Clinical mapping analysis


########################## LOAD PACKAGES ##########################

library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(RMINC)
library(MRIcrotome)
library(reshape2)
library(glue)
library(dplyr)

# Load median data
load("/projects/schoi/MTXcomparison/analysis/MTX17022024/med_MTX.RData")

# Load clinical data
clinical <- read.csv("/projects/schoi/MTXcomparison/analysis/aug2024analysis/QoL_volume_data.csv")

# Load rough mappings
map <- read.csv("/projects/schoi/MTXcomparison/analysis/aug2024analysis/Human-Mouse_Mapping.csv")

# Delta ALL divided by CTRL Vol times 100 to get vol change in human data
# Calculate human changes
# Delta.ALL / VCTL * 100
clinical$volChange <- clinical$Delta.ALL / clinical$VCTL * 100

clinical_df <- med_df_MTX

# merging groups of structures based on map
clinical_df$cerebellum_group <- clinical_df$flocculus_FL + clinical_df$paraflocculus_PFL + clinical_df$copula_pyramis_lobule_8 + 
  clinical_df$crus_1_ansiform_lobule_lobule_6 + clinical_df$crus_2_ansiform_lobule_lobule_7 + clinical_df$paramedian_lobule_lobule_7 +
  clinical_df$simple_lobule_lobule_6 + clinical_df$anterior_lobule_lobules_4_5 + clinical_df$lobule_10_nodulus +
  clinical_df$lobule_3_central_lobule_dorsal + clinical_df$lobule_6_declive + clinical_df$lobule_7_tuber_or_folium +
  clinical_df$lobule_8_pyramis + clinical_df$lobule_9_uvula + clinical_df$lobules_1_2_lingula_and_central_lobule_ventral +
  clinical_df$lobules_4_5_culmen_ventral_and_dorsal + clinical_df$flocculus_white_matter + clinical_df$paraflocculus_white_matter +
  clinical_df$lobule_1_2_white_matter + clinical_df$lobule_10_white_matter + clinical_df$lobule_3_white_matter +
  clinical_df$lobule_8_white_matter + clinical_df$lobule_9_white_matter + clinical_df$lobules_4_5_white_matter +
  clinical_df$lobules_6_7_white_matter + clinical_df$trunk_of_lobules_1_3_white_matter + clinical_df$trunk_of_lobules_6_8_white_matter + 
  clinical_df$copula_white_matter + clinical_df$crus_1_white_matter + clinical_df$crus_2_white_matter +
  clinical_df$paramedian_lobule + clinical_df$simple_lobule_white_matter + clinical_df$trunk_of_arbor_vita +
  clinical_df$trunk_of_crus_2_and_paramedian_white_matter + clinical_df$trunk_of_simple_and_crus_1_white_matter +
  clinical_df$dentate_nucleus + clinical_df$fastigial_nucleus + clinical_df$nucleus_interpositus + clinical_df$anterior_lobule_white_matter

clinical_df$hippocampus_group <- clinical_df$CA10r + clinical_df$CA1Py + clinical_df$CA1Rad + clinical_df$CA20r +
  clinical_df$CA2Py + clinical_df$CA2Rad + clinical_df$CA30r + clinical_df$CA3Py_Inner + clinical_df$CA3Py_Outer + 
  clinical_df$CA3Rad + clinical_df$GrDG + clinical_df$LMol + clinical_df$MoDG + clinical_df$PoDG + clinical_df$SLu + clinical_df$subiculum

clinical_df$thalamus_group <- clinical_df$thalamus

clinical_df$amygdala_group <- clinical_df$amygdala + clinical_df$Medial_amygdala

clinical_df$striatum_group <- clinical_df$striatum

clinical_df$brainstem_group <- clinical_df$pons + clinical_df$medulla + clinical_df$midbrain + clinical_df$colliculus_inferior + clinical_df$colliculus_superior

clinical_df$frontallobe_group <- clinical_df$Frontal_association_cortex + clinical_df$Frontal_cortex_area_3 + clinical_df$Dorsolateral_orbital_cortex +
  clinical_df$Lateral_orbital_cortex + clinical_df$Medial_orbital_cortex + clinical_df$Ventral_orbital_cortex + 
  clinical_df$Primary_motor_cortex + clinical_df$Secondary_motor_cortex

clinical_df$occipitallobe_group <- clinical_df$Primary_visual_cortex + clinical_df$Primary_visual_cortex_binocular_area + 
  clinical_df$Primary_visual_cortex_monocular_area + clinical_df$Secondary_visual_cortex_lateral_area +
  clinical_df$Secondary_visual_cortex_mediolateral_area + clinical_df$Secondary_visual_cortex_mediomedial_area

clinical_df$parietallobe_group <- clinical_df$Lateral_parietal_association_cortex + clinical_df$Medial_parietal_association_cortex +
  clinical_df$Parietal_cortex_posterior_area_rostral_part + clinical_df$Primary_somatosensory_cortex + clinical_df$Primary_somatosensory_cortex_barrel_field +
  clinical_df$Primary_somatosensory_cortex_dysgranular_zone + clinical_df$Primary_somatosensory_cortex_forelimb_region +
  clinical_df$Primary_somatosensory_cortex_hindlimb_region + clinical_df$Primary_somatosensory_cortex_jaw_region +
  clinical_df$Primary_somatosensory_cortex_shoulder_region + clinical_df$Primary_somatosensory_cortex_trunk_region +
  clinical_df$Primary_somatosensory_cortex_upper_lip_region + clinical_df$Secondary_somatosensory_cortex

clinical_df$temporallobe_group <- clinical_df$Primary_auditory_cortex + clinical_df$Secondary_auditory_cortex_dorsal_area +
  clinical_df$Secondary_auditory_cortex_ventral_area + clinical_df$Temporal_association_area

clinical_df$cingulatecortex_group <- clinical_df$Cingulate_cortex_area_24a + clinical_df$Cingulate_cortex_area_24a.1 +
  clinical_df$Cingulate_cortex_area_24b + clinical_df$Cingulate_cortex_area_24b.1 + clinical_df$Cingulate_cortex_area_25 +
  clinical_df$Cingulate_cortex_area_29a + clinical_df$Cingulate_cortex_area_29b + clinical_df$Cingulate_cortex_area_29c +
  clinical_df$Cingulate_cortex_area_30 + clinical_df$Cingulate_cortex_area_32

clinical_df$parahippocampalcortex_group <- clinical_df$Amygdalopiriform_transition_area + clinical_df$Cortex_amygdala_transition_zones +
  clinical_df$Piriform_cortex + clinical_df$Posterolateral_cortical_amygdaloid_area + clinical_df$Posteromedial_cortical_amygdaloid_area +
  clinical_df$Rostral_amygdalopiriform_area + clinical_df$Caudomedial_entorhinal_cortex + clinical_df$Dorsal_intermediate_entorhinal_cortex +
  clinical_df$Dorsal_tenia_tecta + clinical_df$Dorsolateral_entorhinal_cortex + clinical_df$Medial_entorhinal_cortex +
  clinical_df$Ventral_intermediate_entorhinal_cortex + clinical_df$Ventral_tenia_tecta

clinical_df$insula_group <- clinical_df$Insular_region_not_subdivided + clinical_df$Ectorhinal_cortex + clinical_df$Perirhinal_cortex +
  clinical_df$Dorsal_nucleus_of_the_endopiriform + clinical_df$Intermediate_nucleus_of_the_endopiriform_claustrum +
  clinical_df$Ventral_nucleus_of_the_endopiriform_claustrum

clinical_df$cingulum_group <- clinical_df$Cingulum

clinical_df$globuspallidus_group <- clinical_df$globus_pallidus

clinical_df$cerebralpeduncle_group <- clinical_df$cerebral_peduncle

clinical_df$internalcapsule_group <- clinical_df$internal_capsule

clinical_df$corpuscallosum_group <- clinical_df$corpus_callosum

clinical_df$fornix_group <- clinical_df$fornix + clinical_df$fimbria

clinical_df$mediallemniscusmediallongitudinalfasciculus_group <- clinical_df$medial_lemniscus_medial_longitudinal_fasciculus

clinical_df$corticospinaltractpyramids_group <- clinical_df$corticospinal_tract_pyramids

clinical_df$cerebellarpeduncle_group <- clinical_df$cerebellar_peduncle_inferior + clinical_df$cerebellar_peduncle_middle + clinical_df$cerebellar_peduncle_superior

clinical_df$ventricles_group <- clinical_df$cerebral_aqueduct + clinical_df$fourth_ventricle + clinical_df$lateral_ventricle + clinical_df$third_ventricle

##############################################

### Human volumes ### have to do these manually for structure of interest
#amygdala
amygdalaVCTL <- clinical[18, "VCTL"] + clinical[32, "VCTL"]
amygdalaDelta.ALL <- clinical[18, "Delta.ALL"] + clinical[32, "Delta.ALL"]
amygdala_change <- (amygdalaDelta.ALL / amygdalaVCTL) * 100

#corpus callosum
CCVCTL <- clinical[21, "VCTL"] + clinical[34, "VCTL"] + clinical[38, "VCTL"] + clinical[3, "VCTL"] + clinical[4, "VCTL"]
CCDelta.ALL <- clinical[21, "Delta.ALL"] + clinical[34, "Delta.ALL"] + clinical[38, "Delta.ALL"] + clinical[3, "Delta.ALL"] + clinical[4, "Delta.ALL"]
CC_change <- (CCDelta.ALL / CCVCTL) * 100

#hippocampus
hippocampusVCTL <- clinical[80, "VCTL"] + clinical[82, "VCTL"]
hippocampusDelta.ALL <- clinical[80, "Delta.ALL"] + clinical[82, "Delta.ALL"]
hippocampus_change <- (hippocampusDelta.ALL / hippocampusVCTL) * 100

#thalamus
thalamusVCTL <- clinical[46, "VCTL"] + clinical[10, "VCTL"]
thalamusDelta.ALL <- clinical[46, "Delta.ALL"] + clinical[10, "Delta.ALL"]
thalamus_change <- (thalamusDelta.ALL / thalamusVCTL) * 100

#cerebellum
cerebellumVCTL <- clinical[16, "VCTL"] + clinical[44, "VCTL"]
cerebellumDelta.ALL <- clinical[16, "Delta.ALL"] + clinical[44, "Delta.ALL"]
cerebellum_change <- (cerebellumDelta.ALL / cerebellumVCTL) * 100

#striatum
striatumVCTL <- clinical[40, "VCTL"] + clinical[52, "VCTL"] + clinical[59, "VCTL"] + clinical[74, "VCTL"] + clinical[22, "VCTL"] + clinical[37, "VCTL"]
striatumDelta.ALL <- clinical[40, "Delta.ALL"] + clinical[52, "Delta.ALL"] + clinical[59, "Delta.ALL"] + clinical[74, "Delta.ALL"] + clinical[22, "Delta.ALL"] + clinical[37, "Delta.ALL"]
striatum_change <- (striatumDelta.ALL / striatumVCTL) * 100

#brainstem
brainstemVCTL <- clinical[71, "VCTL"]
brainstemDelta.ALL <- clinical[71, "Delta.ALL"]
brainstem_change <- (brainstemDelta.ALL / brainstemVCTL) * 100

#frontal lobe
frontallobeVCTL <- clinical[39, "VCTL"] + clinical[47, "VCTL"]
frontallobeDelta.ALL <- clinical[39, "Delta.ALL"] + clinical[47, "Delta.ALL"]
frontallobe_change <- (frontallobeDelta.ALL / frontallobeVCTL) * 100

#occipital lobe
occipitallobeVCTL <- clinical[75, "VCTL"] + clinical[81, "VCTL"]
occipitallobeDelta.ALL <- clinical[75, "Delta.ALL"] + clinical[81, "Delta.ALL"]
occipitallobe_change <- (occipitallobeDelta.ALL / occipitallobeVCTL) * 100

#parietal lobe
parietallobeVCTL <- clinical[41, "VCTL"] + clinical[55, "VCTL"]
parietallobeDelta.ALL <- clinical[41, "Delta.ALL"] + clinical[55, "Delta.ALL"]
parietallobe_change <- (parietallobeDelta.ALL / parietallobeVCTL) * 100

#temporal lobe
temporallobeVCTL <- clinical[48, "VCTL"] + clinical[56, "VCTL"]
temporallobeDelta.ALL <- clinical[48, "Delta.ALL"] + clinical[56, "Delta.ALL"]
temporallobe_change <- (temporallobeDelta.ALL / temporallobeVCTL) * 100

#cingulate cortex
cingulatecortexVCTL <- clinical[14, "VCTL"] + clinical[33, "VCTL"]
cingulatecortexDelta.ALL <- clinical[14, "Delta.ALL"] + clinical[33, "Delta.ALL"]
cingulatecortex_change <- (cingulatecortexDelta.ALL / cingulatecortexVCTL) * 100

#parahippocampal cortex
parahippocampalcortexVCTL <- clinical[42, "VCTL"] + clinical[45, "VCTL"]
parahippocampalcortexDelta.ALL <- clinical[42, "Delta.ALL"] + clinical[45, "Delta.ALL"]
parahippocampalcortex_change <- (parahippocampalcortexDelta.ALL / parahippocampalcortexVCTL) * 100

#insula
insulaVCTL <- clinical[17, "VCTL"] + clinical[54, "VCTL"]
insulaDelta.ALL <- clinical[17, "Delta.ALL"] + clinical[54, "Delta.ALL"]
insula_change <- (insulaDelta.ALL / insulaVCTL) * 100

#cingulum
cingulumVCTL <- clinical[31, "VCTL"] + clinical[58, "VCTL"] + clinical[77, "VCTL"] + clinical[78, "VCTL"]
cingulumDelta.ALL <- clinical[31, "Delta.ALL"] + clinical[58, "Delta.ALL"] + clinical[77, "Delta.ALL"] + clinical[78, "Delta.ALL"]
cingulum_change <- (cingulumDelta.ALL / cingulumVCTL) * 100

#globus pallidus
globuspallidusVCTL <- clinical[20, "VCTL"] + clinical[25, "VCTL"]
globuspallidusDelta.ALL <- clinical[20, "Delta.ALL"] + clinical[25, "Delta.ALL"]
globuspallidus_change <- (globuspallidusDelta.ALL / globuspallidusVCTL) * 100

#cerebral peduncle
cerebralpeduncleVCTL <- clinical[12, "VCTL"] + clinical[27, "VCTL"]
cerebralpeduncleDelta.ALL <- clinical[12, "Delta.ALL"] + clinical[27, "Delta.ALL"]
cerebralpeduncle_change <- (cerebralpeduncleDelta.ALL / cerebralpeduncleVCTL) * 100

#internal capsule
internalcapsuleVCTL <- clinical[1, "VCTL"] + clinical[2, "VCTL"] + clinical[15, "VCTL"] + clinical[30, "VCTL"]
internalcapsuleDelta.ALL <- clinical[1, "Delta.ALL"] + clinical[2, "Delta.ALL"] + clinical[15, "Delta.ALL"] + clinical[30, "Delta.ALL"]
internalcapsule_change <- (internalcapsuleDelta.ALL / internalcapsuleVCTL) * 100

#fornix
fornixVCTL <- clinical[13, "VCTL"] + clinical[26, "VCTL"]
fornixDelta.ALL <- clinical[13, "Delta.ALL"] + clinical[26, "Delta.ALL"]
fornix_change <- (fornixDelta.ALL / fornixVCTL) * 100

#medial lemniscus
mediallemniscusVCTL <- clinical[70, "VCTL"] + clinical[73, "VCTL"]
mediallemniscusDelta.ALL <- clinical[70, "Delta.ALL"] + clinical[73, "Delta.ALL"]
mediallemniscus_change <- (mediallemniscusDelta.ALL / mediallemniscusVCTL) * 100

#corticospinal tract
corticospinaltractVCTL <- clinical[57, "VCTL"] + clinical[63, "VCTL"]
corticospinaltractDelta.ALL <- clinical[57, "Delta.ALL"] + clinical[63, "Delta.ALL"]
corticospinaltract_change <- (corticospinaltractDelta.ALL / corticospinaltractVCTL) * 100

#cerebellar peduncles
cerebellarpedunclesVCTL <- clinical[65, "VCTL"] + clinical[67, "VCTL"] + clinical[68, "VCTL"]
cerebellarpedunclesDelta.ALL <- clinical[65, "Delta.ALL"] + clinical[67, "Delta.ALL"] + clinical[68, "Delta.ALL"]
cerebellarpeduncles_change <- (cerebellarpedunclesDelta.ALL / cerebellarpedunclesVCTL) * 100

#ventricles
ventriclesVCTL <- clinical[87, "VCTL"] + clinical[88, "VCTL"] + clinical[89, "VCTL"] + clinical[90, "VCTL"]
ventriclesDelta.ALL <- clinical[87, "Delta.ALL"] + clinical[88, "Delta.ALL"] + clinical[89, "Delta.ALL"] + clinical[90, "Delta.ALL"]
ventricles_change <- (ventriclesDelta.ALL / ventriclesVCTL) * 100



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





########################## FIGURE 6A: NEUROANATOMICAL MAPS ##########################

anatVolClinical <- mincArray(mincGetVolume('/projects/bjnieman/QoL/figures/refimg/t1_final_avg.mnc'))
mask <- '/projects/bjnieman/QoL/figures/refimg/t1_segment_animal_voted.mnc'

mappings<-'/projects/schoi/MTXcomparison/fig6brian/animalLabelsNew.csv'
dt.struc <- read.csv(mappings)
structnames <- dt.struc$Structure

# 0.5

df0.5 <- read.csv("/projects/schoi/MTXcomparison/fig6brian/df0.5new.csv")

overlay0.5 <- struct_wise_to_vxl_wise(df0.5$Volume, mincGetVolume(mask),
                                      structure.names=structnames, defs=mappings)

map0.5 <- sliceSeries(begin=200, end = 85, ncol=2, nrow=3) %>% 
  anatomy(anatVolClinical, low=100, high=800) %>% 
  overlay(mincArray(overlay0.5), low=1, high=5, symmetric=T) %>%
  addtitle("IT 0.5 mg/kg") %>%
  MRIcrotome::legend('Volume Change (%)') %>%
  contourSliceIndicator(anatVolClinical, c(800, 1400)) %>%
  draw()


# IT 5.0

df5 <- read.csv("/projects/schoi/MTXcomparison/fig6brian/df5new.csv")

overlay5 <- struct_wise_to_vxl_wise(df5$Volume, mincGetVolume(mask),
                                    structure.names=structnames, defs=mappings)

map5 <- sliceSeries(begin=200, end = 85, ncol=2, nrow=3) %>% 
  anatomy(anatVolClinical, low=100, high=800) %>% 
  overlay(mincArray(overlay5), low=1, high=5, symmetric=T) %>%
  addtitle("IT 5.0 mg/kg") %>%
  MRIcrotome::legend('Volume Change (%)') %>%
  contourSliceIndicator(anatVolClinical, c(800, 1400)) %>%
  draw()


# IV 5.0

dfIV <- read.csv("/projects/schoi/MTXcomparison/fig6brian/dfIVnew.csv")

overlayIV <- struct_wise_to_vxl_wise(dfIV$Volume, mincGetVolume(mask),
                                     structure.names=structnames, defs=mappings)

mapIV <- sliceSeries(begin=200, end = 85, ncol=2, nrow=3) %>% 
  anatomy(anatVolClinical, low=100, high=800) %>% 
  overlay(mincArray(overlayIV), low=1, high=5, symmetric=T) %>%
  addtitle("IV 5.0 mg/kg") %>%
  MRIcrotome::legend('Volume Change (%)') %>%
  contourSliceIndicator(anatVolClinical, c(800, 1400)) %>%
  draw()


# IV 5.0

overlayClinical <- struct_wise_to_vxl_wise(dfIV$Clinical_Percent_Change, mincGetVolume(mask),
                                           structure.names=structnames, defs=mappings)

mapIV <- sliceSeries(begin=200, end = 85, ncol=2, nrow=3) %>% 
  anatomy(anatVolClinical, low=100, high=800) %>% 
  overlay(mincArray(overlayClinical), low=1, high=5, symmetric=T) %>%
  addtitle("Clinical") %>%
  MRIcrotome::legend('Volume Change (%)') %>%
  contourSliceIndicator(anatVolClinical, c(800, 1400)) %>%
  draw()


########################## FIGURE 6B: HEATMAP ##########################

##############################################
# MAKING TABLE

# make lists and loop through them to generate tables

# Vectors
structures <- c("Cerebellum", "Hippocampus", "Thalamus", "Amygdala", "Striatum", "Brain Stem", "Frontal Lobe",
                "Occipital Lobe", "Parietal Lobe", "Temporal Lobe", "Cingulate Cortex", "Parahippocampal Cortex",
                "Insula", "Cingulum", "Globus Pallidus", "Cerebral Peduncle", "Internal Capsule", "Corpus Callosum",
                "Fornix", "Medial Lemniscus", "Corticospinal Tract", "Cerebellar Peduncles", "Ventricles")
#IT0.5
IT0.5_df <- clinical_df %>%
  filter(Route == "IT", Dose == 0.5)

# Initialize a vector to store percent changes
IT0.5_changes <- numeric()

# Loop through columns 199 to 221
for (col in 199:221) {
  # Get the column name
  col_name <- colnames(IT0.5_df)[col]
  
  # Calculate mean values for MTX and Saline
  mean_treatment <- mean(IT0.5_df[IT0.5_df$Treatment == "MTX", col], na.rm = TRUE)
  mean_control <- mean(IT0.5_df[IT0.5_df$Treatment == "Saline", col], na.rm = TRUE)
  
  # Calculate percent change and add to the vector
  percent_change <- (mean_treatment - mean_control) / mean_control * 100
  IT0.5_changes <- c(IT0.5_changes, percent_change)
}

# Set names of percent_changes vector to match column names
names(IT0.5_changes) <- colnames(IT0.5_df)[199:221]

# View the result
IT0.5_changes


#IT5
IT5_df <- clinical_df %>%
  filter(Route == "IT", Dose == 5)

# Initialize a vector to store percent changes
IT5_changes <- numeric()

# Loop through columns 199 to 221
for (col in 199:221) {
  # Get the column name
  col_name <- colnames(IT5_df)[col]
  
  # Calculate mean values for MTX and Saline
  mean_treatment <- mean(IT5_df[IT5_df$Treatment == "MTX", col], na.rm = TRUE)
  mean_control <- mean(IT5_df[IT5_df$Treatment == "Saline", col], na.rm = TRUE)
  
  # Calculate percent change and add to the vector
  percent_change <- (mean_treatment - mean_control) / mean_control * 100
  IT5_changes <- c(IT5_changes, percent_change)
}

# Set names of percent_changes vector to match column names
names(IT5_changes) <- colnames(IT5_df)[199:221]


# View the result
IT5_changes

#IV5
IV5_df <- clinical_df %>%
  filter(Route == "IV", Dose == 5)

# Initialize a vector to store percent changes
IV5_changes <- numeric()

# Loop through columns 199 to 221
for (col in 199:221) {
  # Get the column name
  col_name <- colnames(IV5_df)[col]
  
  # Calculate mean values for MTX and Saline
  mean_treatment <- mean(IV5_df[IV5_df$Treatment == "MTX", col], na.rm = TRUE)
  mean_control <- mean(IV5_df[IV5_df$Treatment == "Saline", col], na.rm = TRUE)
  
  # Calculate percent change and add to the vector
  percent_change <- (mean_treatment - mean_control) / mean_control * 100
  IV5_changes <- c(IV5_changes, percent_change)
}

# Set names of percent_changes vector to match column names
names(IV5_changes) <- colnames(IV5_df)[199:221]

# View the result
IV5_changes


clinical_change <- c(cerebellum_change, hippocampus_change, thalamus_change, amygdala_change, striatum_change, brainstem_change, 
                     frontallobe_change, occipitallobe_change, parietallobe_change, temporallobe_change, cingulatecortex_change, 
                     parahippocampalcortex_change, insula_change, cingulum_change, globuspallidus_change, cerebralpeduncle_change, 
                     internalcapsule_change, CC_change, fornix_change, mediallemniscus_change, corticospinaltract_change, 
                     cerebellarpeduncles_change, ventricles_change)

# Create a data frame
table <- data.frame(Structure = structures, IT0.5_Percent_Change = IT0.5_changes, IT5_Percent_Change = IT5_changes, 
                    IV5_Percent_Change = IV5_changes, Clinical_Percent_Change = clinical_change)

table$Type <- c("GM", "GM", "GM", "GM", "GM", "GM", "GM", "GM", "GM", "GM", "GM", "GM", 
                "GM", "WM", "GM", "WM", "WM", "WM", "WM", "WM", "WM", "WM", "NA")



filtered_table <- table %>%
  filter(Structure != "Ventricles")


# Heatmap

library(pheatmap)
df_cluster <- t(filtered_table)
colnames(df_cluster) <- as.character(df_cluster[1,])
df_cluster_new <- df_cluster[c(-1,-6),]
dftest <- data.frame(df_cluster_new)
df_cluster_numeric <-  mutate_all(dftest, function(x) as.numeric(as.character(x)))

dist_mat <- dist(df_cluster_numeric, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')


rownames(df_cluster_numeric) <- c("IT 0.5 mg/kg (Mouse)", "IT 5.0 mg/kg (Mouse)", "IV 5.0 mg/kg (Mouse)", "Clinical")

# Define color breaks including 0
neg_breaks <- seq(-5, -1.01, length.out = 49)
zero_break <- 0  # Exact zero
pos_breaks <- seq(1.01, 5, length.out = 49)
breaks <- c(neg_breaks, zero_break, pos_breaks)

# Colors: blue gradient for negative, white for zero, yellow-red for positive
neg_colors <- colorRampPalette(c("cyan", "blue"))(length(neg_breaks) - 1)
zero_color <- "white"
pos_colors <- colorRampPalette(c("red", "yellow"))(length(pos_breaks) - 1)
color_palette <- c(neg_colors, zero_color, pos_colors)

# Create heatmap
pheatmap(as.matrix(df_cluster_numeric),
         scale = "none",
         fontsize = 14,
         color = color_palette,
         breaks = breaks,
         border_color = 'grey')
