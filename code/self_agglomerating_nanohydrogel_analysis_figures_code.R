###### Pan-Tumor Drug Delivery with Self-Agglomerating Nanohydrogels 

# Version of April 5, 2023.
# Stephen N. Housley 
# housley.nick@gmail.com
# 

# This work is licensed under the licenses 
# Paper: Creative Commons Attribution 3.0 Unported License 
# Code: GPL-3 License 
# Depends: R (>= 3.5.0)
# Version: 0.1
# Description: code to run analytics and graphic functions associated with:
#         tumor agnostic delivery
#
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 
#

## version Hx
# v0.1- original



##### begin TEMPLATE##### 


########################### Figure ZZZZ ###########################
########################### description

########################### load dependencies
########################### Custom Functions
########################### Load Data
########################### Data Wrangling
########################### quick visualization
########################### analyses/modeling
########################### saving data
########################### saving figures
########################### Clean up


##### end TEMPLATE #####



########################### prelims ########################### 

invisible(rm(list = ls()))
invisible(gc())
setwd("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/")

########################### set dirs ########################### 
data_Name<- 'ZZZZZZZZ'   ### fill with name of file
mainDir <- "~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/" ### set this file path to the main directory 
figDir <- "Figures" 
dataDir <- "data"
dataFinalFold <- "Final"
figFolder <-"figFolder" ## this will be created if not already in existence
saveFolder <- "saveFolder" ## this will be created if not already in existence
invisible(  ifelse(!dir.exists(file.path(mainDir, figDir, figFolder)), dir.create(file.path(mainDir, figDir,figFolder)), FALSE))
invisible(  ifelse(!dir.exists(file.path(mainDir, dataDir, saveFolder)), dir.create(file.path(mainDir, dataDir,saveFolder)), FALSE))


########################### load general dependencies ########################### 
source("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/code/load_gen_dependencies.R")
rm(package.check)






########################### Figure 1 dls size distribution over batches ###########################
########################### description

########################### load dependencies

########################### Custom Functions
########################### Load Data

dls_ng_size <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_1/dls_ng_size_batches_shell_core.xlsx",  na = "NA") 


########################### Data Wrangling

dls_ng_size$batch <- factor(dls_ng_size$batch)
dls_ng_size$index <- 'dls'


########################### quick visualization

dls_batch_size_fig<-dls_ng_size %>% filter(use == "yes") %>%
  ggplot(aes(x = batch, y = Peak2_rad_nm*2, fill = core_shell)) +
  introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_x_discrete(name = "Batch", labels = c("Batch 1", "Batch 2", "Batch 3")) +
  scale_y_continuous(name = "Hydrodynamic Size (nm)",
                     breaks = seq(0, 160, 20), 
                     limits = c(60, 160)) +
  scale_fill_brewer(palette = "Blues", name = "Batch") +
  theme_minimal()

########################### analyses/modeling

dls_ng_size %>% 
  filter(use == "yes") %>%
  group_by(core_shell) %>%
  summarise(meanDiameter = mean(Peak2_rad_nm*2),
            sdDiameter = sd(Peak2_rad_nm*2))



########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("dls_batch_size_fig.pdf", width = 6, height = 4)
dls_batch_size_fig
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up



########################### Figure 1 dls size distribution over time ###########################
########################### description

########################### load dependencies

########################### Custom Functions
########################### Load Data

dls_ng_size_overTime <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_1/dls_ng_size_overTime.xlsx",  na = "NA") 


########################### Data Wrangling

dls_ng_size_overTime$temp <- factor(dls_ng_size_overTime$temp)
dls_ng_size_overTime$loadingStatus <- factor(dls_ng_size_overTime$time)


########################### quick visualization

dls_ng_size_overTime_fig<- ggplot(dls_ng_size_overTime, aes(x = time, y = percentChange, colour = loadingStatus)) +
  theme_bw() +
  geom_line() +
  geom_point()+
  facet_grid(. ~temp, scale = "free_y")+
  ylim(-25,25)+
  theme_classic()+
  scale_color_brewer(palette = "Blues")
  


########################### analyses/modeling



########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("dls_ng_size_overTime_fig.pdf", width = 6, height = 4)
dls_ng_size_overTime_fig
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up




########################### Figure 1 cell internalization ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
heyCell_allData <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_5/heyCell_allData.xlsx", 
                              na = "NA")

########################### Data Wrangling
heyCell_allData$ngConcen_mgML <- factor(heyCell_allData$ngConcen_mgML, levels=c('3', '1.5', '0.75', '0.375'))

########################### quick visualization

## concentration-dependent internalization
heyCell_dose_dependentInternalization<-heyCell_allData %>% 
  filter(inhibitor == 'na') %>%
  filter(objective == '20x') %>%
  # filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 4) %>%
  # group_by(inhibitorConcen) %>%

ggbarplot( x = "ngConcen_mgML", y = "mean",
           add = c("mean_se"),
           color = "black",
           fill = "ngConcen_mgML",
           palette = "Blues",
           width = 0.8,
           add.params = list(size = .4),
           position = position_dodge(0.4),
)+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
stat_compare_means(label.y = 7) +                                         # Global p-value
  stat_compare_means(ref.group = "3", label = "p.signif",
                     label.y = c(4.5, 3.5, 2.5))




########################### analyses/modeling
heyCell_allData %>% 
  filter(inhibitor == 'Cpz') %>%
  filter(objective == '20x') %>%
  group_by(inhibitorConcen) %>%
  summarise(n=n())

########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("heyCell_dose_dependentInternalization.pdf", width = 1.5, height = 3)
heyCell_dose_dependentInternalization
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up



########################### Figure 1 viability in vitro ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
viability_data <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_1/viability.xlsx", 
                         na = "NA")

########################### Data Wrangling
viability_data$treatment <- as.factor(viability_data$treatment)
viability_data$group <- as.factor(viability_data$group)
viability_data$sample <- as.factor(viability_data$sample)
viability_data$parameter <- as.factor(viability_data$parameter)

########################### quick visualization
# 
# rtPCR_data %>% 
#   filter(experiment == 'in_vivo') %>%
#   filter(treatment == 'kras') %>%
#   filter(species == 'mouse') %>%
#   ggbarplot( x = "group", y = "fold_change",
#              add = c("mean_se"),
#              color = "black",
#              fill = "dose",
#              palette = "Blues",
#              width = 0.3,
#              add.params = list(size = .4),
#              position = position_dodge(0.4),
#   )


viability_fig <-viability_data %>% 
  # filter(treatment == 'kras') %>%
  ggbarplot( x = "group", y = "values",
             add = c("mean_se"),
             color = "black",
             fill = "group",
             palette = "Blues",
             width = 0.3,
             add.params = list(size = .4),
             position = position_dodge(0.4),
             facet.by = "parameter"
             
  )

viability_fig<- facet(viability_fig, facet.by = "parameter", scales = "free_y")

########################### analyses/modeling
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))


cairo_pdf("viability_fig.pdf", width = 5, height = 4)
viability_fig
invisible(suppressMessages(suppressWarnings(dev.off())))



########################### Clean up







########################### Figure 2 rat biodistribution ###########################


########################### load dependencies
########################### Custom Functions 

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

########################### Load Data
biodist_rat_crc <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_2/raw_organ_biodistribution_rat.xlsx", 
                      na = c("NA", "QNS"))

########################### Data Wrangling
biodist_rat_crc <- biodist_rat_crc %>% filter(use == "yes")
biodist_rat_crc$scaledTotalRadiance <- scale_values(biodist_rat_crc$total_radiance_efficency)
biodist_rat_crc$tissue <- as.factor(biodist_rat_crc$tissue)
biodist_rat_crc$time_hrs <- factor(biodist_rat_crc$time_hrs, levels=c('control', '4', '24', '48', '288'))
biodist_rat_crc_mean<-summarySE(biodist_rat_crc, measurevar="fold_change_within_animal", groupvars=c("tissue","time_hrs"))

########################### quick visualization


# Default bar plot
figure2_rat_biodistribution<-biodist_rat_crc %>% 
  filter(tissue != "healthy_colon") %>%
  filter(time_hrs != "control") %>%
  ggbarplot( x = "tissue", y = "fold_change_within_animal", 
             add = c("mean_se", "dotplot"),
             color = "time_hrs", 
             fill = "time_hrs",
             palette = "Blues",
             add.params = list(size = .2),
             position = position_dodge(0.8))+
  theme(legend.position = "none")
  

# # Default box/whisker plot
# biodist_rat_crc %>% 
#   filter(tissue != "healthy_colon") %>%
#   filter(time_hrs != 0) %>%
#   ggplot( aes(x=as.factor(time_hrs), y=fold_change_within_animal, fill=tissue)) +
#   geom_boxplot(position=position_dodge(1))+
#   geom_dotplot(binaxis='y', stackdir='center',
#                position=position_dodge(1),
#                binwidth = 0.5)+
#   scale_fill_brewer(palette="Dark2") + 
#   theme_classic()

  
  
  
########################### analyses/modeling
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure2_rat_biodistribution.pdf", width = 6, height = 6)
figure2_rat_biodistribution
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up



 

########################### Figure 2 PK and clearance ###########################
########################### description

########################### load dependencies
########################### Custom Functions
########################### Load Data

ng_half_life <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_2/half_life_pk_cathRats.xlsx", 
                           na = c("NA", "QNS"))

########################### Data Wrangling
ng_half_life$rat <- as.factor(ng_half_life$rat)
ng_half_life$dose <- as.factor(ng_half_life$dose)
ng_half_life<-rename(ng_half_life, Time = min)
ng_half_life<-rename(ng_half_life, fluorI = value)
########################### quick visualization
data <- subset(ng_half_life, dose=="high" & rat == 1) 
time <- data$Time
conc <- data$relative 
res1 <- lee(conc=conc, time=time, method='ols', points=2, lt=TRUE)
print(res1$parms)
plot(res1)


########################### analyses/modeling
data <- subset(ng_half_life, dose=="high") 
time <- data$Time
conc <- data$relative 
auc(conc=conc, time=time, method='z', design='ssd')


data <- subset(ng_half_life, dose=="med") 
time <- data$Time
conc <- data$relative 
auc(conc=conc, time=time, method='z', design='ssd')


data <- subset(ng_half_life, dose=="low"| dose == "blank") 
time <- data$Time
conc <- data$relative 
auc(conc=conc, time=time, method='z', design='ssd')

########################### saving data
########################### saving figures
########################### Clean up


########################### Figure 2 mouse biodistribution ###########################


########################### load dependencies
########################### Custom Functions 

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

########################### Load Data
biodist_mouse_oc <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_2/raw_organ_biodistribution_mouse.xlsx", 
                              na = c("NA", "QNS"))

########################### Data Wrangling
biodist_mouse_oc <- biodist_mouse_oc %>% filter(use == "yes")
biodist_mouse_oc$scaledTotalRadiance <- scale_values(biodist_mouse_oc$total_radiance_efficency)
biodist_mouse_oc$tissue <- as.factor(biodist_mouse_oc$tissue)
biodist_mouse_oc$time_hrs <- factor(biodist_mouse_oc$time_hrs, levels=c('control', '4', '24', '48','72', '288'))
biodist_mouse_oc$tissue <- factor(biodist_mouse_oc$tissue, levels=c('ovarian_tumor','spleen_mets', 'liver_mets', 'kidney_mets','gi_mets',
                                                                    'liver', 'kidney','spleen', 'heart',   'lungs','fat', 'muscle',
                                                                     'gi_healthy'))

biodist_mouse_oc$animalType <- factor(biodist_mouse_oc$animalType, levels=c("NG" ,"control", "sirna"))
biodist_mouse_oc$ng_Control <- biodist_mouse_oc$animalType %>% fct_collapse(control = c("control","sirna"))

biodist_mouse_oc_mean<-summarySE(biodist_mouse_oc, measurevar="fold_change_within_animal", groupvars=c("tissue","time_hrs"))


########################### quick visualization


# 1
figure2_mouse_biodistribution_cancer_timeSeries<-biodist_mouse_oc %>% 
  # filter(time_hrs == "72" |time_hrs == "48" |time_hrs == "4" ) %>%
  # filter(animalType != "control") %>%
  # filter(animalType != "sirna") %>%
  filter(tumorMetNorm == "primary") %>%
  
  ggbarplot( x = "tumorMetNorm", y = "fold_change_within_animal",width = 0.3, 
             add = c("mean_se", "dotplot"),
             color = "time_hrs",
             fill = "time_hrs",
             palette = "Blues",
             add.params = list(size = .4),
             position = position_dodge(0.4),
             facet.by = "ng_Control"
  )+
  scale_y_continuous(trans='log2')+
  theme(legend.position = "none")


# 2
figure2_mouse_biodistribution_cancer_metVScontrol<-biodist_mouse_oc %>% 
  # filter(time_hrs == "72" |time_hrs == "48" |time_hrs == "4" ) %>%
  # filter(animalType != "control") %>%
  # filter(animalType != "sirna") %>%
  filter(tumorMetNorm == "mets") %>%
  
  ggbarplot( x = "tumorMetNorm", y = "fold_change_within_animal",width = 0.3, 
             add = c("mean_se", "dotplot"),
             color = "time_hrs",
             fill = "time_hrs",
             palette = "Blues",
             add.params = list(size = .4),
             position = position_dodge(0.4),
             facet.by = "ng_Control"
             )+
  theme(legend.position = "none")

# 3
figure2_mouse_biodistribution_cancer_metCompare<-biodist_mouse_oc %>% 
  # filter(time_hrs == "72" |time_hrs == "48" |time_hrs == "4" ) %>%
  filter(ng_Control != "control") %>%
  # filter(animalType != "sirna") %>%
  filter(tumorMetNorm == "mets") %>%
  
  ggbarplot( x = "tissue", y = "fold_change_within_animal",width = 0.3, 
             add = c("mean_se", "dotplot"),
             color = "time_hrs",
             fill = "time_hrs",
             palette = "Blues",
             add.params = list(size = .4),
             position = position_dodge(0.4)
             )+
  theme(legend.position = "none")

# 4
figure2_mouse_biodistribution_cancer_tissueCompareSANGS<-biodist_mouse_oc %>% 
   # filter(time_hrs == "72" |
   #          time_hrs == "48" |
   #          time_hrs == "24" |
   #          time_hrs == "4") %>%
  filter(ng_Control != "control") %>%
  # filter(animalType != "sirna") %>%
  filter(tumorMetNorm != "mets") %>%
  # 
  ggbarplot( x = "tissue", y = "fold_change_within_animal",width = 0.3, 
             add = c("mean_se", "dotplot"),
             color = "time_hrs", 
             fill = "time_hrs",
             palette = "Blues",
             add.params = list(size = .4),
             position = position_dodge(0.4),
             # facet.by = "ng_Control",
             ylim = c(0, 400)
  )+
  theme(legend.position = "none")


#5 
figure2_mouse_biodistribution_cancer_tissueCompareSIRNA<-biodist_mouse_oc %>% 
  # filter(time_hrs == "72" |
  #          time_hrs == "48" |
  #          time_hrs == "24" |
  #          time_hrs == "4") %>%
  filter(ng_Control == "control") %>%
  # filter(animalType != "sirna") %>%
  filter(tumorMetNorm != "mets") %>%
  group_by(tissue, time_hrs)%>%
  summarise(meanFC=mean(fold_change_within_animal),
            sdFC=sd(fold_change_within_animal)) 
  # 
  ggbarplot( x = "tissue", y = "fold_change_within_animal",width = 0.3, 
             add = c("mean_se", "dotplot"),
             color = "time_hrs", 
             fill = "time_hrs",
             # palette = "Purples",
             add.params = list(size = .4),
             position = position_dodge(0.4),
             # facet.by = "ng_Control",
             ylim = c(0, 400)
  )+
# figure2_mouse_biodistribution_cancer_tissueCompareSIRNA
  scale_colour_brewer(palette="Greens", direction=-1)+
  scale_fill_brewer(palette="Greens", direction=-1)
########################### analyses/modeling


biodist_mouse_oc %>%
  filter(tumorMetNorm == "primary") %>%
  # filter(animalType != "control") %>%
  group_by(ng_Control, time_hrs) %>%
  summarise( mean = mean(fold_change_within_animal),
             sd = sd(fold_change_within_animal))


########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure2_mouse_biodistribution_cancer_timeSeries.pdf", width = 3, height = 4)
figure2_mouse_biodistribution_cancer_timeSeries
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_pdf("figure2_mouse_biodistribution_cancer_metVScontrol.pdf", width = 3, height = 4)
figure2_mouse_biodistribution_cancer_metVScontrol
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_pdf("figure2_mouse_biodistribution_cancer_metCompare.pdf", width = 3, height = 4)
figure2_mouse_biodistribution_cancer_metCompare
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_pdf("figure2_mouse_biodistribution_cancer_tissueCompareSANGS.pdf", width = 8, height = 4)
figure2_mouse_biodistribution_cancer_tissueCompareSANGS
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_pdf("figure2_mouse_biodistribution_cancer_tissueCompareSIRNA.pdf", width = 8, height = 4)
figure2_mouse_biodistribution_cancer_tissueCompareSIRNA
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up





########################### Figure 2 met targeting ###########################

  ########################### load dependencies
  ########################### Custom Functions
  ########################### Load Data
met__targeting_b_ndg_OC <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_2/met_targeting_cytofluorogram.xlsx", 
                                na = c("NA", "QNS"))

  met__targeting_b_ndg_OC<- as.data.table(met__targeting_b_ndg_OC)
  
  
  ### https://www.youtube.com/watch?v=rRyJnFo57xU 
  
  ########################### Data Wrangling
  scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
  met__targeting_b_ndg_OC$bli_intensity_scaled <- scale_values(met__targeting_b_ndg_OC$bli_intensity)
  met__targeting_b_ndg_OC$ng_intensity_scaled <- scale_values(met__targeting_b_ndg_OC$ng_intensity)
  met__targeting_b_ndg_OC[, diff := bli_intensity_scaled -ng_intensity_scaled]
  met__targeting_b_ndg_OC[, diff_cat := ifelse(bli_intensity_scaled > 0.5 & ng_intensity_scaled>0.5 & abs(diff)>.15, "offTargetNG_1",
                          ifelse(bli_intensity_scaled<0.5 & ng_intensity_scaled < 0.5 & abs(diff)>0.15, "low",
                                 ifelse(ng_intensity_scaled>0.5 & bli_intensity_scaled<0.5 & abs(diff)>0.15, "offTargetNG_2",
                                        ifelse(ng_intensity_scaled<0.5 & bli_intensity_scaled>0.5 & abs(diff)>0.15, "untargeted_cancer", "onTargetNG"))))]
  ### https://stackoverflow.com/questions/52397363/r-ggplot2-ggrepel-label-a-subset-of-points-while-being-aware-of-all-points 
  
  # make plot

  color_list <- c("#3a34eb","#eb3434","#eb3434", "#000000", "#edae49")
  color_list <- c("#808080","#eb3434","#eb3434", "#2171B5", "#edae49")
  color_list <- c("#2171B5","#eb3434","#eb3434", "#6BAED6", "#edae49") ### for paper
  color_list <- c("#2171B5","#2171B5","#2171B5", "#6BAED6", "#edae49") ### for fundamental grant
  
  
  figure2_b_ndg_OC_met_targeting <- ggplot(met__targeting_b_ndg_OC, aes(x=bli_intensity_scaled,y=ng_intensity_scaled, color=diff_cat))+
    geom_point(alpha=0.7)+
    theme_classic()+
    geom_abline(intercept = 0, slope = 1, linetype="dashed")+
    geom_abline(intercept = 0.15, slope = 1)+
    geom_abline(intercept = -0.15, slope = 1)+
    scale_fill_manual(values=color_list) + 
    scale_color_manual(values=color_list)+
    theme(legend.position = c(0.85, 0.15),
          legend.background = element_rect(fill = "white", color = "black"))
  
  ### sig test
  
  ggscatter(met__targeting_b_ndg_OC, x = "bli_intensity_scaled", y = "ng_intensity_scaled", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "bli", ylab = "SANGS")
  
  ### get percentage in bins
  met__targeting_b_ndg_OC_md <- met__targeting_b_ndg_OC %>%
    group_by (diff_cat) %>%
    dplyr::summarize(count = n())%>%
    mutate(rel.freq = paste0(round(100 * count/sum(count), 4), "%"))
  

  ########################### quick visualization
  
  met__targeting_b_ndg_OC$testColorScale <- log10((met__targeting_b_ndg_OC$bli_intensity_scaled+
                                               met__targeting_b_ndg_OC$ng_intensity_scaled*10))

    ggplot(met__targeting_b_ndg_OC, aes(x=bli_intensity_scaled,y=ng_intensity_scaled, color=testColorScale))+
    geom_point(alpha=0.5)+
    theme_classic()+
    geom_abline(intercept = 0, slope = 1, linetype="dashed")+
      scale_color_viridis_c(option = "plasma",
                            direction = -1)+
    geom_abline(intercept = 0.15, slope = 1)+
    geom_abline(intercept = -0.15, slope = 1)
    
    theme(legend.position = c(0.85, 0.15),
          legend.background = element_rect(fill = "white", color = "black"))
  
  
  figure2_b_ndg_OC_met_targeting+scale_color_gradientn(colours = rainbow(5))
  ########################### analyses/modeling
  ########################### saving data
  ########################### saving figures
  
  
  setwd(file.path(mainDir,figDir,figFolder))
  
  cairo_pdf("figure2_b_ndg_OC_met_targeting.pdf", width = 6, height = 4)
  figure2_b_ndg_OC_met_targeting
  invisible(suppressMessages(suppressWarnings(dev.off())))
  
  ggsave("figure2_b_ndg_OC_met_targeting.pdf", width = 6, height = 4)
  figure2_b_ndg_OC_met_targeting
  invisible(suppressMessages(suppressWarnings(dev.off())))
  
  ggsave(figure2_b_ndg_OC_met_targeting, file="figure2_b_ndg_OC_met_targeting.png", scale=1)
  
  
  
  ########################### Clean up
  
  

  ########################### Figure 3 RNA knockdown in vitro ###########################
  
  ########################### load dependencies
  ########################### Custom Functions
  ########################### Load Data
  rtPCR_data <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_3/rtPCR_data.xlsx", 
                                na = "NA")
  
  ########################### Data Wrangling

  ########################### quick visualization
  
  ## concentration-dependent internalization
rtPCR_data %>% 
    filter(experiment == 'culture') %>%
    filter(treatment == 'kras') %>%
    ggbarplot( x = "group", y = "fold_change",
               add = c("mean_se"),
               color = "black",
               fill = "group",
               palette = "Blues",
               width = 0.8,
               add.params = list(size = .4),
               position = position_dodge(0.4),
    )+
    theme(legend.position = "none")+
    stat_compare_means(label.y = 1.2, method = "t.test") 

  
rtPCR_data %>% 
  filter(experiment == 'culture') %>%
  filter(treatment == 'egfr') %>%
  ggbarplot( x = "group", y = "fold_change",
             add = c("mean_se"),
             color = "black",
             fill = "group",
             palette = "Blues",
             width = 0.8,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  stat_compare_means(label.y = 1.2, method = "t.test") 
  
rtPCR_data %>% 
  filter(experiment == 'culture') %>%
  filter(treatment == 'zeb1') %>%
  ggbarplot( x = "group", y = "fold_change",
             add = c("mean_se"),
             color = "black",
             fill = "group",
             palette = "Blues",
             width = 0.8,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  stat_compare_means(label.y = 1.2, method = "t.test") 


inVitroSANGexpression<-rtPCR_data %>% 
  filter(experiment == 'culture') %>%
  # filter(treatment == 'zeb1') %>%
  ggbarplot( x = "group", y = "fold_change",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "group",
             palette = "Blues",
             width = 0.8,
             add.params = list(size = .4),
             position = position_dodge(0.4),
             facet.by = "treatment"
  )+
  theme(legend.position = "none")+
  stat_compare_means(label.y = 1.2, method = "t.test") 

  
  ########################### analyses/modeling

rtPCR_data %>% 
  filter(experiment == 'in_vivo') %>%
  group_by(group)%>%
  summarise(mean= mean(fold_change),
            sd = sd(fold_change))
  
cohensD = (1-0.540)/(sqrt((0.116^2+0.124^2)/2))

  ########################### saving data
  ########################### saving figures
  setwd(file.path(mainDir,figDir,figFolder))
  
  cairo_pdf("inVitroSANGexpression.pdf", width = 4, height = 4)
  inVitroSANGexpression
  invisible(suppressMessages(suppressWarnings(dev.off())))
  
  
  ########################### Clean up
  
  
  
  
  ########################### Figure 3 RNA knockdown in vivo ###########################
  
  ########################### load dependencies
  ########################### Custom Functions
  ########################### Load Data
  rtPCR_data <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_3/rtPCR_data.xlsx", 
                           na = "NA")
  
  ########################### Data Wrangling
  rtPCR_data$dose <- as.factor(rtPCR_data$dose)
  ########################### quick visualization
  # 
  # rtPCR_data %>% 
  #   filter(experiment == 'in_vivo') %>%
  #   filter(treatment == 'kras') %>%
  #   filter(species == 'mouse') %>%
  #   ggbarplot( x = "group", y = "fold_change",
  #              add = c("mean_se"),
  #              color = "black",
  #              fill = "dose",
  #              palette = "Blues",
  #              width = 0.3,
  #              add.params = list(size = .4),
  #              position = position_dodge(0.4),
  #   )

  
  inVivoSANGexpressionMouse <-rtPCR_data %>% 
    filter(experiment == 'in_vivo') %>%
    # filter(treatment == 'kras') %>%
    filter(species == 'mouse') %>%
    ggbarplot( x = "group", y = "fold_change",
               add = c("mean_se"),
               color = "black",
               fill = "dose",
               palette = "Blues",
               width = 0.3,
               add.params = list(size = .4),
               position = position_dodge(0.4),
               facet.by = "treatment"
               
    )
  
  inVivoSANGexpressionRat <-rtPCR_data %>% 
    filter(experiment == 'in_vivo') %>%
    # filter(treatment == 'kras') %>%
    filter(species == 'rat') %>%
    ggbarplot( x = "group", y = "fold_change",
               add = c("mean_se"),
               color = "black",
               fill = "group",
               palette = "Blues",
               width = 0.3,
               ylim = c(0, 1.3),
               add.params = list(size = .4),
               position = position_dodge(0.4),
               facet.by = "treatment"
               

    )
  ########################### analyses/modeling
  ########################### saving data
  ########################### saving figures
  setwd(file.path(mainDir,figDir,figFolder))
  
  cairo_pdf("inVivoSANGexpressionMouse.pdf", width = 6, height = 4)
  inVivoSANGexpressionMouse
  invisible(suppressMessages(suppressWarnings(dev.off())))
  
  cairo_pdf("inVivoSANGexpressionRat.pdf", width = 4, height = 4)
  inVivoSANGexpressionRat
  invisible(suppressMessages(suppressWarnings(dev.off())))
  
  
  
  ########################### Clean up
  
  
  
  
  
  
########################### Figure 4 toxicity cd1###########################


########################### load dependencies
########################### Custom Functions
########################### Load Data
tox_cd1 <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_4_tox/cd_1_weights.xlsx", 
                             na = c("NA", "QNS"))
colnames(tox_cd1)

########################### Data Wrangling
### group summary stats
tox_cd1_mean <- tox_cd1 %>%
  filter(!is.na(weight)) %>%
  group_by(day, treatment) %>%
  summarise(n = n(),
            mean = mean(weight),
            median = median(weight),
            sd = sd(weight)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

### % change from baseline for individual
tox_cd1 <- tox_cd1 %>% group_by(animalNum) %>% arrange(animalNum, day) %>% 
  mutate(change = weight - first(weight))

tox_cd1_percent_change <- tox_cd1 %>%
  filter(!is.na(change)) %>%
  group_by(day, treatment) %>%
  summarise(n = n(),
            mean = mean(change),
            median = median(change),
            sd = sd(change)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)


########################### quick visualization
color_list <- c("#eb3434","#3a34eb")


figure4_cd1_mean <-ggplot(tox_cd1_mean, aes(x=day, y=mean, color = treatment)) +
  geom_line(aes(x=day, y=mean, color=treatment)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=treatment),color=NA,alpha=0.4)+
  # theme_light(base_size = 16) + 
  theme_classic()+
  xlim(0,15) + 
  ylim(16,36) + 
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values=color_list) + 
  scale_color_manual(values=color_list)+
  labs(title = "weight agross time")+
  xlab("Days post final infusion ")+
  ylab("weight (g)")+
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))


setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure4_cd1_mean.pdf", width = 6, height = 4)
figure4_cd1_mean
invisible(suppressMessages(suppressWarnings(dev.off())))



figure4_cd1_percentChange <- ggplot(tox_cd1_percent_change, aes(x=day, y=mean, color = treatment)) +
  geom_line(aes(x=day, y=mean, color=treatment)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=treatment),color=NA,alpha=0.4)+
  # theme_light(base_size = 16) + 
  theme_classic()+
  xlim(0,15) + 
  ylim(-5, 5)+
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values=color_list) + 
  scale_color_manual(values=color_list)+
  labs(title = "% weight change")+
  xlab("Days post final infusion ")+
  ylab("% weight change from baseline")+
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))
  


########################### analyses/modeling
########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure4_cd1_percentChange.pdf", width = 6, height = 4)
figure4_cd1_percentChange
invisible(suppressMessages(suppressWarnings(dev.off())))




########################### Clean up



########################### Figure 4 repeat dose tolerability rat###########################
########################### load dependencies
########################### Custom Functions
########################### Load Data
tox_rat <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_4_tox/rat_repeatDose_weights.xlsx", 
                      na = c("NA", "QNS"))
colnames(tox_rat)


########################### Data Wrangling
### group summary stats

tox_rat_df<- tox_rat %>% 
  filter(day < 15 )  %>%
  filter(animalNum != "R0024_1") %>%
  filter(animalNum != "R0025_0")
  
tox_rat_mean <- tox_rat_df %>%
  select(day, weight,treatment) %>% 
  filter(weight != "NA") %>% 
  group_by(treatment,day) %>%
  summarise(n = n(),
            mean = mean(weight),
            median = median(weight),
            sd = sd(weight)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

### % change from baseline for individual
tox_rat_df <- tox_rat_df %>% group_by(animalNum) %>% arrange(animalNum, day) %>% 
  mutate(change = weight - first(weight))

tox_rat_percent_change <- tox_rat_df %>%
  filter(!is.na(change)) %>%
  group_by(day, treatment) %>%
  summarise(n = n(),
            mean = mean(change),
            median = median(change),
            sd = sd(change)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)


########################### quick visualization
color_list <- c("#3a34eb", "#808080", "#eb3434")

figure4_rat_percentChange <- 
  tox_rat_percent_change %>%
  ggplot(aes(x=day, y=mean, color = treatment)) +
  geom_line(aes(x=day, y=mean, color=treatment)) +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=treatment),color=NA,alpha=0.4)+
  # theme_light(base_size = 16) + 
  theme_classic()+
  # xlim(0,15) + 
  ylim(-20, 20)+
  # geom_vline(xintercept = 1)+
  scale_fill_manual(values=color_list) + 
  scale_color_manual(values=color_list)+
  labs(title = "% weight change")+
  xlab("Days post final infusion ")+
  ylab("% weight change from baseline")+
  theme(legend.position = c(0.17, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure4_rat_percentChange.pdf", width = 6, height = 4)
figure4_rat_percentChange
invisible(suppressMessages(suppressWarnings(dev.off())))





########################### Figure 4 MTD rat ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data

ng_tox_mtd_rat <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_4_tox/ng_tox_mtd_rat.xlsx", 
                             na = c("NA", "QNS"))

########################### Data Wrangling

ng_tox_mtd_rat$collection_time_hr <- as.factor(ng_tox_mtd_rat$collection_time_hr)
ng_tox_mtd_rat$sample_number <- as.factor(ng_tox_mtd_rat$sample_number)
ng_tox_mtd_rat$group <- as.factor(ng_tox_mtd_rat$group)
ng_tox_mtd_rat$animal_num_glp <- as.factor(ng_tox_mtd_rat$animal_num_glp)

ng_tox_mtd_rat$dose[ng_tox_mtd_rat$group==1] <- 'control'
ng_tox_mtd_rat$dose[ng_tox_mtd_rat$group==2] <- '7mg/kg'
ng_tox_mtd_rat$dose[ng_tox_mtd_rat$group==3] <- '12mg/kg'
ng_tox_mtd_rat$dose[ng_tox_mtd_rat$group==4] <- '17mg/kg'
ng_tox_mtd_rat$dose <- factor(ng_tox_mtd_rat$dose, levels=c('control', '7mg/kg', '12mg/kg', '17mg/kg'))

#### 
ng_tox_mtd_rat_df<-
  ng_tox_mtd_rat %>% 
  filter( animal_num_glp != "986") %>%
  filter( animal_num_glp != "983") %>%
  filter( animal_num_glp != "960") %>%
  filter( animal_num_glp != "984")


ng_tox_mtd_rat_sum <- ng_tox_mtd_rat_df %>%
  # select(day, weight,treatment) %>% 
  filter(value != "NA") %>% 
  group_by(dose,collection_time_hr, parameter) %>%
  summarise(n = n(),
            sd = sd(value),
            value = mean(value)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = value + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = value - qt((1-0.95)/2, n - 1) * sem)


########################### quick visualization
#### bun ###

color_list <- c("#3a34eb","#eb3434")

figure4_rat_serum_tox<-ggplot(ng_tox_mtd_rat_df, aes(dose, value, color = collection_time_hr)) +
  facet_wrap(~parameter,nrow=2, scales="free")+
  geom_col(data = ng_tox_mtd_rat_sum, position = position_dodge(0.8), 
           width = 0.7, alpha=0.4, aes(fill = collection_time_hr)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +

  geom_errorbar(
    aes(ymin = value-sd, ymax = value+sd), data = ng_tox_mtd_rat_sum, 
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values=color_list) + 
  scale_color_manual(values=color_list)+
  theme_classic()+
  stat_n_text(y.pos = -0.01, size = 3 )+
  stat_compare_means(label = "p.format", size = 2)
 

########################### analyses/modeling
########################### saving data

########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure4_rat_serum_tox.pdf", width = 10, height = 6)
figure4_rat_serum_tox
invisible(suppressMessages(suppressWarnings(dev.off())))




########################### Figure 4 single dose tolerability cd1 ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data

ng_tox_cd1 <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_4_tox/ng_tox_cd_1.xlsx", 
                             na = c("NA", "QNS"))

########################### Data Wrangling

ng_tox_cd1$collection_time_hr <- as.factor(ng_tox_cd1$collection_time_hr)
ng_tox_cd1$sample_number <- as.factor(ng_tox_cd1$sample_number)
ng_tox_cd1$group <- as.factor(ng_tox_cd1$group)
ng_tox_cd1$animal_num_glp <- as.factor(ng_tox_cd1$animal_num_glp)


#### 
ng_tox_cd1_df<-
  ng_tox_cd1
  
  
ng_tox_cd1_sum <- ng_tox_cd1_df %>%
  # select(day, weight,treatment) %>% 
  filter(value != "NA") %>% 
  group_by(group,collection_time_hr, parameter) %>%
  dplyr::summarise(
            sd = sd(value),
            value = mean(value),
            n=n()) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = value + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = value - qt((1-0.95)/2, n - 1) * sem)


########################### quick visualization
#### bun ###

color_list <- c("#3a34eb","#eb3434")

figure4_mouse_serum_tox<-ggplot(ng_tox_cd1_df, aes(group, value, color = group)) +
  facet_wrap(~parameter,nrow=2, scales="free")+
  geom_col(data = ng_tox_cd1_sum, position = position_dodge(0.2), 
           width = 0.3, alpha=0.4, aes(fill = group)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  
  geom_errorbar(
    aes(ymin = value-sd, ymax = value+sd), data = ng_tox_cd1_sum, 
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values=color_list) + 
  scale_color_manual(values=color_list)+
  theme_classic()+
  stat_n_text(y.pos = -0.01, size = 3 )+
  stat_compare_means(label = "p.format", size = 2)+
  theme(text=element_text(family="Arial"))

########################### analyses/modeling
########################### saving data

########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure4_mouse_serum_tox.pdf", width = 5, height = 4)
figure4_mouse_serum_tox
invisible(suppressMessages(suppressWarnings(dev.off())))






########################### Figure 4 single dose tox porcine ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data

ng_tox_pig <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_4_tox/ng_tox_porcine.xlsx", 
                         na = c("NA", "QNS"))

########################### Data Wrangling

ng_tox_pig$parameter <- as.factor(ng_tox_pig$parameter)




########################### quick visualization
#### bun ###

color_list <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
# 
# rect2 <- data.frame(xmin = c(0,0,0,0),
#                     xmax = c(360,360,360,360),
#                     ymin = c(15,22,160,0.8),
#                     ymax = c(55,47,450,2.3),
#                     parameter = c("ast", "alt", "ldh", "creatinine"))
# 
# figure4_pig_serum_tox<-ggplot() + 
#   geom_rect(data = rect2 , aes(xmin = xmin,
#                                xmax = xmax,
#                                ymin = ymin,
#                                ymax = ymax,
#                                fill = parameter),
#             alpha = 0.2) +
#   geom_line(data = ng_tox_pig, aes(x = collection_time_hr, y = value, colour=parameter), size = 1)+
#   facet_wrap(~parameter, scales="free",nrow=4)+
#   labs(x = "time (mins)",
#        y = "Serum Value",
#        title = "Acute Toxicity in Porcine"
#   ) +
#   scale_fill_manual(values=color_list) + 
#   scale_color_manual(values=color_list)+
#   theme(text=element_text(family="Arial"))+
#   theme_classic()+
#   theme(legend.position = "none")
  


g0 <- ggplot(ng_tox_pig, aes(collection_time_hr, value, colour=parameter)) + 
  geom_line(size =1) + 
  facet_grid(rows = vars(parameter), scales = "free")

facet_bounds <- read.table(header=TRUE,
                           text=                           
                             "parameter ymin ymax breaks
alt     0     100    20
ast     0     100    20
creatinine     0    5    1
ldh     0    500    100",
                           stringsAsFactors=FALSE)

ff <- with(facet_bounds,
           data.frame(value=c(ymin,ymax),
                      parameter=c(parameter,parameter)))


figure4_pig_serum_tox<- g0 + geom_point(data=ff,x=NA)+
  labs(x = "time (mins)",
       y = "Serum Value",
       title = "Acute Toxicity in Porcine"
  ) +
  scale_fill_manual(values=color_list) + 
  scale_color_manual(values=color_list)+
  theme(text=element_text(family="Arial"))+
  theme_classic()+
  theme(legend.position = "none")




########################### analyses/modeling
########################### saving data

########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("figure4_pig_serum_tox.pdf", width = 5, height = 4)
figure4_pig_serum_tox
invisible(suppressMessages(suppressWarnings(dev.off())))







########################### Figure 5 mechanism of entry ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
heyCell_allData <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_5/heyCell_allData.xlsx", 
                              na = "NA")

########################### Data Wrangling
heyCell_allData$ngConcen_mgML <- factor(heyCell_allData$ngConcen_mgML, levels=c('3', '1.5', '0.75', '0.375'))

########################### quick visualization

## sodium-azide-dependent internalization
heyCell_ATPdependent<-heyCell_allData %>% 
  filter(inhibitor == 'Sodium Azide') %>%
  filter(objective == '63x') %>%
  # filter(inhibitorConcen == 0 & mean > 2 | inhibitorConcen == 50) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
  stat_compare_means(label.y = 7) +                                         # Global p-value
  stat_compare_means(ref.group = "0", label = "p.signif",
                     label.y = c(5.5, 3.5))


heyCell_LatA<-heyCell_allData %>% 
  filter(inhibitor == 'LatA') %>%
  filter(objective == '20x') %>%
  filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 1|  inhibitorConcen == 2) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
  stat_compare_means(label.y = 15) +                                         # Global p-value
  stat_compare_means(ref.group = "0", label = "p.signif",
                     label.y = c(10, 10))


heyCell_cytoD<-heyCell_allData %>% 
  filter(inhibitor == 'cytoD') %>%
  filter(objective == '20x') %>%
  filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 2|  inhibitorConcen == 4) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
  stat_compare_means(label.y = 15) +                                         # Global p-value
  stat_compare_means(ref.group = "0", label = "p.signif",
                     label.y = c(10, 10))


heyCell_Cpz<-heyCell_allData %>% 
  filter(inhibitor == 'Cpz') %>%
  filter(objective == '20x') %>%
  filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 50|  inhibitorConcen == 100) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
  stat_compare_means(label.y = 15) +                                         # Global p-value
  stat_compare_means(ref.group = "0", label = "p.signif",
                     label.y = c(11, 11))


heyCell_Cpz<-heyCell_allData %>% 
  filter(inhibitor == 'Cpz') %>%
  filter(objective == '20x') %>%
  mutate_at(vars("intDen"), funs(./1000)) %>%
  
  # filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 100|  inhibitorConcen == 200) %>%
  ggbarplot( x = "inhibitorConcen", y = "intDen",
             add = c("mean_se"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
  stat_compare_means(label.y = 15) +                                         # Global p-value
  stat_compare_means(ref.group = "0", label = "p.signif",
                     label.y = c(11, 11))

heyCell_MBcd<-heyCell_allData %>% 
  filter(inhibitor == 'MBcd') %>%
  filter(objective == '20x') %>%
  filter(intDen < 30000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 0.5| inhibitorConcen == 1) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  # stat_compare_means(aes(group = ngConcen_mgML), label = "p.signif", label.y = 20)
  stat_compare_means(label.y = 15) +                                         # Global p-value
  stat_compare_means(ref.group = "0", label = "p.signif",
                     label.y = c(11, 11))



########################### analyses/modeling
########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,suppFigFolder))


cairo_pdf("heyCell_ATPdependent.pdf", width = 1.5, height = 3)
heyCell_ATPdependent
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_pdf("heyCell_LatA.pdf", width = 1.5, height = 3)
heyCell_LatA
invisible(suppressMessages(suppressWarnings(dev.off())))

cairo_pdf("heyCell_cytoD.pdf", width = 1.5, height = 3)
heyCell_cytoD
invisible(suppressMessages(suppressWarnings(dev.off())))

cairo_pdf("heyCell_Cpz.pdf", width = 1.5, height = 3)
heyCell_Cpz
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_pdf("heyCell_MBcd.pdf", width = 1.5, height = 3)
heyCell_MBcd
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up



########################### Figure 5 siRNA release ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
sirnaEndosomColoc_allData <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_5/heyCellsiRNAendosomeCOLOC.xlsx", 
                              na = "NA")

########################### Data Wrangling
sirnaEndosomColoc_allData$cell <- factor(sirnaEndosomColoc_allData$cell)
sirnaEndosomColoc_allData$metric <- factor(sirnaEndosomColoc_allData$metric)
########################### quick visualization
siRNA_endosome_Fig<- sirnaEndosomColoc_allData %>% 
  ggplot(aes(x = factor(metric), y = pearsonCorr, fill = factor(metric))) +
  
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    # ploting on left side
    side = "left",
    # adjusting position
    justification = 1.1,
    # adjust grouping (binning) of observations
    binwidth = 0.025
  ) +
  # Themes and Labels
  scale_fill_tq() +
  theme_tq() +
  labs(
    title = "colocalization of siRNA and Endosomes",
    x = "",
    y = "Pearson R",
    fill = ""
  ) +
  coord_flip()+
  theme(legend.position = "none")


########################### analyses/modeling
sirnaEndosomColoc_allData %>%
  summarise(meanColoc = mean(pearsonCorr),
            sdColoc = sd(pearsonCorr))
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))


cairo_pdf("siRNA_endosome_Fig.pdf", width = 5, height = 3)
siRNA_endosome_Fig
invisible(suppressMessages(suppressWarnings(dev.off())))



########################### Clean up




########################### Figure 5 siRNA release ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
heyCellsiRNAsangsCOLOC_allData <- read_excel("data/figure_5/heyCellsiRNAsangsCOLOC.xlsx", 
                                        na = "NA")

########################### Data Wrangling
heyCellsiRNAsangsCOLOC_allData$cell <- factor(heyCellsiRNAsangsCOLOC_allData$cell)
heyCellsiRNAsangsCOLOC_allData$metric <- factor(heyCellsiRNAsangsCOLOC_allData$metric)
heyCellsiRNAsangsCOLOC_allData$stage <- factor(heyCellsiRNAsangsCOLOC_allData$stage)

########################### quick visualization
siRNA_sang_Fig<- heyCellsiRNAsangsCOLOC_allData %>% 
  # filter(metric != "pearsonCorr") %>%
  filter(metric == "ColocCoeff1") %>%
  ggbarplot( x = "metric", y = "value",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "stage",
             palette = "Blues",
             width = 0.3,
             ylim = c(0, 1.3),
             add.params = list(size = .4),
             position = position_dodge(0.4)
            
  )

########################### analyses/modeling



########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))


cairo_pdf("siRNA_sang_Fig.pdf", width = 4, height = 4)
siRNA_sang_Fig
invisible(suppressMessages(suppressWarnings(dev.off())))



########################### Clean up





########################### Figure 5 tumor weights ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
combo_allData <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_5/combo_data.xlsx", 
                                             na = "NA")

########################### Data Wrangling
combo_allData$group <- factor(combo_allData$group)


########################### quick visualization
tumor_weight_Fig<- combo_allData %>% 
  filter(group == "cis" | group == "simi" | group == "nc") %>%
  ggbarplot( x = "group", y = "tumor_weight",
             add = c("mean_se"),
             color = "black",
             fill = "group",
             palette = "Blues",
             width = 0.3,
             ylim = c(0, 2),
             add.params = list(size = .4),
             position = position_dodge(0.4)
             
  )

########################### analyses/modeling
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))


cairo_pdf("tumor_weight_Fig.pdf", width = 3, height = 4)
tumor_weight_Fig
invisible(suppressMessages(suppressWarnings(dev.off())))



########################### Clean up






########################### Figure 6 in vivo colocalization size ###########################


########################### load dependencies
########################### Custom Functions
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

########################### Load Data
Results_colocalization <- read_csv("~/GaTech Dropbox/Nick Housley/papers_dropbox/tumor_agnostic_delivery/data/figure_6/Results_colocalization.csv", 
                                   na = c("NA", "QNS", "NaN"))


########################### Data Wrangling

########################### quick visualization

# ggplot(Results_colocalization, aes(x = colocalized, y = Area, fill = colocalized)) +
#   geom_violindot() +
#   theme_modern()+
#   scale_y_continuous(trans='log2')


sang_colocalization_size_Fig<-Results_colocalization %>%
  ggplot(aes(x= factor(colocalized), y = Area, fill = factor(colocalized)))+
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -0.2,
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.colour = NA,
    alpha = 0.5
  ) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1,
    dotsize = 2,
    binwidth =0.06
  ) +
  scale_fill_tq()+
  theme_tq()+
  coord_flip()+
  scale_y_continuous(trans='log2',
                     # labels = function(x) format(x, scientific = TRUE)
                     )+
  theme(legend.position = "none")+
  xlab( "colocalized")



########################### analyses/modeling
Results_colocalization %>%
  group_by(colocalized) %>%
  summarise(mean = min(Area))

  ggplot(aes(x= factor(colocalized), y = Area, fill = factor(colocalized)))+
  
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))


cairo_pdf("sang_colocalization_size_Fig.pdf", width = 6, height = 17)
sang_colocalization_size_Fig
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up





########################### Figure 6 DOSY concentration
########################### Figure 6 dls concentration###########################
########################### description

########################### load dependencies

########################### Custom Functions
########################### Load Data

dls_ng_wyatt_concentration<- read_excel("~/Dropbox-GaTech/Nick Housley/papers_dropbox/tumor_agnostic_delivery/data/figure_6/dls_wyatt_concentration.xlsx",  na = "NA") 
dls_ng_zetasizer_concentration<- read_excel("~/Dropbox-GaTech/Nick Housley/papers_dropbox/tumor_agnostic_delivery/data/figure_6/dls_zetaSizer_concentration.xlsx",  na = "NA") 
dls_ng_concentration<- read_excel("~/Dropbox-GaTech/Nick Housley/papers_dropbox/tumor_agnostic_delivery/data/figure_6/dls_concentration_both.xlsx",  na = "NA") 

########################### Data Wrangling

dls_ng_wyatt_concentration$acquisition <- factor(dls_ng_wyatt_concentration$acquisition)
dls_ng_wyatt_concentration$concentration <- factor(dls_ng_wyatt_concentration$concentration)

dls_ng_zetasizer_concentration$acquisition <- factor(dls_ng_zetasizer_concentration$acquisition)
dls_ng_zetasizer_concentration$concentration <- factor(dls_ng_zetasizer_concentration$concentration)

dls_ng_concentration$acquisition <- factor(dls_ng_concentration$acquisition)
dls_ng_concentration$concentration <- factor(dls_ng_concentration$concentration)
dls_ng_concentration$instrument <- factor(dls_ng_concentration$instrument)

########################### quick visualization

dls_ng_wyatt_concentration %>% 
  group_by(concentration, diameter) %>%
  summarise(intensity = mean(intensity)) %>%
  ggplot(aes(x = diameter, y = intensity, color = concentration)) +
  theme_bw() +
  geom_line() +
  # geom_point()+
  theme_classic()+
  # scale_x_continuous(trans='log10')+
  # scale_x_continuous(trans='log10',breaks = c(1,10,100,1000, 10000), limits = c(0, 10000))+
  theme(legend.position = "none")+
  scale_color_brewer(palette = "Blues")+
  xlab( "diameter")+
  scale_x_log10(limits = c(1,1e5))
  
  
dls_ng_zetasizer_fig<-dls_ng_zetasizer_concentration %>% 
    group_by(concentration, diameter) %>%
    summarise(intensity = mean(intensity)) %>%
    ggplot(aes(x = diameter, y = intensity, color = concentration)) +
    theme_bw() +
    geom_line() +
    # geom_point()+
    theme_classic()+
    scale_x_continuous(trans='log10', breaks = c(1,10,100,1000, 10000))+
    theme(legend.position = "none")+
    scale_color_brewer(palette = "Blues")+  
    xlab( "diameter")+
    scale_x_log10(limits = c(1,1e5))
  
  
  

  
  dls_ng_concentration %>% 
    group_by(instrument, concentration, diameter) %>%
    summarise(intensity = mean(intensity)) %>%
    ggplot(aes(x = diameter, y = intensity, color = interaction(instrument, concentration, sep=':'))) +
    theme_bw() +
    geom_line() +
    # geom_point()+
    theme_classic()+
    scale_x_continuous(trans='log10')+
    theme(legend.position = "none")+
    scale_color_brewer(palette = "Blues")
  xlab( "diameter")
  
  

########################### analyses/modeling



########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("dls_ng_zetasizer_fig.pdf", width = 6, height = 4)
dls_ng_zetasizer_fig
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up





########################### Figure 6 dosy ###########################
########################### description

########################### load dependencies

########################### Custom Functions
########################### Load Data

dosy_data <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_6/dosy.xlsx",  na = "NA") 


########################### Data Wrangling

dosy_data$metric <- factor(dosy_data$metric)


########################### quick visualization

dls_ng_size_overTime_fig<- ggplot(dosy_data, aes(x = nanomolar, y = diffusion_coefficient, colour = metric)) +
  theme_bw() +
  geom_line() +
  geom_point()+
  facet_grid(. ~metric, scale = "free_y")+
  ylim(-25,25)+
  theme_classic()+
  scale_color_brewer(palette = "Blues")


dosy_fig<-dosy_data %>% 
  # filter(treatment == 'kras') %>%
  ggline( x = "nanomolar", y = "diffusion_coefficient",
             # add = c("mean_se"),
             color = "black",
             fill = "metric",
             palette = "Blues",
             width = 0.3,
             add.params = list(size = .4),
             position = position_dodge(0.4),
             facet.by = "metric"
             
  )

dosy_fig<- facet(dosy_fig, facet.by = "metric", scales = "free_y")


########################### analyses/modeling



########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("dosy_fig.pdf", width = 6, height = 4)
dosy_fig

invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up







########################### Figure 7 dls size concentrations full range ###########################
########################### description

########################### load dependencies

########################### Custom Functions
########################### Load Data

dls_ng_size <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/data/figure_6/dls_concentration_full_curve.xlsx",  na = "NA") 


########################### Data Wrangling

dls_ng_size$acquisition <- factor(dls_ng_size$acquisition)
dls_ng_size$radius_nm <- factor(dls_ng_size$radius_nm)
dls_ng_size$Item <- factor(dls_ng_size$Item)
dls_ng_size$acquisition <- factor(dls_ng_size$acquisition, levels=c('ng_0.4', 'ng_0.9', 'ng_1.875', 'ng_3.75', 'ng_7.5', 'ng_15', 'ng_60', 'ng_60_redlute_15', 'ng_60_redlute_3'))


########################### quick visualization

dls_batch_size_fig<-dls_ng_size %>% 
  # filter(radius_nm == "radius_range_2") %>%
  ggplot(aes(x = acquisition, y = value*2, fill = radius_nm)) +
  introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  # geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  # scale_x_discrete(name = "Batch", labels = c("Batch 1", "Batch 2", "Batch 3")) +
  # scale_y_continuous(name = "Hydrodynamic Size (nm)",
  #                    breaks = seq(0, 160, 20), 
  #                    limits = c(60, 160)) +
  scale_y_log10(limits = c(10,1e6))+
  scale_fill_brewer(palette = "Blues", name = "Batch") +
  theme_minimal()

########################### analyses/modeling

dls_summary<-dls_ng_size %>% 
  group_by(acquisition, radius_nm) %>%
  # summarise(meanDiameter = median(value*2),
  #           sdDiameter = confint(value))
  summarise(mean.diam = mean(value, na.rm = TRUE),
          sd.diam = sd(value, na.rm = TRUE),
          n.diam = n()) %>%
  mutate(se.diam = sd.diam / sqrt(n.diam),
         lower.ci.diam = mean.diam - qt(1 - (0.05 / 2), n.diam - 1) * se.diam,
         upper.ci.diam = mean.diam + qt(1 - (0.05 / 2), n.diam - 1) * se.diam)


dls_concentration_size_fig<- ggplot(dls_summary, aes(x=acquisition, y=mean.diam, group=radius_nm, color=radius_nm)) + 
  # geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean.diam-se.diam, ymax=mean.diam+se.diam), width=.2,
                position=position_dodge(0.05))+
  scale_color_brewer(palette = "Blues", name = "radius_nm") +
  theme_classic()+
  theme(legend.position="none")
  
  


########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("dls_concentration_size_fig.pdf", width = 6, height = 4)
dls_concentration_size_fig
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up




########################### Supplemental Figure  TEM vs DLS  ###########################


########################### Load Data

## read tem data 
tem_ng_size <- read_excel("~/Dropbox-GaTech/Nick Housley/papers_dropbox/tumor_agnostic_delivery/data/figure_1/tem_ng_size.xlsx", 
                          na = "NA")


########################### Data Wrangling
tem_ng_size$Radius_nm <- tem_ng_size$Area/2
tem_ng_size$index <- 'tem'
########################### analyses/modeling

tem_ng_size %>% 
  group_by(index) %>%
  summarise(meanSize = mean(Radius_nm*2),
            sdSize = sd(Radius_nm*2))


########################### saving data
########################### saving figures

setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("dls_batch_size_fig.pdf", width = 6, height = 4)
dls_batch_size_fig
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up





########################### Supp Figure 5 cell internalization ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
heyCell_allData <- read_excel("~/Dropbox (GaTech)/papers_dropbox/tumor_agnostic_delivery/figures/supplemental/Supplementary Figure 5. Internalization/time/Book1.xlsx", 
                              na = "NA")

########################### Data Wrangling
heyCell_allData$slideNum <- factor(heyCell_allData$slideNum, levels=c('3', '1.5', '0.75', '0.375'))

########################### quick visualization

## concentration-dependent internalization
heyCell_allData %>% 
  filter(slideNum != 'background') %>%

  ggbarplot( x = "incubationTime", y = "mean",
             add = c("mean_se",  "jitter"),
             color = "black",
             fill = "incubationTime",
             palette = "Blues",
             width = 0.8,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")+
  stat_compare_means(label.y = 27) +                                         # Global p-value
  stat_compare_means(ref.group = "1", label = "p.signif",
                     label.y = c(10, 13, 15, 22, 18))


my_comparisons <- list( c("0", "1"),
                        c("1", "3"), 
                        c("3", "6"), 
                        c("6", "12"), c("6", "18"), 
                        c("12", "18"),
                        c("12", "24"),
                        c("18", "24")
)

heyCell_time_dependentInternalization<- heyCell_allData %>% 
  filter(slideNum != 'background') %>%
ggboxplot(x = "incubationTime", 
          y = "mean",
          color = "incubationTime", 
          palette = "Blues",
          add = "jitter")+
          # shape = "incubationTime")+ 
  theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" )+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50, method = "anova")  

########################### analyses/modeling
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))

cairo_pdf("heyCell_time_dependentInternalization.pdf", width = 6, height = 12)
heyCell_time_dependentInternalization
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up


