###### Housley et al 2025. Tumor Agnostic Drug Delivery with Dynamic Nanohydrogels

########################### ReadME begin ########################### 
# Version of April 5, 2023.
# Stephen N. Housley 
# nickhousley@gatech.edu

# This work is licensed under the licenses 
# Paper: Creative Commons Attribution 3.0 Unported License 
# Code: GPL-3 License 
# Depends: R (>= 3.5.0)
# Version: 0.1
# Description: code to run analytics and graphic functions associated with:
#         Housley et al 2025. Tumor Agnostic Drug Delivery with Dynamic Nanohydrogels Nature Communications
#         
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 

## version Hx
# v0.1- original April 5, 2023.
# v1.1- original May 5, 2024. - publish to github
# v1.2- original Oct 21, 2025. - publish to github

# paths must be changed to accommodate end user file structure (e.g. line 59)
# run on MacOS 14.6.1 

########################### ReadME end ########################### 

##### begin TEMPLATE##### 
########################### Figure ZZZZ ###########################
########################### description
########################### load dependencies
########################### custom functions
########################### load data
########################### data wrangling
########################### quick visualization
########################### analyses/modeling
########################### saving data
########################### saving figures
########################### clean up
##### end TEMPLATE #####


########################### prelims ########################### 

invisible(rm(list = ls()))
invisible(gc())
setwd("sangs_2025/")

########################### load general dependencies ########################### 
source("code/load_gen_dependencies.R")
rm(package.check)
########################### custom Functions 
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

########################### Figure 3c - 3g   ###########################
########################### mouse biodistribution
########################### load data
biodist_mouse_oc <- read_excel("processed_data/fig_3/raw_organ_biodistribution_mouse.xlsx", 
                               na = c("NA", "QNS"))
########################### data wrangling
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

###########################  visualization
figure_3c<-biodist_mouse_oc %>% 
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

figure_3g<-biodist_mouse_oc %>% 
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

figure_3d<-biodist_mouse_oc %>% 
  filter(tissue == "liver" |
           tissue == "kidney" |
           tissue == "spleen" |
           tissue == "heart" |
           tissue == "lungs" |
           tissue == "gi_healthy") %>%
  filter(ng_Control != "control") %>%
  filter(tumorMetNorm != "mets") %>%
   
  ggbarplot( x = "tissue", y = "fold_change_within_animal",width = 0.3, 
             add = c("mean_se", "dotplot"),
             color = "time_hrs", 
             fill = "time_hrs",
             palette = "Blues",
             add.params = list(size = .4),
             position = position_dodge(0.4),
             # facet.by = "ng_Control",
             ylim = c(0, 100)
  )+
  theme(legend.position = "none")
########################### saving figures
ggsave(figure_3c, file = "figure_3c.pdf", width = 3, height = 4, units = "in", path = "figures/fig_3")
ggsave(figure_3d, file = "figure_3d.pdf", width = 8, height = 4, units = "in", path = "figures/fig_3")
ggsave(figure_3g, file = "figure_3g.pdf", width = 3, height = 4, units = "in", path = "figures/fig_3")
########################### clean up
rm(figure_3c, figure_3d, figure_3g, biodist_mouse_oc, biodist_mouse_oc_mean)


