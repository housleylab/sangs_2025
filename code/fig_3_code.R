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
########################### Clean up
##### end TEMPLATE #####


########################### prelims ########################### 

invisible(rm(list = ls()))
invisible(gc())
setwd("~/GaTech Dropbox/CoS/BioSci/BioSci-Housley_Lab/04-papers/nature_comm/sangs_2025/")

########################### load general dependencies ########################### 
source("code/load_gen_dependencies.R")
rm(package.check)
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


########################### Figure 2a-c & Supp Fig 7   ###########################
########################### description cell internalization and mechanism of entry
########################### load data
biodist_mouse_oc <- read_excel("data/fig_3/raw_organ_biodistribution_mouse.xlsx", 
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


# 2
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
             ylim = c(0, 100)
  )+
  theme(legend.position = "none")



########################### saving figures
ggsave(figure_3c, file = "figure_3c.pdf", width = 3, height = 4, units = "in", path = "figures/fig_3")
ggsave(figure_3d, file = "figure_3d.pdf", width = 8, height = 4, units = "in", path = "figures/fig_3")
ggsave(figure_3g, file = "figure_3g.pdf", width = 3, height = 4, units = "in", path = "figures/fig_3")
ggsave(heyCell_cytoD, file = "heyCell_cytoD_2c.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_MBcd, file = "heyCell_MBcd_supp_fig_5.pdf", width = 1.5, height = 3, units = "in", path = "figures/supp_figs/")
ggsave(heyCell_LatA, file = "heyCell_LatA__supp_fig_5.pdf", width = 1.5, height = 3, units = "in", path = "figures/supp_figs/")
########################### Clean up
rm(heyCell_ATPdependent, heyCell_Cpz, heyCell_cytoD, heyCell_LatA, heyCell_MBcd, heyCell_allData)


########################### Figure 2g  ###########################
########################### description: siRNA release and colocalization of siRNA and endosomes

########################### load data
sirnaEndosomColoc_allData <- read_excel("data/fig_2/heyCellsiRNAendosomeCOLOC.xlsx", 
                                        na = "NA")
########################### data wrangling
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
ggsave(siRNA_endosome_Fig, file = "siRNA_endosome_Fig_2g.pdf", width = 5, height = 3, units = "in", path = "figures/fig_2")
########################### Clean up
rm(siRNA_endosome_Fig, sirnaEndosomColoc_allData)




########################### Figure 2h siRNA release ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
heyCellsiRNAsangsCOLOC_allData <- read_excel("data/fig_2/heyCellsiRNAsangsCOLOC.xlsx", 
                                             na = "NA")

########################### Data Wrangling
heyCellsiRNAsangsCOLOC_allData$cell <- factor(heyCellsiRNAsangsCOLOC_allData$cell)
heyCellsiRNAsangsCOLOC_allData$metric <- factor(heyCellsiRNAsangsCOLOC_allData$metric)
heyCellsiRNAsangsCOLOC_allData$stage <- factor(heyCellsiRNAsangsCOLOC_allData$stage)

########################### quick visualization
siRNA_sang_Fig<- heyCellsiRNAsangsCOLOC_allData %>% 
  # filter(metric != "pearsonCorr") %>%
  filter(metric == "ColocCoeff1") %>%
  ggbarplot( x = "stage", y = "value",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "stage",
             palette = "Blues",
             width = 0.3,
             ylim = c(0, 1.3),
             add.params = list(size = .4),
             position = position_dodge(0.4)
             
  )
########################### saving figures
ggsave(siRNA_sang_Fig, file = "siRNA_sang_Fig_2h.pdf", width = 3, height = 4, units = "in", path = "figures/fig_2")
########################### Clean up
rm(siRNA_endosome_Fig, heyCellsiRNAsangsCOLOC_allData, siRNA_sang_Fig)



########################### Figure 2i RNA knockdown in vitro ###########################

########################### load data
rtPCR_data <- read_excel("data/fig_2/rtPCR_data.xlsx", 
                         na = "NA")
########################### quick visualization
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
  )
########################### analyses/modeling
rtPCR_data %>% 
  filter(experiment == 'in_vivo') %>%
  group_by(group)%>%
  summarise(mean= mean(fold_change),
            sd = sd(fold_change))

cohensD = (1-0.540)/(sqrt((0.116^2+0.124^2)/2))
########################### saving figures
ggsave(inVitroSANGexpression, file = "inVitroSANGexpression.pdf", width = 4, height = 4, units = "in", path = "figures/fig_2")
########################### Clean up
rm(inVitroSANGexpression, rtPCR_data)




