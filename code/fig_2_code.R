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

########################### Figure 2a-c & Supp Fig 7   ###########################
########################### description cell internalization and mechanism of entry
########################### load data
heyCell_allData <- read_excel("data/fig_2/heyCell_allData.xlsx", 
                              na = "NA")

########################### data wrangling
heyCell_allData$ngConcen_mgML <- factor(heyCell_allData$ngConcen_mgML, levels=c('3', '1.5', '0.75', '0.375'))
###########################  visualization

## sodium-azide-dependent internalization
heyCell_ATPdependent<-heyCell_allData %>% 
  filter(inhibitor == 'Sodium Azide') %>%
  # filter(objective == '63x') %>%
  filter(inhibitorConcen == 0 & mean > 2 | inhibitorConcen == 50) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")


heyCell_LatA<-heyCell_allData %>% 
  filter(inhibitor == 'LatA') %>%
  filter(objective == '20x') %>%
  filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 1|  inhibitorConcen == 2) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")


heyCell_cytoD<-heyCell_allData %>% 
  filter(inhibitor == 'cytoD') %>%
  filter(objective == '20x') %>%
  filter(intDen < 40000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 2|  inhibitorConcen == 4) %>%
  ggbarplot( x = "inhibitorConcen", y = "intDen",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")


heyCell_Cpz<-heyCell_allData %>% 
  filter(inhibitor == 'Cpz') %>%
  filter(objective == '20x') %>%
  mutate_at(vars("intDen"), funs(./1000)) %>%
    ggbarplot( x = "inhibitorConcen", y = "intDen",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")

heyCell_MBcd<-heyCell_allData %>% 
  filter(inhibitor == 'MBcd') %>%
  filter(objective == '20x') %>%
  filter(intDen < 30000) %>%
  # filter(inhibitorConcen == 0 | inhibitorConcen == 0.5| inhibitorConcen == 1) %>%
  ggbarplot( x = "inhibitorConcen", y = "mean",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "inhibitorConcen",
             palette = "Blues",
             width = 0.41,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")

########################### saving figures
ggsave(heyCell_ATPdependent, file = "heyCell_ATPdependent_2a.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_Cpz, file = "heyCell_Cpz_2b.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_cytoD, file = "heyCell_cytoD_2c.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_MBcd, file = "heyCell_MBcd_supp_fig_5.pdf", width = 1.5, height = 3, units = "in", path = "figures/supp_figs/")
ggsave(heyCell_LatA, file = "heyCell_LatA__supp_fig_5.pdf", width = 1.5, height = 3, units = "in", path = "figures/supp_figs/")
########################### clean up
rm(heyCell_ATPdependent, heyCell_Cpz, heyCell_cytoD, heyCell_LatA, heyCell_MBcd, heyCell_allData)


########################### Figure 2g  ###########################
########################### description: siRNA release and colocalization of siRNA and endosomes

########################### load data
sirnaEndosomColoc_allData <- read_excel("data/fig_2/heyCellsiRNAendosomeCOLOC.xlsx", 
                                        na = "NA")
########################### data wrangling
sirnaEndosomColoc_allData$cell <- factor(sirnaEndosomColoc_allData$cell)
sirnaEndosomColoc_allData$metric <- factor(sirnaEndosomColoc_allData$metric)
########################### visualization
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
########################### saving figures
ggsave(siRNA_endosome_Fig, file = "siRNA_endosome_Fig_2g.pdf", width = 5, height = 3, units = "in", path = "figures/fig_2")
########################### clean up
rm(siRNA_endosome_Fig, sirnaEndosomColoc_allData)




########################### Figure 2h siRNA release ###########################

########################### load data
heyCellsiRNAsangsCOLOC_allData <- read_excel("data/fig_2/heyCellsiRNAsangsCOLOC.xlsx", 
                                             na = "NA")
########################### data wrangling
heyCellsiRNAsangsCOLOC_allData$cell <- factor(heyCellsiRNAsangsCOLOC_allData$cell)
heyCellsiRNAsangsCOLOC_allData$metric <- factor(heyCellsiRNAsangsCOLOC_allData$metric)
heyCellsiRNAsangsCOLOC_allData$stage <- factor(heyCellsiRNAsangsCOLOC_allData$stage)
########################### visualization
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
########################### clean up
rm(siRNA_endosome_Fig, heyCellsiRNAsangsCOLOC_allData, siRNA_sang_Fig)



########################### Figure 2i RNA knockdown in vitro ###########################

########################### load data
rtPCR_data <- read_excel("data/fig_2/rtPCR_data.xlsx", 
                         na = "NA")
########################### visualization
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
########################### clean up
rm(inVitroSANGexpression, rtPCR_data)




