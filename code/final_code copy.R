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


########################### Figure 1d  ###########################
########################### description
# dls size distribution over batches
########################### Load Data
dls_ng_size <- read_excel("data/fig_1/dls_ng_size_batches_shell_core.xlsx",  na = "NA") 

########################### data wrangling
dls_ng_size$batch <- factor(dls_ng_size$batch)
dls_ng_size$index <- 'dls'
########################### visualization
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
            sdDiameter = sd(Peak2_rad_nm*2)
            )
########################### saving figures
ggsave(dls_batch_size_fig, file = "dls_batch_size_fig_1d.pdf", width = 6, height = 4, units = "in", path = "figures/fig_1/")
########################### Clean up
rm(dls_ng_size, dls_batch_size_fig)



########################### Figure 1e  ###########################
########################### description
# dls size distribution over time and temps
########################### Load Data
dls_ng_size_overTime <- read_excel("data/fig_1/dls_ng_size_overTime.xlsx",  na = "NA") 
########################### Data Wrangling
dls_ng_size_overTime$temp <- factor(dls_ng_size_overTime$temp)
dls_ng_size_overTime$loadingStatus <- factor(dls_ng_size_overTime$loadingStatus)
########################### visualization
dls_ng_size_overTime_fig<- ggplot(dls_ng_size_overTime, aes(x = time, y = percentChange, group = loadingStatus, colour = temp)) +
  theme_bw() +
  # geom_line() +
  geom_line(aes(linetype = loadingStatus))+
  geom_point()+
  facet_grid(. ~temp, scale = "free_y")+
  ylim(-25,25)+
  theme_classic()+
  scale_color_brewer(palette = "Blues")
########################### analyses/modeling
dls_ng_size_overTime %>% 
  group_by(loadingStatus, time) %>%
  summarise(n= n()
  )
########################### saving figures
ggsave(dls_ng_size_overTime_fig, file = "dls_ng_size_overTime_fig_1e.pdf", width = 6, height = 4, units = "in", path = "figures/fig_1/")
########################### Clean up
rm(dls_ng_size_overTime, dls_ng_size_overTime_fig)


########################### Figure 1f ###########################
########################### description
# viability in vitro 
########################### load data
viability_data <- read_excel("data/fig_1/viability.xlsx", 
                             na = "NA")
########################### data wrangling
viability_data$treatment <- as.factor(viability_data$treatment)
viability_data$group <- as.factor(viability_data$group)
viability_data$sample <- as.factor(viability_data$sample)
viability_data$parameter <- as.factor(viability_data$parameter)
########################### visualization
viability_fig <-viability_data %>% 
  filter(parameter == 'viability') %>%
  ggbarplot( x = "group", y = "values",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "group",
             palette = "Blues",
             width = 0.3,
             add.params = list(size = .4),
             position = position_dodge(0.4),
             facet.by = "parameter"
             
  )
# viability_fig<- facet(viability_fig, facet.by = "parameter", scales = "free_y")
########################### analyses/modeling
viability_data %>% 
  group_by(group, parameter) %>%
  summarise(n= n()
  )
########################### saving figures
ggsave(viability_fig, file = "viability_fig_1f.pdf", width = 4, height = 4, units = "in", path = "figures/fig_1/")
########################### Clean up
rm(viability_data, viability_fig)



########################### Figure 1h  ###########################
########################### description
# cell internalization
########################### load data
heyCell_allData <- read_excel("data/fig_1/heyCell_allData.xlsx", 
                              na = "NA")
########################### data wrangling
heyCell_allData$ngConcen_mgML <- factor(heyCell_allData$ngConcen_mgML, levels=c('3', '1.5', '0.75', '0.375'))
########################### visualization
heyCell_dose_dependentInternalization<-heyCell_allData %>% 
  filter(inhibitor == 'na') %>%
  filter(objective == '20x') %>%
  ggbarplot( x = "ngConcen_mgML", y = "mean",
             add = c("mean_se", "jitter"),
             color = "black",
             fill = "ngConcen_mgML",
             palette = "Blues",
             width = 0.8,
             add.params = list(size = .4),
             position = position_dodge(0.4),
  )+
  theme(legend.position = "none")
  # stat_compare_means(label.y = 7) +                                         # Global p-value
  # stat_compare_means(ref.group = "3", label = "p.signif",
  #                    label.y = c(4.5, 3.5, 2.5))
########################### analyses/modeling
heyCell_allData %>% 
  filter(inhibitor == 'na') %>%
  filter(objective == '20x') %>%
  group_by(ngConcen_mgML) %>%
  summarise(n=n())
########################### saving data
########################### saving figures
ggsave(heyCell_dose_dependentInternalization, file = "heyCell_dose_dependentInternalization_1h.pdf", width = 3, height = 5, units = "in", path = "figures/fig_1/")
########################### Clean up




########################### Figure 2a-c & Supp Fig 7   ###########################
########################### description cell internalization and mechanism of entry
########################### load data
heyCell_allData <- read_excel("data/fig_2/heyCell_allData.xlsx", 
                              na = "NA")

########################### Data Wrangling
heyCell_allData$ngConcen_mgML <- factor(heyCell_allData$ngConcen_mgML, levels=c('3', '1.5', '0.75', '0.375'))
########################### quick visualization

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
########################### analyses/modeling
heyCell_allData %>% 
  filter(inhibitor == 'Sodium Azide') %>%
  # filter(objective == '63x') %>%
  filter(inhibitorConcen == 0 & mean > 2 | inhibitorConcen == 50) %>%
  group_by(inhibitorConcen) %>%
  summarise(n=n())
  
heyCell_allData %>% 
  filter(inhibitor == 'Cpz') %>%
  # filter(objective == '63x') %>%
  # filter(inhibitorConcen == 0 & mean > 2 | inhibitorConcen == 50) %>%
  group_by(inhibitorConcen) %>%
  summarise(n=n())

heyCell_allData %>% 
  filter(inhibitor == 'cytoD') %>%
  # filter(objective == '63x') %>%
  # filter(inhibitorConcen == 0 & mean > 2 | inhibitorConcen == 50) %>%
  group_by(inhibitorConcen) %>%
  summarise(n=n())
########################### saving figures
ggsave(heyCell_ATPdependent, file = "heyCell_ATPdependent_2a.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_Cpz, file = "heyCell_Cpz_2b.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_cytoD, file = "heyCell_cytoD_2c.pdf", width = 1.5, height = 3, units = "in", path = "figures/fig_2")
ggsave(heyCell_MBcd, file = "heyCell_MBcd_supp_fig_5.pdf", width = 1.5, height = 3, units = "in", path = "figures/supp_figs/")
ggsave(heyCell_LatA, file = "heyCell_LatA__supp_fig_5.pdf", width = 1.5, height = 3, units = "in", path = "figures/supp_figs/")
########################### Clean up










########################### Figure 2g  ###########################
########################### description: siRNA release and colocalization of siRNA and endosomes

########################### Load Data
sirnaEndosomColoc_allData <- read_excel("data/figure_5/heyCellsiRNAendosomeCOLOC.xlsx", 
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
  # coord_flip()+
  theme(legend.position = "none")


########################### analyses/modeling
sirnaEndosomColoc_allData %>%
  summarise(meanColoc = mean(pearsonCorr),
            sdColoc = sd(pearsonCorr))
########################### saving data
########################### saving figures
setwd(file.path(mainDir,figDir,figFolder))


cairo_pdf("siRNA_endosome_Fig.pdf", width = 5, height = 3)
siRNA_endosome_Fig+coord_flip()
invisible(suppressMessages(suppressWarnings(dev.off())))



########################### Clean up











########################### Figure 5b ###########################
########################### description
## TEM modeling

########################### load dependencies
########################### Custom Functions
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
########################### Load Data
tem_allData <- read_excel("data/figure_6/tem_ng_data_raw.xlsx", 
                          na = "NA")

########################### Data Wrangling
tem_allData$concentrationFactor <- factor(tem_allData$concentrationFactor, levels=c('ng_0.46', 'ng_0.93', 'ng_1.875', 'ng_3.75', 'ng_7.5', 'ng_30'))

# tem_allDataScaled<-tem_allData %>% mutate_at(c("area"), ~(scale(.) %>% as.vector))
temSummary <- summarySE(tem_allData, measurevar="area", groupvars=c("concentrationFactor", "concentrationNum"))


pd <- position_dodge(0.1) # move them .05 to the left and right

########################### quick visualization

fig_6g<-temSummary %>%
  filter(concentrationNum != 'NA') %>%
  ggplot( aes(x=concentrationNum, y=area, group=1)) + 
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)+
  geom_line(data = temSummary, aes(x = concentrationNum, y = (area)**0+3136, group = 1), 
            inherit.aes = FALSE,linetype = "dashed")+
  geom_line(data = temSummary, aes(x = concentrationNum, y = c(66.3044, 134.0502,272.4246, 540.525, 1081.05, 4324.2)*100, group = 1),
            inherit.aes = FALSE, linetype = "dashed", color="red")+
  geom_line(data = temSummary, aes(x = concentrationNum, y = (N)*100, group = 1), 
            inherit.aes = FALSE,linetype = "solid", color="red")+
  theme_classic()+
  scale_y_continuous("NG particle area", trans='log2',   sec.axis = sec_axis(~.x * .01, name = "NG particle #"))+
  scale_x_continuous(trans='log2')+
  theme(axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red"))




########################### analyses/modeling
setwd("figures/")
cairo_ps("fig_6g.eps")
fig_6g
invisible(suppressMessages(suppressWarnings(dev.off())))

########################### saving data
########################### saving figures
########################### Clean up



















