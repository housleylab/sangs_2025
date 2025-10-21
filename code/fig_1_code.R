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
rm(heyCell_allData, heyCell_dose_dependentInternalization)













