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



########################### Figure 5a and Figure 5b  ###########################

########################### load data
rtPCR_data <- read_excel("data/fig_5/rtPCR_data.xlsx", 
                         na = "NA")
########################### data wrangling
rtPCR_data$dose <- as.factor(rtPCR_data$dose)
########################### visualization
########################### visualization
figure_5a<-rtPCR_data %>% 
  filter(experiment == 'in_vivo') %>%
  filter(species == 'mouse') %>%
  ggbarplot( x = "group", y = "fold_change", 
             color = "dose", 
             fill = "dose",
             add = c("mean_se", "jitter"), 
             palette = "Blues",
             position = position_dodge(.8),
             facet.by = "treatment") +
  scale_color_manual(values = c(V = "#000000", Z = "#000000"))

figure_5b<-rtPCR_data %>% 
  filter(experiment == 'in_vivo') %>%
  filter(species == 'rat') %>%
  ggbarplot( x = "group", y = "fold_change", 
             color = "dose", 
             fill = "dose",
             add = c("mean_se", "jitter"), 
             palette = "Blues",
             position = position_dodge(.8),
             facet.by = "treatment") +
  scale_color_manual(values = c(V = "#000000", Z = "#000000"))




########################### analyses/modeling
rtPCR_data %>% 
  filter(experiment == 'in_vivo') %>%
  group_by(group)%>%
  summarise(mean= mean(fold_change),
            sd = sd(fold_change))

cohensD = (1-0.540)/(sqrt((0.116^2+0.124^2)/2))
########################### saving figures
ggsave(figure_5a, file = "figure_5a.pdf", width = 6, height = 4, units = "in", path = "figures/fig_5")
ggsave(figure_5b, file = "figure_5b.pdf", width = 6, height = 4, units = "in", path = "figures/fig_5")

########################### clean up
rm(inVitroSANGexpression, rtPCR_data)




