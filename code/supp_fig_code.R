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
# v1.2- original Oct 21, 2025. - publish to public github

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

########################### Supp Fig 3   ###########################
########################### 
########################### load data
heyCell_allData <- read_excel("processed_data/supp_figs/time_dependent_internal.xlsx", 
                              na = "NA")
########################### data wrangling
###########################  visualization
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

########################### saving figures
ggsave(heyCell_time_dependentInternalization, file = "heyCell_time_dependentInternalization.pdf", width = 6, height = 12, units = "in", path = "figures/supp_figs/")

########################### clean up
invisible(rm(list = ls()))
invisible(gc())

