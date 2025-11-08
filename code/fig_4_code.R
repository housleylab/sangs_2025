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

########################### Figure 4c  ###########################
########################### rat biodistribution
########################### load data
biodist_rat_crc <- read_excel("processed_data/fig_4/raw_organ_biodistribution_rat.xlsx", 
                               na = c("NA", "QNS"))
########################### data wrangling
biodist_rat_crc <- biodist_rat_crc %>% filter(use == "yes")
biodist_rat_crc$scaledTotalRadiance <- scale_values(biodist_rat_crc$total_radiance_efficency)
biodist_rat_crc$tissue <- as.factor(biodist_rat_crc$tissue)
biodist_rat_crc$time_hrs <- factor(biodist_rat_crc$time_hrs, levels=c('control', '4', '24', '48', '288'))
biodist_rat_crc_mean<-summarySE(biodist_rat_crc, measurevar="fold_change_within_animal", groupvars=c("tissue","time_hrs"))
###########################  visualization
figure4c<- biodist_rat_crc %>% 
  filter(tissue != "healthy_colon") %>%
  filter(time_hrs != "control") %>%
  ggbarplot(  x = "tissue", y = "fold_change_within_animal", 
          color = "time_hrs", 
          fill = "time_hrs",
          add = c("mean_se", "jitter"), palette = "Blues",
          # add.params = list(fill = "black"),
          position = position_dodge(.8)) +
  scale_color_manual(values = c(V = "#000000", Z = "#000000"))


########################### saving figures
ggsave(figure4c, file = "figure4c.pdf", width = 6, height = 6, units = "in", path = "figures/fig_4/")
ggsave(figure4f, file = "figure4f.pdf", width = 6, height = 6, units = "in", path = "figures/fig_4/")
########################### clean up
rm(figure4c, figure4f, biodist_rat_crc_mean, biodist_rat_crc)



########################### Figure 4f   ###########################
########################### rat biodistribution
########################### load data
biodist_rat_lung <- read_excel("processed_data/fig_4/lung_biodist_database.xlsx", 
                              na = c("NA", "QNS"))
########################### data wrangling
biodist_rat_lung$condition <- as.factor(biodist_rat_lung$condition)
biodist_rat_lung$channel <- as.factor(biodist_rat_lung$channel)
biodist_rat_lung$time_factor <- factor(biodist_rat_lung$time_factor, levels=c('control','4', '24', '48'))
biodist_rat_lung$scaledTotalRadiance <- scale_values(biodist_rat_lung$total_radiance_efficency)

###########################  visualization
figure4f<- 
  biodist_rat_lung %>% 
  filter(channel != "mda_luc2") %>%
  filter(channel != "mda_gfp") %>%
  ggbarplot( x = "time_factor", y = "fold_change_flux", 
             color = "time_factor", 
             fill = "time_factor",
             add = c("mean_se", "jitter"), palette = "Blues",
             # add.params = list(fill = "black"),
             position = position_dodge(.8)) +
  scale_color_manual(values = c(V = "#000000", Z = "#000000"))

########################### saving figures
ggsave(figure4c, file = "figure4c.pdf", width = 6, height = 6, units = "in", path = "figures/fig_4/")
ggsave(figure4f, file = "figure4f.pdf", width = 6, height = 6, units = "in", path = "figures/fig_4/")
########################### clean up
rm(figure4c, figure4f, biodist_rat_crc_mean, biodist_rat_crc, biodist_rat_lung)
