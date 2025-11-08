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
########################### Figure 6   ###########################
########################### description: in vivo coloc
########################### load data
Results_colocalization <- read_csv("processed_data/fig_6/Results_colocalization.csv", 
                                   na = c("NA", "QNS", "NaN"))
###########################  visualization
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
########################### saving figures
ggsave(sang_colocalization_size_Fig, file = "sang_colocalization_size_Fig.pdf", width = 6, height = 17, units = "in", path = "figures/fig_6/")
########################### clean up
rm(sang_colocalization_size_Fig)


########################### Figure 6  ###########################
########################### description: DLS
########################### load data
dls_ng_zetasizer_concentration<- read_excel("processed_data/fig_6/dls_zetaSizer_concentration.xlsx",  na = "NA") 
########################### data wrangling
dls_ng_zetasizer_concentration$acquisition <- factor(dls_ng_zetasizer_concentration$acquisition)
dls_ng_zetasizer_concentration$concentration <- factor(dls_ng_zetasizer_concentration$concentration)
########################### visualization
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
########################### saving figures
ggsave(dls_ng_zetasizer_fig, file = "dls_ng_zetasizer_fig.pdf", width = 6, height = 4, units = "in", path = "figures/fig_6/")
########################### clean up
rm(siRNA_endosome_Fig, sirnaEndosomColoc_allData)



########################### Figure 6  ###########################
########################### description: dosy
########################### load data
dosy_data<- read_excel("processed_data/fig_6/dosy.xlsx",  na = "NA") 
########################### data wrangling
dosy_data$metric <- factor(dosy_data$metric)
########################### visualization
dls_ng_size_overTime_fig<- ggplot(dosy_data, aes(x = nanomolar, y = diffusion_coefficient, colour = metric)) +
  theme_bw() +
  geom_line() +
  geom_point()+
  facet_grid(. ~metric, scale = "free_y")+
  ylim(-25,25)+
  theme_classic()+
  scale_color_brewer(palette = "Blues")

dosy_fig<-dosy_data %>% 
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
########################### saving figures
ggsave(dosy_fig, file = "dosy_fig.pdf", width = 6, height = 4, units = "in", path = "figures/fig_6/")
########################### clean up
rm(siRNA_endosome_Fig, sirnaEndosomColoc_allData)



########################### Figure 6  ###########################
########################### description: tem 
########################### custom functions
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
########################### load data
tem_allData <- read_excel("processed_data/fig_6/tem_ng_data_raw.xlsx", 
                          na = "NA")
########################### data wrangling
tem_allData$concentrationFactor <- factor(tem_allData$concentrationFactor, levels=c('ng_0.46', 'ng_0.93', 'ng_1.875', 'ng_3.75', 'ng_7.5', 'ng_30'))
temSummary <- summarySE(tem_allData, measurevar="area", groupvars=c("concentrationFactor", "concentrationNum"))
pd <- position_dodge(0.1) # move them .05 to the left and right
########################### visualization
tem_fig<-temSummary %>%
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
########################### saving figures
ggsave(tem_fig, file = "tem_fig.pdf", width = 6, height = 4, units = "in", path = "figures/fig_6/")
########################### clean up
rm(siRNA_endosome_Fig, sirnaEndosomColoc_allData)

