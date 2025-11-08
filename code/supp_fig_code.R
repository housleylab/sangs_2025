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
# setwd("sangs_2025/")

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


########################### Supp Fig 20 single dose tox porcine ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data

ng_tox_pig <- read_excel("processed_data/supp_figs/ng_tox_porcine.xlsx", 
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

########################### saving figures
ggsave(figure4_pig_serum_tox, file = "supp_fig_20_pig_serum_tox.pdf", width = 5, height = 4, units = "in", path = "figures/supp_figs/")

########################### clean up
invisible(rm(list = ls()))
invisible(gc())




########################### Supp Fig 19 single dose tolerability cd1 ###########################
########################### Load Data
ng_tox_cd1 <- read_excel("processed_data/supp_figs/ng_tox_cd_1.xlsx", 
                         na = c("NA", "QNS"))
########################### Data Wrangling
ng_tox_cd1$collection_time_hr <- as.factor(ng_tox_cd1$collection_time_hr)
ng_tox_cd1$sample_number <- as.factor(ng_tox_cd1$sample_number)
ng_tox_cd1$group <- as.factor(ng_tox_cd1$group)
ng_tox_cd1$animal_num_glp <- as.factor(ng_tox_cd1$animal_num_glp)
#### 
ng_tox_cd1_df<- ng_tox_cd1
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
########################### saving figures
ggsave(figure4_mouse_serum_tox, file = "supp_fig_19_mouse_serum_tox.pdf", width = 5, height = 4, units = "in", path = "figures/supp_figs/")

########################### clean up
invisible(rm(list = ls()))
invisible(gc())




########################### Supp Fig 20 MTD rat ###########################

########################### load dependencies
########################### Custom Functions
########################### Load Data
ng_tox_mtd_rat <- read_excel("processed_data/supp_figs/ng_tox_mtd_rat.xlsx", 
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
########################### saving figures
ggsave(figure4_rat_serum_tox, file = "supp_fig_20_rat_serum_tox.pdf", width = 5, height = 4, units = "in", path = "figures/supp_figs/")

########################### clean up
invisible(rm(list = ls()))
invisible(gc())




########################### Supp Fig 20 repeat dose tolerability rat###########################
########################### load dependencies
########################### Custom Functions
########################### Load Data
tox_rat <- read_excel("processed_data/supp_figs/rat_repeatDose_weights.xlsx", 
                      na = c("NA", "QNS"))
########################### Data Wrangling
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
########################### saving figures
ggsave(figure4_rat_percentChange, file = "supp_fig_20_rat_percentChange.pdf", width = 6, height = 4, units = "in", path = "figures/supp_figs/")

########################### clean up
invisible(rm(list = ls()))
invisible(gc())









########################### Supp Fig 19 toxicity cd1 ###########################


########################### load dependencies
########################### Custom Functions
########################### Load Data
tox_cd1 <- read_excel("processed_data/supp_figs/cd_1_weights.xlsx", 
                      na = c("NA", "QNS"))
########################### Data Wrangling
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
########################### saving figures
ggsave(figure4_cd1_percentChange, file = "supp_fig_19_cd1_percentChange.pdf", width = 6, height = 4, units = "in", path = "figures/supp_figs/")
ggsave(figure4_cd1_mean, file = "supp_fig_19_cd1_mean.pdf", width = 6, height = 4, units = "in", path = "figures/supp_figs/")

########################### clean up
invisible(rm(list = ls()))
invisible(gc())






