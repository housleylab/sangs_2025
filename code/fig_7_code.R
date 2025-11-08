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
source("code/load_gen_dependencies_corona.R")
rm(package.check)
cols <- colorRampPalette(brewer.pal(10, "PiYG"))(256)
cols <- colorRampPalette(brewer.pal(9, "PRGn"))(256)
cols <- colorRampPalette(viridis(12))

########################### Figure 7   ###########################
########################### description: 
########################### load data
df <- read_excel("processed_data/fig_7/protein_corona_raw.xlsx")
########################### data wrangling
newnames <- c(fdrCombined = "Protein FDR Confidence: Combined", 
              numUniquePeptides = "# Unique Peptides",
              numPeptides = "# Peptides",
              lp01_s1 = "Abundance: F26: Sample, LP01",
              lp01_s2 = "Abundance: F27: Sample, LP01",
              lp01_s3 = "Abundance: F28: Sample, LP01",
              lp01_s4 = "Abundance: F29: Sample, LP01",
              lp01_s5 = "Abundance: F30: Sample, LP01",
              plasma_s1 = "Abundance: F31: Sample, PlasmaOnly", 
              plasma_s2 = "Abundance: F32: Sample, PlasmaOnly",
              plasma_s3 = "Abundance: F33: Sample, PlasmaOnly",
              plasma_s4 = "Abundance: F34: Sample, PlasmaOnly",
              plasma_s5 = "Abundance: F35: Sample, PlasmaOnly",
              plga_s1 = "Abundance: F36: Sample, PLGA",
              plga_s2 = "Abundance: F37: Sample, PLGA",
              plga_s3 = "Abundance: F38: Sample, PLGA",
              plga_s4 = "Abundance: F39: Sample, PLGA",
              plga_s5 = "Abundance: F40: Sample, PLGA",
              sang_s1 =  "Abundance: F41: Sample, SANG",
              sang_s2 =  "Abundance: F42: Sample, SANG",
              sang_s3 =  "Abundance: F43: Sample, SANG",
              sang_s4 =  "Abundance: F44: Sample, SANG",
              sang_s5 =  "Abundance: F45: Sample, SANG",
              vlp_s1 = "Abundance: F46: Sample, VLP", 
              vlp_s2 = "Abundance: F47: Sample, VLP",
              vlp_s3 = "Abundance: F48: Sample, VLP",
              vlp_s4 = "Abundance: F49: Sample, VLP",
              vlp_s5 = "Abundance: F50: Sample, VLP",
              lp01_mean = "Abundances (Grouped): LP01",
              plasma_mean = "Abundances (Grouped): PlasmaOnly",
              plga_mean = "Abundances (Grouped): PLGA",
              sang_mean = "Abundances (Grouped): SANG",
              vlp_mean = "Abundances (Grouped): VLP"
)

df<-dplyr::rename(df, all_of(newnames))
rm(newnames)

df<-df %>% separate_wider_delim(Description, " OS=Rattus norvegicus", names = c("description", "oldNames")) 
df<-df%>% arrange(Accession)

######### biological processes analysis
x <- strsplit(as.character(df$`Cellular Component`), ";")
l <- lengths(x)  ## R 3.3.0 onward
m <- max(l)
x <- t(sapply(x[as.logical(l)], function(a) c(a, rep("",m-length(a)))))
cellComp <- as.data.frame(x)
names(cellComp) <- paste0('cellComp_', seq_len(length(cellComp)))

x <- strsplit(as.character(df$`Biological Process`), ";")
l <- lengths(x)  ## R 3.3.0 onward
m <- max(l)
x <- t(sapply(x[as.logical(l)], function(a) c(a, rep("",m-length(a)))))
bioProcess <- as.data.frame(x)
names(bioProcess) <- paste0('bioProcess_', seq_len(length(bioProcess)))

x <- strsplit(as.character(df$`Molecular Function`), ";")
l <- lengths(x)  ## R 3.3.0 onward
m <- max(l)
x <- t(sapply(x[as.logical(l)], function(a) c(a, rep("",m-length(a)))))
molFun <- as.data.frame(x)
names(molFun) <- paste0('molFun_', seq_len(length(molFun)))

topProcessDF<- cbind(cellComp[1], bioProcess[1], molFun[1])

df_heatmap<-cbind(df,topProcessDF)
rm(x,molFun,cellComp,bioProcess)

######### data imputation

dfLong<- df_heatmap %>%
  select(Accession,lp01_s1,lp01_s2,lp01_s3,lp01_s4,lp01_s5,
         sang_s1,sang_s2,sang_s3,sang_s4,sang_s5,
         plga_s1,plga_s2,plga_s3,plga_s4,plga_s5,
         vlp_s1, vlp_s2, vlp_s3, vlp_s4, vlp_s5,
         plasma_s1,plasma_s2,plasma_s3,plasma_s4,plasma_s5, 
         cellComp_1,bioProcess_1, molFun_1) %>%   
  group_by(Accession) %>%
  gather(variable, value, 'lp01_s1':'plasma_s5', factor_key=TRUE)

experimentNames_df<-data.frame(str_split(dfLong$variable, "_", simplify=TRUE))
oldnames = c("X1","X2")
newnames = c("cohort","sample")
experimentNames_df<- experimentNames_df %>% rename_at(vars(oldnames), ~ newnames)

dfLong<-cbind(experimentNames_df,dfLong)
rm(experimentNames_df)

dfLong<-dfLong%>%
  group_by(cohort,Accession)%>%
  mutate(impute = ifelse(is.na(value), mean(value, na.rm = TRUE), value))

dfwide<-dfLong[,-1] %>% select(Accession,variable, impute) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(impute)
  ) %>% arrange(Accession)


dfwideMean<-dfwide %>%
  pivot_longer(cols=lp01_s1:plasma_s5, names_pattern = "(.*)_(.*)", names_to=c("particle","sample")) %>%
  dplyr::group_by(Accession, particle) %>%
  # mutate(scaledVal = scale_this(value)) %>%
  dplyr::summarise(mean_val=mean(value, na.rm = TRUE)) %>% 
  pivot_wider(
    names_from = particle,
    values_from = c(mean_val)
  ) %>% arrange(Accession)


dfwideMean<-cbind(dfwideMean[1],scale(dfwideMean[-1]))

### bind imputed data original
drop.cols <- c("Accession","lp01_s1","lp01_s2","lp01_s3","lp01_s4","lp01_s5",
               "sang_s1","sang_s2","sang_s3","sang_s4","sang_s5",
               "plga_s1","plga_s2","plga_s3","plga_s4","plga_s5",
               "vlp_s1", "vlp_s2", "vlp_s3", "vlp_s4", "vlp_s5",
               "plasma_s1","plasma_s2","plasma_s3","plasma_s4","plasma_s5")
df_heatmap<-df_heatmap %>% select(-one_of(drop.cols))


df_heatmap<-cbind(dfwideMean,dfwide,df_heatmap)











###########################  visualization
########################   venn Diagram
sang_sig<-df %>% 
  filter(sang_plasma_sig == 'yes' & `Abundance Ratio: (SANG) / (PlasmaOnly)` >1.2) %>%
  select(Accession, description, `Abundance Ratio: (SANG) / (PlasmaOnly)`)

lp01_sig<-df %>% 
  filter(lp01_plasma_sig == 'yes' & `Abundance Ratio: (LP01) / (PlasmaOnly)` >1.2) %>%
  select(Accession, description, `Abundance Ratio: (LP01) / (PlasmaOnly)`)

plga_sig<-df %>% 
  filter(plga_plasma_sig == 'yes' & `Abundance Ratio: (PLGA) / (PlasmaOnly)` >1.2) %>%
  select(Accession, description, `Abundance Ratio: (PLGA) / (PlasmaOnly)`)

vlp_sig<-df %>%
  filter(vlp_plasma_sig == 'yes' & `Abundance Ratio: (VLP) / (PlasmaOnly)` >1.2) %>%
  select(Accession, description, `Abundance Ratio: (VLP) / (PlasmaOnly)`)

vennDB<- list(SANG = sang_sig$Accession,
              # PLGA = plga_sig$Accession,
              LP01 = lp01_sig$Accession
              # VLP = vlp_sig$Accession
)

venn_fig<-ggvenn(
  vennDB, 
  fill_color = c("#a655aa", "#0073C2"),
  stroke_size = 0.5, set_name_size = 4
)



########################   heatmap

df_heatmap1 <- df_heatmap %>%
  filter( lp01_plasma_sig == 'yes' & `Abundance Ratio: (LP01) / (PlasmaOnly)` >1.5 |
            sang_plasma_sig == 'yes' & `Abundance Ratio: (SANG) / (PlasmaOnly)` >1.5 |
            plga_plasma_sig == 'yes' & `Abundance Ratio: (PLGA) / (PlasmaOnly)` >1.5|
            vlp_plasma_sig == 'yes' & `Abundance Ratio: (VLP) / (PlasmaOnly)` >1.5)%>%
  # tidyr::unite(rowname, description, Accession) %>%
  column_to_rownames(var="Accession...1") 

# generate the subsetted DF
df_heatmap2<-df_heatmap1 %>%
  drop_na(lp01_s1,lp01_s2,lp01_s3,lp01_s4,lp01_s5,
          sang_s1,sang_s2,sang_s3,sang_s4,sang_s5,
          # plga_s1,plga_s2,plga_s3,plga_s4,plga_s5,
          # vlp_s1, vlp_s2, vlp_s3, vlp_s4, vlp_s5,
          plasma_s1,plasma_s2,plasma_s3,plasma_s4,plasma_s5) %>%
  select(lp01_s1,lp01_s2,lp01_s3,lp01_s4,lp01_s5,
         sang_s1,sang_s2,sang_s3,sang_s4,sang_s5,
         # plga_s1,plga_s2,plga_s3,plga_s4,plga_s5,
         # vlp_s1, vlp_s2, vlp_s3, vlp_s4, vlp_s5,
         plasma_s1,plasma_s2,plasma_s3,plasma_s4,plasma_s5, 
         cellComp_1,bioProcess_1, molFun_1) 

### coerce to none if no BP/CC/MF process identified
df_heatmap2$molFun_1[is.na(df_heatmap2$molFun_1)] <- "a_none"
df_heatmap2$bioProcess_1[is.na(df_heatmap2$bioProcess_1)] <- "a_none"
df_heatmap2$cellComp_1[is.na(df_heatmap2$cellComp_1)] <- "a_none"

# generate the matrix
mat <- as.matrix(df_heatmap2[,1:(length(df_heatmap2)-3)])
scaled_mat = t(scale(t(mat)))
type <- gsub("_s\\d", "", colnames(mat))




# Create the heatmap annotation
col_fun = list(c("#868686", "#cd534c", "#efc000", "#0073c2", "#a655aa"))

ha = HeatmapAnnotation(
  df = data.frame(NP = type),
  annotation_height = unit(4, "mm"),
  col = list(NP = c("lp01" = "#0073c2",
                    "sang" = "#a655aa",
                    "plga" = "#efc000",
                    "vlp" = "#cd534c",
                    "plasma" = "#868686"
  )))


### find optimal k for kmeans = 8 based on elbow method
fviz_nbclust(scaled_mat, kmeans, method = "wss", nboot = 1000, k.max = 20)
# f1 = colorRamp2(seq(min(scaled_mat), max(scaled_mat), length = 3), c("blue", "#EEEEEE", "red"))
# f2 = colorRamp2(seq(min(scaled_mat), max(scaled_mat), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

set.seed(12) # must use for reproducibility
hm = Heatmap(scaled_mat, name = "expression", km = 8, top_annotation = ha,
             show_row_names = F, show_column_names = F) +
  Heatmap(df_heatmap2$molFun_1, name = "MF", width = unit(5, "mm"),
          col = brewer.pal(n = length(unique(df_heatmap2$molFun_1)), name = "Paired"))+
  Heatmap(df_heatmap2$cellComp_1, name = "CC", width = unit(5, "mm"),
          col = brewer.pal(n = length(unique(df_heatmap2$cellComp_1)), name = "Paired")) +
  Heatmap(df_heatmap2$bioProcess_1, name = "BP", width = unit(5, "mm"),
          col = brewer.pal(n = length(unique(df_heatmap2$bioProcess_1)), name = "Paired")
          # col = circlize::rand_color(length(unique(df_heatmap2$bioProcess_1)))
  )
### fig 7d
heatmap = draw(hm) # or Hm2236 = make_layout(Hm2236)
r.dend <- row_dend(heatmap)  #Extract row dendrogram
rcl.list <- row_order(heatmap)  #Extract clusters (output is a list)
# lapply(rcl.list, function(x ) length(x)) 
# rownames(scaled_mat)[ rcl.list[[1]] ]

# loop to extract genes for each cluster.
for (i in 1:length(row_order(heatmap))){
  if (i == 1) {
    clu <- t(t(row.names(scaled_mat[row_order(heatmap)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("Accession...1", "Cluster")
  } else {
    clu <- t(t(row.names(scaled_mat[row_order(heatmap)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}


### fig 7b
fig_densityHeatmap<-densityHeatmap(scaled_mat)



########################### analyses/modeling
unique_LP01<-anti_join(lp01_sig,sang_sig) %>% arrange(desc(`Abundance Ratio: (LP01) / (PlasmaOnly)`)) %>% dplyr::rename(abundanceRatio = `Abundance Ratio: (LP01) / (PlasmaOnly)`) %>% mutate(group = "lp01_sig")
unique_SANG <- anti_join(sang_sig,lp01_sig) %>% arrange(desc(`Abundance Ratio: (SANG) / (PlasmaOnly)`)) %>% dplyr::rename(abundanceRatio = `Abundance Ratio: (SANG) / (PlasmaOnly)`) %>% mutate(group = "sang_sig")
shared_SANG_LPO1<- inner_join(sang_sig,lp01_sig) %>% select(!`Abundance Ratio: (SANG) / (PlasmaOnly)`) %>%arrange(desc(`Abundance Ratio: (LP01) / (PlasmaOnly)`)) %>% dplyr::rename(abundanceRatio = `Abundance Ratio: (LP01) / (PlasmaOnly)`) %>% mutate(group = "shared_sig")



cluster1<-out %>% as.data.frame() %>% filter(Cluster == "cluster1") %>% select(`Accession...1`) %>% as.vector() 
# cluster1= c("P23562", "A0A8I6GER1", "P06765", "A0A0G2K4H7", "A0A8L2Q996") 
# cluster1= out[,1]
cluster1<-out %>% as.data.frame() %>% filter(Cluster == "cluster1") %>% select(`Accession...1`) %>% as.vector() 
c1<-df_heatmap %>% filter(`Accession...1` %in% cluster1$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster1")

cluster2<-out %>% as.data.frame() %>% filter(Cluster == "cluster2") %>% select(`Accession...1`) %>% as.vector() 
c2<-df_heatmap %>% filter(`Accession...1` %in% cluster2$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster2")

cluster3<-out %>% as.data.frame() %>% filter(Cluster == "cluster3") %>% select(`Accession...1`) %>% as.vector() 
c3<-df_heatmap %>% filter(`Accession...1` %in% cluster3$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster3")

cluster4<-out %>% as.data.frame() %>% filter(Cluster == "cluster4") %>% select(`Accession...1`) %>% as.vector() 
c4<-df_heatmap %>% filter(`Accession...1` %in% cluster4$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster4")

cluster5<-out %>% as.data.frame() %>% filter(Cluster == "cluster5") %>% select(`Accession...1`) %>% as.vector() 
c5<-df_heatmap %>% filter(`Accession...1` %in% cluster5$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster5")

cluster6<-out %>% as.data.frame() %>% filter(Cluster == "cluster6") %>% select(`Accession...1`) %>% as.vector() 
c6<-df_heatmap %>% filter(`Accession...1` %in% cluster6$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster6")

cluster7<-out %>% as.data.frame() %>% filter(Cluster == "cluster7") %>% select(`Accession...1`) %>% as.vector() 
c7<-df_heatmap %>% filter(`Accession...1` %in% cluster7$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster7")

cluster8<-out %>% as.data.frame() %>% filter(Cluster == "cluster8") %>% select(`Accession...1`) %>% as.vector() 
c8<-df_heatmap %>% filter(`Accession...1` %in% cluster8$`Accession...1`) %>% select(description) %>% mutate(cluster = "cluster8")



########################### saving figures
ggsave(venn_fig, file = "venn_fig.pdf", width = 10, height = 10, units = "in", path = "figures/fig_7/")

########################### saving data
write.csv(rbind(unique_LP01,unique_SANG,shared_SANG_LPO1),file = file.path("processed_data/fig_7/","vennDiagramTable"))
write.csv(rbind(c1,c2,c3,c4,c5,c6,c7,c8),file = file.path("processed_data/fig_7/","heatmapTable"))

########################### clean up
invisible(rm(list = ls()))
invisible(gc())
