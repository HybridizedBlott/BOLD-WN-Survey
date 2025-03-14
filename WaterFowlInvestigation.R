##*****************************
##Thomas Tekle - 1079669
##
##University of Guelph
## 
##2024-09-24
##
##******************************

#Mosquito are a major vector/reservoir of the west nile virus, and so plenty of efforts have been made in order to prevent their propagation and migration, therefore limiting the spread of disease and mitigating the risk of an outbreak #occurring. However, west nile virus has been observed to have infected various types of birds including the waterfowl, a highly migratory bird, which presents another taxonomic group pf interest to monitor, as to prevent the further spread of the zoonotic disease. Although there have been studies which address the role waterfowls #play in the spread of WNV, none have used comprehensive databases to pinpoint and #correlate the spread of virus to the migratory behavior of the waterfowl. This #project explores the global distribution of the waterfowl and its correlation to the #distribution of WNV outbreaks globally. This may also serve as a proof of concept #for the use of BOLD as a species monitoring platform.....

#installing required packages

packageList <- c("tidyverse","viridisLite","dplyr","viridis","vegan","ggplot2","remotes","pheatmap")

#install.packages(packageList)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("sangerseqR")

#remotes::install_github("ropensci/bold")

#initializing required libraries
for (i in 1:length(packageList)){
  lapply(packageList[i],library, character.only = TRUE)
}
library(bold)

#Installing waterfowl data frame from BOLD database, note cleanData argument will replace missing cells with N/A.

dfWaterfowl <- bold_specimens(taxon = 'Anseriformes', cleanData = TRUE)

#Subsetting data frame and filtering out NA's. The filtering of this data included all countries which provided samples in the BOLD system. All data frames using all countries will be labelled with either a '.full' or '.F' to denote a full data frame.

dfWaterfowl.filter.full <- dfWaterfowl %>% 
  select(c(country,bin_uri)) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(bin_uri))

dfWF.filter.by.country.F <- dfWaterfowl.filter.full %>% 
  group_by(country,bin_uri) %>% 
  summarise(count = n()) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(bin_uri))

dfWF.spread.by.country.F <- dfWF.filter.by.country.F %>% 
  pivot_wider(names_from = bin_uri, values_from = count) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'country')

dfWF.spread.by.country.F[is.na(dfWF.spread.by.country.F)] <- 0

AccumCurveWF.F <- specaccum(dfWF.spread.by.country.F)

#Second round of filtering, However, Only 7 countries with relatively high cases of WNV were used to narrow the focus in the heat map. An additional 6 countries with no cases of WNV were included to act as a control. The data frames formed using the selected countries were denoted with '.relevant' or '.R'.

dfWaterfowl.filter.relevant <- dfWaterfowl %>% 
  select(country,bin_uri) %>%
  filter(country %in% c("Italy","Canada","Mexico","India","Germany","Russia","Sweden","Norway","Japan","United Kingdom","New Zealand")) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(bin_uri))

dfWF.filter.by.country.R <- dfWaterfowl.filter.relevant %>% 
  group_by(country,bin_uri) %>% 
  summarise(count = n()) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(bin_uri))

dfWF.spread.by.country.R <- dfWF.filter.by.country.R %>% 
  pivot_wider(names_from = bin_uri, values_from = count) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'country')

dfWF.spread.by.country.R[is.na(dfWF.spread.by.country.R)] <- 0

AccumCurveWF.R <- specaccum(dfWF.spread.by.country.R)

diss.vegdist.output <- vegdist(dfWF.spread.by.country.R)
diss.matrix <- as.matrix(diss.vegdist.output)
diss.df <- as.data.frame(diss.matrix)

#Plotting a bar graph for the number of BINs associated in each country

#Bar graph

ggplot(data = dfWF.filter.by.country, aes(x = country, y = count, fill = country)) +
  geom_bar(stat = "identity")+
  labs(title = "TITLE", y = "Total number of BINS per Country", x = 'Countries ')+
  coord_flip()+
  scale_x_discrete(limits = c(sort(unique(x = dfWF.filter.by.country$country), decreasing = T)))


#Sample Completeness Curve
plot(AccumCurveWF.F, main = "Sampling completeness of the taxanomic order Anseriformes globally", ylab = 'BIN Richness', xlab = 'Countries sampled',,, col = "blue")
plot(AccumCurveWF.R, main = "Sampling completeness of the taxanomic order Anseriformes globally", ylab = 'BIN Richness', xlab = 'Countries sampled',,, col = "red", add = TRUE)
legend(x = "bottomright", legend = c("Full Data", "Relevent Data"), fill=c("blue","red"))

#Heatmap
pheatmap(diss.df, main = "Heatmap comparing Anseriformes BIN Composition similarity of different countries")
