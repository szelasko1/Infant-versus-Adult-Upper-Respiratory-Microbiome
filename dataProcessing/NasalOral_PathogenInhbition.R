#WISC Fractional Inhibition Data Plots
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)

#import excel ("Active" column as numeric)
WISC_all_bioassays <- read_excel("WISC_all_bioassays.xlsx", 
                                 col_types = c("text", "text", "text", 
                                               "text", "text", "numeric", "numeric", 
                                               "text", "text"))

data_select <- WISC_all_bioassays %>% 
  filter(Active == 1 | Active == 0) #remove NA values

#summary stats by Body Site and Pathogen Type
ds_sum <- data_select %>% 
  group_by(Site, Type) %>% 
  summarise(
    mean = mean(Active), 
    sd = sd(Active), 
    n = n(), 
    SE = sd(Active) / sqrt(n()))

#boxplot for bioactivity data by body site and pathogen type
# x: Pathogen Type | y: Fractional Inhibition | dodge: Body Site
ggplot(ds_sum, aes(x=Type, y=mean, fill=Site)) +
          geom_bar(stat="identity", position = "dodge")+
          geom_errorbar(aes(ymin=mean - SE, ymax=mean + SE), position="dodge")
 

#boxplot for nasal isolate bioactivity data by genus ID and pathogen class
# x: isolate genus ID | y: Fractional Inhibition | facet: Pathogen Class
ds_sum_nasal <- data_select %>% 
  filter(Site == "Nasal") %>% 
  group_by(Type, Genus) %>% 
  summarise(
    mean = mean(Active), 
    sd = sd(Active), 
    n = n(), 
    SE = sd(Active) / sqrt(n()))

ggplot(ds_sum_nasal, aes(x=Genus, y=mean, fill=Genus)) +
  geom_bar(stat="identity", position = "dodge")+
  geom_errorbar(aes(ymin=mean - SE, ymax=mean + SE), position="dodge", width = 0.25)+
  facet_wrap(~Type)+
  scale_x_discrete(limits = c("Coryneb", "Kocuria", "Microco", 
                              "Microba", "Moraxel", "Staphyl", "Strepto"),
                   labels = c("Coryneb" = "Corynebacterium",
                              "Microco" = "Micrococcus",
                              "Microba" = "Microbacterium",
                              "Moraxel" = "Moraxella",
                              "Staphyl" = "Staphylococcus",
                              "Strepto" = "Streptococcus"))

#boxplot for oral isolate bioactivity data by genus ID and pathogen class
# x: isolate genus ID | y: Fractional Inhibition | facet: Pathogen Class
ds_sum_oral <- data_select %>%
  filter(Site == "Oral") %>% 
  group_by(Type, Genus) %>% 
  summarise(
    mean = mean(Active), 
    sd = sd(Active), 
    n = n(), 
    SE = sd(Active) / sqrt(n()))

ggplot(ds_sum_oral, aes(x=Genus, y=mean, fill=Genus)) +
  geom_bar(stat="identity", position = "dodge")+
  geom_errorbar(aes(ymin=mean - SE, ymax=mean + SE), position="dodge", width = 0.25)+
  facet_wrap(~Type)+
  scale_x_discrete(limits = c("Bacillu", "Enteroc", "Granuli", "Lysinib",
                              "Microco", "Rothia", "Strepto"),
                   labels = c("Bacillu" = "Bacillus",
                              "Enteroc" = "Enterococcus",
                              "Granuli" = "Granulicatella",
                              "Lysinib" = "Lysinibacillus",
                              "Microco" = "Micrococcus",
                              "Rothia" = "Rothia",
                              "Strepto" = "Streptococcus"))

