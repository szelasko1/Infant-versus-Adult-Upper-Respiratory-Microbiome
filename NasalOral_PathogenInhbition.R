#WISC Fractional Inhibition Data Plots
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

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

