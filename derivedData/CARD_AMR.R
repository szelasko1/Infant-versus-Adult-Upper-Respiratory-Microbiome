library(data.table)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)
library(edgeR)

#assemble dataframe
mypath <- "~/CARD_all/"
list_of_files <- list.files(path = mypath,
                            recursive = TRUE,
                            pattern = "\\.txt$", 
                            full.names = TRUE)
DT <- rbindlist(sapply(list_of_files, fread, simplify = FALSE),
                use.names = TRUE, idcol = "FileName")
DT$name  <- gsub(mypath, "", DT$FileName)
DT_select <- DT %>% select(name, `ARO Term`, `Resistomes & Variants: Observed Pathogen(s)`,
                            `Completely Mapped Reads`, `Average MAPQ (Completely Mapped Reads)`, 
                           `Drug Class`, `Resistance Mechanism`, `AMR Gene Family`, `Reference Length`)
colnames(DT_select)[1:9] <- c("Sample", "ARO", "Pathogen", "Reads", "MAPQ", "Drug", "Mechanism", "AMR_Family", "Gene_Length")

#Filter by MAPQ > 10
DT_MAPQ <- DT_select %>% filter(MAPQ > 10)

#FPKM calculation
  fpkm = function (counts, lengths) {
    (counts * 1E9) / (lengths * sum(counts))
  }
  
  DT_FPKM = DT_MAPQ %>%
    mutate(FPKM = fpkm(Reads, Gene_Length))

write.csv(DT_FPKM, file = "~\\DT_select_FPKM.csv")


#heatmap 
data <- read.csv("DT_select_FPKM_heatmap.csv", header=TRUE, sep = ",", row.names = 1, check.names = FALSE)
metadata_col <-  read.csv("DT_select_FPKM_metadata_col.csv", header=TRUE, sep = ",", row.names = 1)
metadata_row <-  read.csv("DT_select_FPKM_metadata_row.csv", header=TRUE, sep = ",", row.names = 1)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_norm <- t(apply(data, 1, cal_z_score))
data_norm <- na.omit(data_norm)

pheatmap(data.matrix(data_norm), cluster_rows=TRUE, cluster_cols=TRUE)


#differential abundance
data_edgeR <- read.csv("DT_select_FPKM_heatmap_edgeR.csv", 
                 header=TRUE, sep = ",", row.names = 1, check.names = FALSE)
metadata <- metadata_col$Site

data_edgeR <- DGEList(counts=data_edgeR, group=metadata)
apply(data_edgeR$counts, 2, sum) 
keep <- rowSums(cpm(data_edgeR)>100) >= 2
data_edgeR <- data_edgeR[keep,]
data_edgeR$samples$lib.size <- colSums(data_edgeR$counts)
data_edgeR <- calcNormFactors(data_edgeR)

d1 <- estimateCommonDisp(data_edgeR)
estimateTagwiseDisp(d1)
list <- estimateDisp(d1)
plotBCV(list)

      #Nasal by Age group
      nasal <- exactTest(list, pair=c("InfantNasal", "AdultNasal"))
      results_nasal <- topTags(nasal, n=nrow(nasal$table), adjust.method="fdr")$table
      results_nasal$topDE <- "NA"
      results_nasal$topDE[results_nasal$logFC > 1 & results_nasal$FDR < 0.05] <- "UP"
      results_nasal$topDE[results_nasal$logFC < -1 & results_nasal$FDR < 0.05] <- "DOWN"
      
      ggplot(data=results_nasal, aes(x=logFC, y=-log10(FDR), color=topDE)) + 
        geom_point(size = 5, alpha = 0.5) +
        scale_colour_discrete(type = c("UP" = "#ff8000", "NA" = "grey", "DOWN" = "#FAD999"),
                                breaks = c("Up", "NA", "Down"))
      
      #Oral by Age group
      oral <- exactTest(list, pair=c("InfantOral", "AdultOral"))
      results_oral <- topTags(oral, n=nrow(oral$table), adjust.method="fdr")$table
      results_oral$topDE <- "NA"
      results_oral$topDE[results_oral$logFC > 1 & results_oral$FDR < 0.05] <- "UP"
      results_oral$topDE[results_oral$logFC < -1 & results_oral$FDR < 0.05] <- "DOWN"
      
      ggplot(data=results_oral, aes(x=logFC, y=-log10(FDR), color=topDE)) + 
        geom_point(size = 5, alpha = 0.5) +
        scale_colour_discrete(type = c("UP" = "#ff8000", "NA" = "grey", "DOWN" = "#FAD999"),
                              breaks = c("Up", "NA", "Down"))

