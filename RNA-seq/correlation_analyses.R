# Correlation analyses

library(dplyr)
library(ggplot2)
library(ggpubr)

# Load RNA-seq day 9 and day 12 files
d9_bulk_all <- read.csv("/scratch/pve232/bptf_analyses/bulk_RNA_X1_unc_bptf/day9_files/D9_unc_bptf_bulk_X1_DEG_NO_FILTER.csv",header = TRUE,sep = ",")

d12_bulk_all <- read.csv("/scratch/pve232/bptf_analyses/bulk_RNA_X1_unc_bptf/day12_files/D12_unc_bptf_bulk_X1_DEG_NO_FILTER.csv",header = TRUE,sep = ",")

# Get significant genes for each time point
d9_sig <- d9_bulk_all %>% 
  filter(padj < 0.05) %>% 
  rename("d9_LFC" = log2FoldChange)

d12_sig <- d12_bulk_all %>% 
  filter(padj < 0.05) %>% 
  rename("d12_LFC" = log2FoldChange)

# get common d9 and d12 gene IDs
d9_d12_corr <- inner_join(d9_sig,d12_sig,by = "gene_ID")

#check if data follows normal distribution
par(mfrow = c(1,2))  

# QQ plot for d9_LFC
qqnorm(d9_d12_corr$d9_LFC, main = "QQ Plot - d9_LFC")
qqline(d9_d12_corr$d9_LFC, col = "red")

# QQ plot for d12_LFC
qqnorm(d9_d12_corr$d12_LFC, main = "QQ Plot - d12_LFC")
qqline(d9_d12_corr$d12_LFC, col = "blue")

# correlation
cor.test(d9_d12_corr$d9_LFC,
         d9_d12_corr$d12_LFC,
         method = "spearman",
         exact = FALSE)

# plot and correlation 
d9_d12_corr %>% 
  ggplot(aes(x=d9_LFC,y=d12_LFC)) +
  geom_point(color="#1A50E5",alpha = 0.4)+
  geom_smooth(method = "lm",se = TRUE,level = 0.95,col="black")+
  stat_cor(method = "spearman")+
  theme_classic()+
  labs(x="Day 9 Log2FC",
       y="Day 12 Log2FC",
       title = "Correlation of day 9 vs day 12 DEGs; significant genes only")


# csaw promoter vs RNA-seq with no significance filter correlation
#csaw all
csaw_all <- read.table("/scratch/pve232/bptf_analyses/atac-seq-analyses/csaw_annotated_NO_FILTER_ALL.csv",header = TRUE)

## get gene id col and log2Fc for each gene
csaw_prmtr_all <- csaw_all %>% 
  filter(annotation == "Promoter (<=1kb)") %>% 
  select("gene_ID" = geneId,
         "csaw_LFC" = logFC,
         FDR,annotation,FDR)

d9_comp <- inner_join(x = csaw_prmtr_all,
                      y = d9_bulk_all,by = "gene_ID")

d12_comp <- inner_join(x = csaw_prmtr_all,
                       y = d12_bulk_all,by = "gene_ID")

# plot and correlation 
# day 9 plot
d9_comp %>% 
  ggplot(aes(x=csaw_LFC,y=log2FoldChange)) +
  geom_point(color="#32CDB4",alpha = 0.4)+
  geom_smooth(method = "lm",se = TRUE,level = 0.95,col="black")+
  stat_cor(method = "spearman")+
  theme_classic()+
  labs(x="csaw Log2FC",
       y="Day 9 Log2FC",
       title = "Correlation of csaw promoters vs day 9 DEGs; no sig filter")

#day 12 plot
d12_comp %>% 
  ggplot(aes(x=csaw_LFC,y=log2FoldChange)) +
  geom_point(color="#CD324B",alpha = 0.4)+
  geom_smooth(method = "lm",se = TRUE,level = 0.95,col="black")+
  stat_cor(method = "spearman")+
  theme_classic()+
  labs(x="csaw Log2FC",
       y="Day 12 Log2FC",
       title = "Correlation of csaw promoters vs day 12 DEGs; no sig filter")




