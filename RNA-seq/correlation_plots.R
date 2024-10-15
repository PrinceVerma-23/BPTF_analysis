library(ggpubr)
library(ggplot2)

bulk_comm_d9_d12 %>% 
  ggplot(aes(x=d9_LFC,y=d12_LFC)) +
  geom_point(color="#0078e0",alpha = 0.4)+
  geom_smooth(method = "lm",se = TRUE,col="black")+
  stat_cor(method = "pearson")+theme_classic()+
  labs(x="Day 9 Log2FC",y="Day 12 Log2FC",
       title = "Correlation of RNA-seq changes in RNAi stem cells, day 9 vs day 12")
