# Author: Hu He (hehu@niglas.ac.cn)
# 1. prepare work-----------------------
# set your own work place
setwd("C:/Users/hehu/Documents/data/2019_global_zoobenthos/analysis/data_refined/ISAR") 

# loading packages
library(dplyr)
library(plyr)
library(tidyr)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(data.table)
library(glmmTMB)
library(ggeffects)

# define a theme
mytheme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 10,color = "black"),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title = element_text(size = 10),
        panel.border=element_rect(size = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        aspect.ratio = 2/2)

# 2. upload data -------
data <- readRDS("diversity_inext2.Rds") 

# 2. upload data -------
data2 <-  mutate(data,
                 Area = log10(lakearea*100), # transform km2 to ha
                 Depth = log10(waterdepth),
                 TN = log10(tn),
                 TP = log10(tp),
                 s_a = log10(s_area),
                 s_n = log10(Sn),
                 s_pie = log10(Spie),
                 s_error = (s_area_UCI - s_area_LCI),
                 sn_error = log10(Sn_UCI - Sn_LCI+1),
                 spie_error = (Spie_UCI - Spie_LCI),
                 density = log10(abundance)) %>%
  filter(elev < 2000,source != "Joseline") %>%  # remove 6 reservoirs
  setDT()

data2$climate <- factor(data2$climate,labels = c("Cool","Warm"))

count(data2$sampler_area)

# 3. relationship between nutrients and density--------------------
# glmm
fit_n_np <- glmmTMB(density ~  TN + TP + (1|source) +
                      (1|lat) + (1|long),
                    family = gaussian(),
                    data2)

summary(fit_n_np)

data_n_tn <- as.data.frame(ggpredict(fit_n_np,terms = c("TN")))

data_n_tp <- as.data.frame(ggpredict(fit_n_np,terms = c("TP")))


# plot

plot_n_tn <-  ggplot(data= data2, 
                  aes(TN,density))+ 
  geom_point(aes(color = climate),
             size = 1.5,alpha = 0.4,
             shape = 16)+
  geom_line(data = data_n_tn, aes(x, predicted),size = 1,
            show.legend = F,
            linetype = "dashed") +
  geom_ribbon(data = data_n_tn,aes(x, predicted,
                                   ymin=conf.low, 
                                   ymax=conf.high,
                                    fill = group),
              show.legend = F,
              alpha=0.2) +  
  mytheme +
  labs(x = expression(Log[10](TN)~(mg~L^{-1})), 
       y =  expression(Log[10](Density)~(ind~m^{-2}))) +
  theme(legend.position = c(0.5,0.1),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  scale_color_manual(values = c("mediumpurple","darkorange1"),
                     labels = c("Cool","Warm")) +
  annotate("text",x = 0.25,y = 5.2,size = 3,
           label = expression(italic(P)== '0.739'~(italic(Z)=="-0.333"))) 

plot_n_tn

plot_n_tp <-  ggplot(data= data2, 
                     aes(TP,density))+ 
  geom_point(aes(color = climate),
             size = 1.5,alpha = 0.4,
             shape = 16)+
  geom_line(data = data_n_tp, aes(x, predicted),size = 1,
            show.legend = F,
            linetype = "dashed") +
  geom_ribbon(data = data_n_tp,aes(x, predicted,
                                   ymin=conf.low, 
                                   ymax=conf.high,
                                   fill = group),
              show.legend = F,
              alpha=0.2) +  
  mytheme +
  labs(x = expression(Log[10](TP)~(mg~L^{-1})), 
       y =  expression(Log[10](Density)~(ind~m^{-2}))) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  scale_color_manual(values = c("mediumpurple","darkorange1"),
                     labels = c("Cool","Warm")) +
  annotate("text",x = -1,y = 5.2,size = 3,
           label = expression(italic(P)== '0.560'~(italic(Z)=="0.582"))) 

plot_n_tp



tiff("np_density.tiff",res = 900,compression = "lzw",
     width = 14,height = 8,units = "cm")
ggarrange(plot_n_tn,plot_n_tp,ncol = 2,nrow = 1,labels = "auto",
          common.legend = T,legend = "top")
dev.off()


## Figure.S2 (check the normality of GLMMs)
res <- data <- readRDS("glmm_residuals.Rds") 

tiff("residuals.tiff",res = 900,compression = "lzw",width = 18,height = 18,
     units = "cm")
ggplot(res, aes(x=res)) + 
  geom_histogram(aes(y =..density..), 
                 color="black", fill="white",
                 size = 0.2) +
  mytheme +
  facet_wrap(~index,labeller = label_parsed,scales = "free",ncol = 3) + 
  geom_density(color = 'red',size = 0.3) +
  labs(x = 'Residuals',y = "Density")
dev.off()




