# Author: Hu He (hehu@niglas.ac.cn)
# plot site map----
require(maps)
library(maptools)
library(ggplot2)
library(dplyr)
library(kgc)

# 1. set your own workplace------------------------
setwd("C:/Users/hehu/Documents/data/2019_global_zoobenthos/analysis/data_refined/ISAR") 

# 2. upload data -------
data <- readRDS("diversity_inext2.Rds") 
# define data
data$climate <- factor(data$climate,levels = c("cool","warm"),
                       labels = c("Cool","Warm"))

# 3. world map-----------
mapworld <- borders(database = "world", 
                    colour="gray70", 
                    fill = "gray70",
                    size = 0) #world map
mp_world <- ggplot() + mapworld  

# 4. added points to world map------
tiff("map.tiff",
     compression = "lzw",
     res = 900,
     width = 14,
     height = 9,
     units = "cm")

mp_world + geom_point(data = data, aes(x = long, 
                                          y = lat,
                                       color = climate),
                      size = 1.8,alpha = 0.5,shape = 16) +
  labs(x="Longtitude (¡ãE)",y="Latitude (¡ãN)") + 
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        legend.spacing.x = unit(0.03,"cm"),
        legend.position = c(0.1,0.45),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))  +
  scale_color_manual(values =  c("mediumpurple","darkorange1"))
dev.off()



















