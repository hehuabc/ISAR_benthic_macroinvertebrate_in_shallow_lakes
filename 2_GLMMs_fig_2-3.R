# Author: Hu He (hehu@niglas.ac.cn)
# 0. loading packages --------------
library(dplyr)
library(plyr)
library(tidyr)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(data.table)
library(glmmTMB)
library(ggeffects)


# 1. prepare work-----------------------
# set your own work place
setwd("C:/Users/hehu/Documents/data/2019_global_zoobenthos/analysis/data_refined/ISAR") 

# 2. upload data -------
data <- readRDS("diversity_inext2.Rds") 

# 3. fit Sn, SPIE and interaction ---------------------------------------------
# log10-transformed data
data2 <-  mutate(data,
                Area = log10(lakearea*100), 
                # transform km2 to ha
                Depth = log10(waterdepth),
                TN = log10(tn),
                TP = log10(tp),
                s = log10(s),
                s_a = log10(s_area),
                s_n = log10(Sn),
                s_pie = log10(Spie),
                s_error = 1/((s_area_UCI - s_area_LCI)/3.92+1), 
                # reversal SE for weights, SE = CI/3.92,
                sn_error = 1/((Sn_UCI - Sn_LCI)/3.92+1), 
                # reversal SE for weights
                spie_error = 1/((Spie_UCI - Spie_LCI)/3.92+1), 
                # reversal SE for weights
                density = log10(abundance)) %>% setDT()

data2$climate <- factor(data2$climate,labels = c("Cool","Warm"))

# fit_sarea
fit_s <- glmmTMB(s_a ~  Area  + (Area|source),
                  weights = s_error,
                  data2)
summary(fit_s)

data_s <- as.data.frame(ggpredict(fit_s,terms = "Area"))

# fit_sn
fit_sn <- glmmTMB(s_n ~  Area  + (Area|source),
                  weights = sn_error,
                  data2)
summary(fit_sn)

data_sn <- as.data.frame(ggpredict(fit_sn,terms = "Area"))

# fit_spie
fit_spie <- glmmTMB(s_pie ~ Area + (Area|source),
                  weights = spie_error,
                  data2)

summary(fit_spie)

data_spie <- as.data.frame(ggpredict(fit_spie,terms = "Area"))


# fit area * climate interaction

# s_int
fit_s_int <- glmmTMB(s_a ~ Area * climate + (Area|source),
                      weights = s_error,
                      data2)
summary(fit_s_int)

data_s_int <- as.data.frame(ggpredict(fit_s_int,terms = c("Area","climate")))


# sn_int
fit_sn_int <- glmmTMB(s_n ~ Area * climate + (Area|source),
                      weights = sn_error,
                    data2)
summary(fit_sn_int)

data_sn_int <- as.data.frame(ggpredict(fit_sn_int,terms = c("Area","climate")))

# Spie-int (Spie is too small, so we multiply 1000 to avoid error)
fit_spie_int <- glmmTMB(s_pie*1000 ~ Area * climate +(Area|source),
                     weights = spie_error, 
                     data2)

summary(fit_spie_int)

data_spie_int <- as.data.frame(ggpredict(fit_spie_int,
                                         terms = c("Area","climate"))) %>%
  mutate(predicted = predicted/1000,std.error = std.error/1000,
         conf.low = conf.low/1000,conf.high = conf.high/1000)


# 4. fit productivity and predation strength-----------
# TP 
## here we used data soure as a random intercept,as the AIC value is lower 
### compared with when used source as random slopes, 
### Moreover, there are convergence problems if we include source as random slopes

fit_tp_int <- glmmTMB(TP ~ Area*climate + (1|source),data2)

summary(fit_tp_int)

data_tp_int <- as.data.frame(ggpredict(fit_tp_int,terms = c("Area","climate")))


# tn
fit_tn_int <- glmmTMB(TN ~ Area * climate  + (1|source), 
                  data2)

summary(fit_tn_int)

data_tn_int <- as.data.frame(ggpredict(fit_tn_int,terms = c("Area","climate")))

# density
# here we used both TP and TN as random intercepts, as zoobenthos density was
# determined by both predation (top-down) and lake productivity (bottom-up)

fit_n_int <- glmmTMB(density ~  Area * climate + (1|TN) + (1|TP) ,
               data2)

summary(fit_n_int)

data_n_int <- as.data.frame(ggpredict(fit_n_int,terms = c("Area","climate")))

# 5. get residuals of each GLMM for normality check----------------
res <- cbind(s = residuals(fit_s),
             sn = residuals(fit_sn),
             spie = residuals(fit_spie),
             int_s = residuals(fit_s_int),
      int_sn = residuals(fit_sn_int),
      int_spie = residuals(fit_spie_int)/1000,
      int_tp = residuals(fit_tp_int),
      int_tn = residuals(fit_tn_int),
      int_n = residuals(fit_n_int)) %>% 
  as.data.frame() %>%
  gather(index,res,1:9)

res$index <- factor(res$index,levels = c("s","sn","spie",
                                         "int_s","int_sn","int_spie",
                                         "int_tp","int_tn","int_n"),
                    labels = c(expression(italic(S)~'~ Area'),
                               expression(italic(S)[n]~'~ Area'),
                               expression(italic(S)[PIE]~'~ Area'),
                               expression(italic(S)~'~ Area * Climate'),
                               expression(italic(S)[n]~'~ Area * Climate'),
                               expression(italic(S)[PIE]~'~ Area * Climate'),
                               expression(TP~'~ Area * Climate'),
                               expression(TN~'~ Area * Climate'),
                               expression(Density~'~Area * Climate')))

# save residuals
saveRDS(res,file = "./glmm_residuals.Rds")



# 6. plot Fig. 2-3  ---------------

data2$climate <- factor(data2$climate,labels = c("Cool","Warm"))

# define a plot theme
mytheme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.title = element_text(size = 12),
        panel.border=element_rect(size = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        aspect.ratio = 2/2)


# plot_s
summary(fit_s)

plot_s <-  ggplot(data= data2,aes(Area,s_a)) +
  geom_line(data = data_s, aes(x, predicted),size = 0.5,linetype = "dashed") +
  geom_ribbon(data = data_s, aes(x, predicted,ymin=conf.low, ymax=conf.high), 
              alpha=0.2) +
  geom_point(aes(shape = climate,color = source),size = 1,alpha = 0.8) + 
  geom_errorbar(aes(ymin = log10(s_area_LCI),ymax = log10(s_area_UCI),
                    color = source),
                alpha = 0.3,size = 0.3,show.legend = F)+
  geom_smooth(aes(color = source),show.legend = F,
              method = "lm",
              se = F,
              size = 0.3)  +
  mytheme  +  guides(color = F) +
  labs(x = "", 
       y =  expression(log[10](italic(S)))) +
  theme(legend.position = c(0.8,0.2),
        legend.box.background = element_blank(),
        legend.background = element_blank()) +
  annotate("text",x = 3, y = 2, size = 4,
           label = expression(beta~' = 0.021'))+
  guides(shape=guide_legend(override.aes = list(size=2)))

plot_s


# plot_sn
summary(fit_sn)

plot_sn <-  ggplot(data= data2,aes(Area,s_n)) +
  geom_line(data = data_sn, aes(x, predicted),size = 0.5,linetype = "dashed") +
  geom_ribbon(data = data_sn, aes(x, predicted,ymin=conf.low, ymax=conf.high), 
              alpha=0.2) +
  geom_point(aes(shape = climate,color = source),size = 1,alpha = 0.8) + 
  geom_errorbar(aes(ymin = log10(Sn_LCI),ymax = log10(Sn_UCI),
                    color = source),
                alpha = 0.3,size = 0.3,show.legend = F)+
  geom_smooth(aes(color = source),show.legend = F,
              method = "lm",
              se = F,
              size = 0.3)  +
  mytheme  +  guides(color = F) +
  labs(x = "", 
       y =  expression(log[10](italic(S)[n]))) +
  theme(legend.position = "none",
        legend.box.background = element_blank(),
        legend.background = element_blank()) +
  annotate("text",x = 3, y = 2, size = 4,
           label = expression(beta~' = 0.048'))
  
plot_sn

## plot Spie
summary(fit_spie)
plot_spie <-  ggplot(data= data2,aes(Area,s_pie)) +
  geom_line(data = data_spie, aes(x, predicted),size = 0.5,linetype = "dashed") +
  geom_ribbon(data = data_spie, aes(x, predicted,ymin=conf.low, ymax=conf.high), 
              alpha=0.2) +
  geom_point(aes(shape = climate,color = source),size = 1,alpha = 0.8) + 
  geom_errorbar(aes(ymin = log10(Spie_LCI),ymax = log10(Spie_UCI),
                    color = source),
                alpha = 0.3,size = 0.3)+
  geom_smooth(aes(color = source),
              method = "lm",
              se = F,
              size = 0.3)  +
  mytheme +
  labs(x = expression(log[10](Area)~(ha)), 
       y =  expression(log[10](italic(S)[PIE]))) +
  theme(legend.position = "none",
        legend.background = element_blank()) +
  annotate("text",x = 3, y = 1.5,size = 4, 
           label = expression(beta~' = 0.035'))


plot_spie



## plot s_interaction
summary(fit_s_int)
plot_s_int <-  ggplot(data = data_s_int, aes(x, predicted))+
  geom_line(aes(color = group),linetype = "dashed",size = 0.8) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,fill = group), 
              alpha = 0.3,show.legend = F) +  
  mytheme +
  labs(x = "", 
       y ="") +
  theme(legend.position = c(0.5,0.7),
        legend.direction = "horizontal",
        legend.background = element_blank()) +
  scale_color_manual(values = c("mediumpurple","darkorange1"),
                     labels = c("Cool","Warm"))+
  scale_fill_manual(values = c("mediumpurple","darkorange1")) +
  scale_y_continuous(limits = c(0,2)) +
  annotate("text",x = 3, y = 2,size = 4, 
           label = expression(beta[Warm]~'= 0.050'~~~beta[Cool]~'= - 0.017')) 
plot_s_int



## plot sn_interaction
summary(fit_sn_int)
plot_sn_int <-  ggplot(data = data_sn_int, aes(x, predicted))+
  geom_line(aes(color = group),linetype = "solid",size = 0.8) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,fill = group), 
              alpha = 0.3,show.legend = F) +  
  mytheme +
  labs(x = "", 
       y ="") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("mediumpurple","darkorange1"),
                     labels = c("Cool","Warm"))+
  scale_fill_manual(values = c("mediumpurple","darkorange1")) +
  scale_y_continuous(limits = c(0,2)) +
  annotate("text",x = 3, y = 2,size = 4, 
           label = expression(beta[Warm]~'= 0.151'~~~beta[Cool]~'= - 0.039'))
plot_sn_int

## plot spie_interaction
summary(fit_spie_int)

plot_spie_int <-  ggplot(data = data_spie_int, aes(x, predicted))+
  geom_line(aes(color = group),size = 0.8) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,fill = group), 
              alpha=0.3) +  
  mytheme +
  labs(x = expression(log[10](Area)~(ha)), 
       y ="") +
  theme(legend.position = "none",
        legend.background = element_blank()) +
  scale_color_manual(values = c("mediumpurple","darkorange1"))+
  scale_fill_manual(values = c("mediumpurple","darkorange1")) +
  scale_y_continuous(limits = c(0,1.5)) +
  annotate("text",x = 3, y = 1.5,size = 4, 
           label = expression(beta[Warm]~'= 0.078'~~~beta[Cool]~'= -0.030'))
plot_spie_int

## Figure 2 arrange
tiff("Figure2.tiff",res = 900,compression = "lzw",
     width = 18,height = 24,units = "cm")
ggarrange(plot_s,
          plot_s_int,
          plot_sn, 
          plot_sn_int,
          plot_spie,
          plot_spie_int,
          ncol = 2,nrow = 3,
          labels = c("a","d","b",
                     "e","c","f"),
          align = "hv")
dev.off()


## plot tn,tp,density
summary(fit_n_int)
plot_n <-  ggplot(data= data2, 
                  aes(Area,density))+ 
  geom_line(data = data_n_int, aes(x, predicted,color = group),size = 1) +
  geom_ribbon(data = data_n_int,aes(x, predicted,ymin=conf.low, ymax=conf.high,
                  fill = group), show.legend = F,
              alpha=0.2) +  
  geom_point(aes(color = climate),
             size = 1,alpha = 0.5, shape = 16) +
  mytheme +
  labs(x = expression(log[10](Area)~(ha)), 
       y =  expression(log[10](Density)~(ind~m^{-2}))) +
  theme(legend.position = c(0.5,0.1),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  scale_color_manual(values = c( "mediumpurple","darkorange1")) +
  scale_fill_manual(values = c("mediumpurple","darkorange1")) +
  scale_y_continuous(limits = c(0,6)) +
  annotate("text",x = 3, y = 5.9,size = 3, 
           label = expression(beta[Warm]~'= 0.123'~~~beta[Cool]~'= -0.263'))
plot_n

summary(fit_tp_int)
plot_area_tp <-  ggplot(data= data2, aes(Area,TP))+ 
  geom_line(data = data_tp_int, aes(x, predicted,color = group),size = 1,
            linetype= "dashed") +
  geom_point(aes(color = climate),
             size = 1,alpha = 0.5, shape = 16) +
  geom_ribbon(data = data_tp_int,aes(x, predicted,ymin=conf.low, ymax=conf.high,
                                fill = group), show.legend = F,
              alpha=0.2) +
  mytheme +
  labs(x = "", 
       y =  expression(log[10](TP)~(mg~L^{-1}))) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "vertical") +
  scale_y_continuous(limits = c(-3,1.5)) +
  scale_color_manual(values = c("mediumpurple","darkorange1")) +
  scale_fill_manual(values = c("mediumpurple","darkorange1")) +
  annotate("text",x = 3, y = 1.4,size = 3, 
           label = expression(beta[Warm]~'= -0.012'~~~beta[Cool]~'= 0.073'))
plot_area_tp


summary(fit_tn_int)
plot_area_tn <-  ggplot(data= data2, aes(Area,TN))+ 
  geom_line(data = data_tn_int, aes(x, predicted,color = group),size = 1,
            linetype= "dashed") +
  geom_point(aes(color = climate),size = 1,alpha = 0.5,shape = 16) +
  geom_ribbon(data = data_tn_int,
              aes(x, predicted,
                  ymin=conf.low, 
                  ymax=conf.high,
                  fill = group), show.legend = F,
              alpha=0.2) +
  mytheme +
  labs(x = "", 
       y =  expression(log[10](TN)~(mg~L^{-1}))) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "vertical") +
  scale_color_manual(values = c("mediumpurple","darkorange1")) +
  scale_fill_manual(values = c("mediumpurple","darkorange1")) +
  scale_y_continuous(limits = c(-1,1.5)) +
  annotate("text",x = 3, y = 1.4,size = 3, 
           label = expression(beta[Warm]~'= 0.034'~~~beta[Cool]~'= 0.044'))
plot_area_tn

# Figure 3 arrange
tiff("Figure3.tiff",res = 900,compression = "lzw",
     width = 9,height = 24,
     units = "cm")
ggarrange(plot_area_tp, 
           plot_area_tn,
          plot_n,
           align = "hv",
           ncol = 1,labels = "auto")

dev.off()
