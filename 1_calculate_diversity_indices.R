#####################################################################################
# Author : Hu He 
# Email: hehu@niglas.ac.cn
#####################################################################################
# Description : This script can be used to calculate the following diversity indices:

# area-based rarefied richness (S)
# individual-based rarefied richness (S_n)
# a measure of evenness emphasis common species (S_PIE)  

# following 'iNEXT' terminology

#####################################################################################
# 1. loading files
# 2. loading data
# 3. calculate biodiversity metrics using iNEXT
# 4. loading env data and combined with diversity metrics
# 5. Save data

##############################################################################
# LOAD R PACKAGES
rm(list = ls())
require(dplyr)
require(tidyr)
library(vegan)
library(tidyr)
library(readr)
library(vegan)
library(rlang)
library(data.table)

# install_package('iNEXT') # 
library(iNEXT)

############################
# Set path and directories #
############################
# set your own workplace
work_dir <- setwd("C:/Users/hehu/Documents/data/2019_global_zoobenthos/analysis/data_refined")

# 1_file loading----
# first set working directory in Menu/Session/Set working directory/to Project
data_path <- paste0(work_dir, "/species_data")
data_path

# Read all data file names
filenames <- list.files(data_path,pattern="*.csv*", full.names = F)
filenames

# Make list of study ids
filename_roots <- gsub(".csv","",filenames) # remove .csv from filenames
study_ids <- unique(filename_roots)
study_ids

# 2_Read in data----
n_files <- length(unique(study_ids))

data_out <- list()

for (i in 1:n_files){
  
  data_file <- paste(data_path,"/", study_ids[i], ".csv", sep ="")
  
  data_out[[i]] <- read.csv(data_file, header = TRUE, 
                            stringsAsFactors = F,
                            fileEncoding = "UTF-8-BOM")
}


# 3_calculate richness----
all_div_out <- list()

for (i in 1:n_files){
  
  gamma_tab <- data_out[[i]]
  
  gamma_tab[is.na(gamma_tab)] <- 0
  
  # combine lakeno and year into one column
  
  gamma_tab <- unite(gamma_tab,"lakeno",lakeno:year,sep = "")
  
  gamma_tab$lakeno <- as.character(gamma_tab$lakeno)
  
  class(gamma_tab) <- ("data.frame")
  
  # estimate reference n for rarefaction and extrapolations
  
  # abundance
  abundance <- rowSums(gamma_tab[,-c(1)])
  
  # observed richness
  s <- specnumber(gamma_tab[,-c(1)])
  
  # set a common N ( median value of abundance in all lakes)
  n_ref <- 2000 
  
  # use iNEXT to calculate error
  # 1. transfer data.frame to list
  gamma_tab3 <- split(as.matrix(gamma_tab[,-1]), 1:length(gamma_tab$lakeno)) 
  
  # rename list names
  names(gamma_tab3) <- gamma_tab$lakeno
  
  #  calculated individual-based rarefied richness 
  out1 <- estimateD(gamma_tab3, q = c(0,2), # Sn(q = 0); Spie(q = 2)
                    datatype = "abundance", 
                    base="size",
                    level = n_ref) %>%
    select(1,4,6:8)
  
  out2 <- dplyr::rename(out1, 'group' = 'Assemblage',
                        "index" = "Order.q",
                        "LCI" = 'qD.LCL',
                        "UCI" = "qD.UCL")
  
  out2$index <- factor(out2$index, labels = c("Sn", "Spie"))
  
  
  # tranfer long dataframe to wide
  out3 <- gather(out2,key = "index2",value = "value",3:5) %>%
    unite(index,2:3,sep = "_") %>% 
    spread(key = 'index',value = "value") %>% 
    cbind(abundance) %>%
    cbind(s)
  out3 <- dplyr::rename(out3,"Sn" = "Sn_qD", "Spie" = "Spie_qD")
  
  
  # calculate area-based rarefied richness 
  # we set a common area of 0.024 m2, then multipled 0.024 to abundance and got a 
  # commom individual, then folled the step of individual-based rarefaction method
  n_ref_fixarea <- (abundance * 0.024) 
  
  ##########################################
  out <- list()
  for (j in 1:length(n_ref_fixarea)) {
    out5 <- estimateD(gamma_tab3[j], q = 0, 
                      datatype = "abundance", 
                      base="size",
                      level = n_ref_fixarea[j])
    out[[j]] <- out5
  }
  out5 <- bind_rows(out)
  
  out6 <- dplyr::rename(out5, 'group' = 'Assemblage',"index" = "Order.q",
                        "s_area" = "qD","s_area_LCI" = 'qD.LCL',
                        "s_area_UCI" = "qD.UCL") %>%
    select(1,2,6:8) %>%
    full_join(out3,by = "group")
  
  all_div_out[[i]] <- bind_rows(out6)
}


# transfer list into a dataframe, seperate lakeno and year
diversity <- bind_rows(all_div_out) %>%
  tidyr::separate(group,c("lakeno","year"),sep = -4) %>%
  distinct()  

diversity$year <- as.numeric(diversity$year)


# 4_loading env data and combined with diversity metrics -------------------------

env_file <- paste0(work_dir, "/env_data")

# Read all data file names
env_filenames <- list.files(env_file, pattern="*.csv*", full.names = F)

# Make list of study ids
env_filename_roots <- gsub(".csv","",env_filenames) 

env_study_ids <- unique(env_filename_roots)
env_study_ids


env_n_files <- length(unique(env_study_ids))

env_data_out <- list()
for (i in 1:n_files){
  data_file <- paste(env_file,"/", env_study_ids[i], ".csv", sep ="")
  env_data_out[[i]] <- read.csv(data_file, header = TRUE,
                                stringsAsFactors = F,
                                fileEncoding = "UTF-8-BOM")
  
}



env_div_out <- list()

for (i in 1:n_files){
  
  env_div_out[[i]]  <- env_data_out[[i]]
}

env_div_out <- bind_rows(env_div_out)
env_div_out = env_div_out %>% as_tibble()


# 5_save data----------------
## we choose shallow lakes with depth < 10 m and area > 1 ha, 
## we remove high elevation lakes which would caused climate difference in a study
## we eliminate data from Brazil because they are reservoirs and samples got in the littoral zones
data_tot <- full_join(env_div_out,diversity,by = c("lakeno","year")) %>%
  filter(waterdepth < 10,
         lakearea > 0.01,
         elev < 2000,
         s_area > 0,
         country != "Brazil") %>% distinct() 


# write.csv(data_tot,"./ISAR/diversity_inext2.csv",row.names = F)

saveRDS(data_tot,file = "./ISAR/diversity_inext2.Rds")



