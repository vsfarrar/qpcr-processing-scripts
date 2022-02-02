### qPCR DATA PROCESSING - FUNCTIONIZED SCRIPT
#Author: Victoria Farrar
#Created: 12/5/2020, Last Updated: 2/2/2022

###################### EDIT HERE ############################
#1.define working directory (folder where data and outputs are stored)
setwd("~/Downloads/")

#2.import raw qPCR data 
  #raw data MUST have the following columns in long format: gene, sample, cq, tissue
  #exclude H2O wells from this dataset as they will throw off code 
raw_data <- read.csv("raw_cq_data.csv", stringsAsFactors = F)

#3.import key with relevant treatment groups for ddCT method 
  #key MUST have the following columns: sample, treatment 
key <-read.csv("sample_treatment_key.csv")

#4. define your control group (level of treatment that matches key above)
control_group <- "ctrl"

#5. define project / sample name 
project_name <- "prl-experience_Hipps"

#6. [OPTIONAL] if you want to normalize to specific reference genes, list names here
#names should match how they are entered in dataframe "hipp" above
custom_refgenes <- c("hprt1", "rpl4")

#7. RUN this whole script! 
######################################################################

#load packages
#use pacman to install and load all packages required
if (!require("pacman")) install.packages("pacman") #install pacman if not already installed
pacman::p_load(tidyverse, dplyr,broom,naniar)

#source functions
source("~/Downloads/qpcr_processing_scripts/qpcr_processing_functions.R")

#clean and join data 
colnames(raw_data) <- tolower(colnames(raw_data))
raw_data2 <- 
  raw_data %>%
  mutate(gene = tolower(gene)) %>% #all lowercase genes 
  replace_with_na(replace = list(gene = "")) %>% #blank to NA
  drop_na(gene) 

#join with project key 
#adds treatment group for normalizing expression 

colnames(key) <- tolower(colnames(key))

dat_joined <- inner_join(raw_data2, key) #joins by "sample"

#clean triplicates
dat_clean <- clean_triplicates(dat_joined)

#average cq across triplicates (returns 1 avg cq value for each sample:gene)
dat_avg <- avg_triplicates(dat_clean)


#get average reference gene and dct
dat_dct <- 
  dat_avg %>%
  group_by(sample) %>%
  mutate(is_ref = ifelse(gene %in% reference_genes, 1, 0),
         ref_gene = mean(mean_cq[gene %in% custom_refgenes], na.rm = T)) %>% #get average reference gene
  mutate(dct = mean_cq - ref_gene) #calculate dct

#get dataframe of controls 
control_genes <- 
  dat_dct %>%
  filter(is_ref == 0) %>% #do not calculate control dct for reference genes
  group_by(gene) %>%
  summarise(control_dct = mean(dct[treatment == control_group], na.rm = T)) # control = average dct for vehicle group

#join with control data by genes
dat_dct_controls <- left_join(dat_dct, abc_controls) #joins by gene

# FINAL PRODUCT: calculate ddct, fold change, and log fold change
dat_qpcr_final <- 
  dat_dct_controls %>%
  group_by(sample, gene) %>%
  mutate(ddct = dct-control_dct, 
         fold_change = 2^-(ddct),
         log_fold = log(fold_change)) %>%
  ungroup()

#####################################################
#optional: save processed data file 
write.csv(dat_qpcr_final, paste0(project_name, "_processed_qpcr_data_", current_date))
