# FUNCTIONS

## CLEAN CQ TRIPLICATES
#removes triplicates that are most distant from the mean 
#removes these triplicates only in those sample:gene combos that have a standard deviation > 1 in the triplicates

#import data frame of cleaned cq data (long format)
#data frame requires a "sample", "gene" and "cq" column 

clean_triplicates <- function(x) { #x = cleaned cq dataframe, long format
  cq <- x %>%
    group_by(tissue,sample,gene)%>%
    mutate(mean_cq = mean(cq, na.rm = TRUE), #average cq value across triplicates
           sd_cq = sd(cq, na.rm = TRUE), #sd across triplicates
           dist_mean = abs(cq - mean(cq))) %>% #for each sample, distance from mean cq 
    ungroup()
  
  #for samples where only one well amplified (sd = NA), keep that cq.
  no_sd <- cq %>%
    filter(is.na(sd_cq)) %>%
    group_by(tissue,sample, gene) %>%
    mutate(new_cq = cq)
  
  #when sd => 1, find well with largest distance from mean and drop that cq.
  large_sd <- cq %>%
    filter(sd_cq >= 1) %>%
    group_by(tissue,sample, gene) %>%
    mutate(new_cq = ifelse(dist_mean == max(dist_mean), NA, cq)) %>%
    ungroup()
  
  #when sd < 1, keep cq. 
  small_sd <- cq %>%
    filter(sd_cq < 1) %>%
    group_by(tissue,sample, gene) %>%
    mutate(new_cq = cq)
  
  #stack all data frames using row bind
  rbind.data.frame(large_sd, small_sd, no_sd)
}

## AVERAGE CQ TRIPLICATES 
# across cleaned triplicates, get average cq for each sample:gene combo. 
# returns one value (average cq) for each sample:gene

avg_triplicates <- function(x) { #input must be data run through clean_triplicates function
  x %>%
    group_by(treatment, tissue, sample, gene) %>% 
    summarise(mean_cq = mean(new_cq, na.rm = TRUE), #new_cq column created by clean_triplicates function
              sd_cq = sd(new_cq, na.rm=TRUE)) 
}

## REFERENCE GENE LIST
# list of all possible names of reference genes 

reference_genes <- c("gapdh", "actb_v2", "actb", "rpl4", "ppia", "18s", "hprt1")

## CURRENT DATE 

#grab current date using Sys.Date()
current_date <- as.Date(Sys.Date(), format = "%m/%d/%Y")

