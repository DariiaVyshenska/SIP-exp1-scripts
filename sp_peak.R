require(plyr)
require(dplyr)

get_maxdensity <- function(in_table, file_name){
  in_table$ex_rep_id <- paste(rawtable$Experiment, rawtable$Replicate, sep = "_")
  exrep_un <- unique(in_table$ex_rep_id)
  pop_df <- data.frame(matrix(NA, nrow = length(exrep_un), ncol = 10))
  colnames(pop_df) <- colnames(in_table[,-c(1,12)])
  rownames(pop_df) <- exrep_un
  
  for(i in 1:length(exrep_un)){
    data_table <- as.matrix(in_table[in_table$ex_rep_id == exrep_un[i],-c(1,12)])
    data_density <- in_table$Density[in_table$ex_rep_id == exrep_un[i]]
    pop_df[i,] <- data_density[apply(data_table, 2, function(x) which(x==max(x)) )]  
  }
  
  write.csv(pop_df, file_name)}

# getting density of the raw read # peak for each spikin
rawtable <- read.csv("RIGHTsummary_out_wd.csv", stringsAsFactors = F)
sp_raw_only <- rawtable[,-c(1:5, 7:10,21,22)]
get_maxdensity(sp_raw_only, "maxDensity_rawreads.csv")

# getting density of the frequency peak for each spikin
sptable_fr <- sp_raw_only
sptable_fr[-c(1,12)] <- sp_raw_only[-c(1,12)]/rawtable$total_spkn_reads
get_maxdensity(sptable_fr, "maxDensity_fr.csv")

