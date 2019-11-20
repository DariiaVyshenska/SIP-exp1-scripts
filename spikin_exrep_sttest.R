require(MESS)
require(pracma)
require(BSDA)

# FUNCTIONS
# do two sample (non-paired) stat test on replicates
funA <- function(in_table, stat_test, file_name){
  # input to the function must be table with first two columns being "Experiment"
  # and "Replicate"
  
  # in-funtion functions
  dispatch <- function(stat_test, x, y){
    if(stat_test == "t"){
      test_res <- t.test(x, y)
    } else if (stat_test == "mw"){
      test_res <- wilcox.test(x,y, exact = F)
    }
    return <- test_res$p.value
  }
  
  # in-function vars
  ex_s <- unique(in_table$Experiment)
  rep <- unique(in_table$Replicate)
  spins <- colnames(in_table[,-c(1,2)])
  
  # building output table
  out_table <- data.frame(matrix(nrow=(ncol(combn(rep,2)) * length(ex_s)), 
                                 ncol=(length(spins))))
  colnames(out_table) <- spins
  rownames(out_table) <- sort(apply(expand.grid(ex_s, apply(combn(rep,2), 2, 
                            function(x) paste(x, collapse = " vs "))), 1, 
                            paste, collapse=" "))
  # populating table
  for(e in ex_s){
    ex_table <- in_table[in_table$Experiment == e,]
    for(s in spins){
      for(a in 1:(length(rep)-1)){
        vec1 <- ex_table[ex_table$Replicate == rep[a],s]
        for(b in (a+1):length(rep)){
          row_id <- paste(e, rep[a], "vs", rep[b], sep = " ")
          vec2 <- ex_table[ex_table$Replicate == rep[b],s]
          out_table[row_id, s] <- dispatch(stat_test, vec1, vec2)
        }
      }
    }
  }
  write.csv(out_table, file_name)
  return(out_table)
}
get_bin_stats <- function(temp_table, file_name, test, method = "NA"){
  
  
  dispatch <- function(x, y, bins , test, method){
    # test options:
    # cor
    # auc
    # sign
    # paired t
    # paired mw
    
    # method options work only on test == "cor":
    # pearson
    # spearman
    
    if(test == "cor"){
      if(method == "pearson"){
        res <- cor.test(x,y, method = "pearson")
      } else if(method == "spearman"){
        res <- cor.test(x,y,method = "spearman", exact = F)
      }
    } else if(test == "auc"){
      if(bins == "NA"){
        print("You need to input 'bins' to dispatch function!")
        exit()
      }
      # insert here test for auc
      
      #AUC1 = trapz(y,bins)
      #AUC2 = trapz(y,bins)
      # some test for AUC
      
    } else if(test == "sign"){
      res <- SIGN.test(x,y, alternative = "two.sided", conf.level = 0.95)
    } else if (test == "paired t"){
      res <- t.test(x,y, paired = T)
    } else if (test == "paired mw"){
      res <- wilcox.test(x,y, paired = T, exact = F)
    }
    return(res$p.value)
  }
  
  
  my_bins <- seq(from=1.70, to=1.76, by=0.005)
  temp_table$bins <- as.numeric(cut(temp_table$Density, c(-Inf, my_bins, Inf), 
                                    include.lowest = T, labels=1:14))
  
  temp_table$temp_label <- paste(temp_table$Experiment, temp_table$Replicate,
                                 temp_table$bins, sep = "_")
  
  agg_table <- aggregate(. ~ temp_label, temp_table[,-c(1:3,14)], mean)
  
  name_l <- strsplit(agg_table$temp_label, "_")
  
  agg_table$Experiment <- sapply(name_l, function(x) x[1])
  agg_table$Replicate <- sapply(name_l, function(x) x[2])
  agg_table$bin <- as.numeric(sapply(name_l, function(x) x[3]))
  
  in_table <- agg_table[,-1]
  write.csv(in_table, "temp_table.csv")
  
  # here I need to loop through spikin and apply next steps per two replicates
  # in one spikin. Then apply stat test
  
  
  
  # in-function vars
  ex_s <- unique(in_table$Experiment)
  rep <- unique(in_table$Replicate)
  spins <- colnames(in_table[,-c(11:13)])
  
  # building output table
  out_table <- data.frame(matrix(nrow=(ncol(combn(rep,2)) * length(ex_s)), 
                                 ncol=(length(spins))))
  colnames(out_table) <- spins
  rownames(out_table) <- sort(apply(expand.grid(ex_s, apply(combn(rep,2), 2, 
                                                            function(x) paste(x, collapse = " vs "))), 1, 
                                    paste, collapse=" "))
  # populating table
  for(e in ex_s){
    
    ex_table <- in_table[in_table$Experiment == e,]
    for(s in spins){
      
      for(a in 1:(length(rep)-1)){
        
        vec1 <- ex_table[ex_table$Replicate == rep[a],s]
        bin1 <- ex_table[ex_table$Replicate == rep[a],"bin"]
        for(b in (a+1):length(rep)){
          
          row_id <- paste(e, rep[a], "vs", rep[b], sep = " ")
          vec2 <- ex_table[ex_table$Replicate == rep[b],s]
          bin2 <- ex_table[ex_table$Replicate == rep[b],"bin"]
          
          # sort vec-s by bin-s
          
          bin1o <- bin1[order(bin1)]
          vec1o <- vec1[order(bin1)]        
          bin2o <- bin2[order(bin2)]
          vec2o <- vec2[order(bin2)]
          
          # select only bins in both replicates (fins duplicated in c(bin1, bin2))
          
          shared_bins <- c(bin1o, bin2o)[duplicated(c(bin1o, bin2o))]
          
          # get only vec-s values that are in both replicates
          
          x <- vec1o[bin1o %in% shared_bins]
          y <- vec2o[bin2o %in% shared_bins]
          
          # input them to statistical test
          
          out_table[row_id, s] <- dispatch(x, y, shared_bins, test = test, method = method)
        }
      }
    }
  }
  
  
  
  write.csv(out_table, file_name)
  #return(out_table)
  
}

# do two sample paired stat test on binned raw read# of replicates
rawtable <- read.csv("RIGHTsummary_out_wd.csv", stringsAsFactors = F, check.names = F)
temp_table <- rawtable[,-c(1:3, 7:10,21,22)]

get_bin_stats(temp_table, file_name = "bin_cor_pearson.csv", test = "cor", method = "pearson")
get_bin_stats(temp_table, file_name = "bin_cor_spearman.csv", test = "cor", method = "spearman")
get_bin_stats(temp_table, file_name = "bin_sign.csv", test = "sign")
get_bin_stats(temp_table, file_name = "bin_paired_t.csv", test = "paired t")
get_bin_stats(temp_table, file_name = "bin_paired_mw.csv", test = "paired mw")



# do same on frequencies

sptable_fr <- temp_table
sptable_fr[-c(1:3)] <- sptable_fr[-c(1:3)]/rawtable$total_spkn_reads


get_bin_stats(sptable_fr, file_name = "frbin_cor_pearson.csv", test = "cor", method = "pearson")
get_bin_stats(sptable_fr, file_name = "frbin_cor_spearman.csv", test = "cor", method = "spearman")
get_bin_stats(sptable_fr, file_name = "frbin_sign.csv", test = "sign")
get_bin_stats(sptable_fr, file_name = "frbin_paired_t.csv", test = "paired t")
get_bin_stats(sptable_fr, file_name = "frbin_paired_mw.csv", test = "paired mw")



# comparing read numbers between all replicates within experiments: non paired test
# on not binned values
rawtable <- read.csv("RIGHTsummary_out_wd.csv", stringsAsFactors = F, check.names = F)
sp_raw_only <- rawtable[,-c(1:3, 6:10,21,22)]

funA(sp_raw_only, "mw", "reads_mw.csv")
funA(sp_raw_only, "t", "reads_welch.csv")

# comparing frequencies of spikin between all replicates within experiments
sptable_fr <- sp_raw_only
sptable_fr[-c(1,2)] <- sp_raw_only[-c(1,2)]/rawtable$total_spkn_reads

funA(sptable_fr, "mw", "fr_mw.csv")
funA(sptable_fr, "t", "fr_welch.csv")





















