
# generating three input tables
# raw read numbers
rawtable <- read.csv("RIGHTsummary_out_wd.csv", stringsAsFactors = F)
speptable <- rawtable[-(c(2,3,7))]
speptable$bac_read_sum <- speptable$ecoli_assigned_reads + speptable$pputida_assigned_reads

# frequency table
speptable_fr <- speptable
speptable_fr$pputida_assigned_reads <- speptable_fr$pputida_assigned_reads/speptable_fr$bac_read_sum
speptable_fr$ecoli_assigned_reads <- speptable_fr$ecoli_assigned_reads/speptable_fr$bac_read_sum

# raw reads table normalized with genome size
norm_gc <- speptable
norm_gc$pputida_assigned_reads <- norm_gc$pputida_assigned_reads/6327024
norm_gc$ecoli_assigned_reads <- norm_gc$ecoli_assigned_reads/4707887

# getting barplot for average # of reads per bin among replicates for each experiment
library(data.table)
require(plyr)

temp_table <- speptable[, c(2,3,4,18, 19)]



my_bins <- seq(from=1.70, to=1.76, by=0.005)
bar_col <- c("red", "blue")
my_xlabels <- c("0-1.700", 
                "1.700-1.705", 
                "1.705-1.710",
                "1.710-1.715",
                "1.715-1.720",
                "1.720-1.725",
                "1.725-1.730",
                "1.730-1.735",
                "1.735-1.740",
                "1.740-1.745",
                "1.745-1.750",
                "1.750-1.755",
                "1.755-1.760",
                "1.760-inf")
get_tablebins <- function(temp_table, my_bins){
  temp_table$exper_group <- paste(temp_table$Experiment, temp_table$Replicate, sep="_")
  temp_table$bins <- cut(temp_table$Density, c(-Inf, my_bins, Inf), 
                         include.lowest = T, labels=1:14)
  temp_table <- temp_table[,c(4,1,6)]
  colnames(temp_table) <- c("value", "group", "bins")
  
  aver_bins <- ddply(temp_table,.(group,bins),summarise,value = mean(value))
  return(aver_bins)
}
mer_prep <- function(in_table, bac){
  in_table$new_group <- paste(in_table$group, in_table$bins, sep = "_")
  new_in_table <- in_table[,c(4,3)]
  colnames(new_in_table) <- c("group", bac)
  return(new_in_table)
}


aver_binsE <- get_tablebins(temp_table[, c(1,2,3,5)], my_bins)
aver_binsP <- get_tablebins(temp_table[, c(1,2,3,4)], my_bins)

aver_binEm <- mer_prep(aver_binsE, "ecoli")
aver_binPm <- mer_prep(aver_binsP, "pputida")
aver_bacs <- merge(aver_binEm, aver_binPm, all = F)
aver_bacs$exp_id <- gsub("_.*","", aver_bacs$group)
aver_bacs$bin <- as.numeric(gsub(".*_","", aver_bacs$group))
ex_id_vec <- unique(aver_bacs$exp_id)

for(i in 1:length(ex_id_vec)){
  temp_data <- aver_bacs[aver_bacs$exp_id == ex_id_vec[i],c(5,2,3)]
  
  mat1 <- as.matrix(temp_data)
  mat <- t(mat1[ nrow(mat1):1,])
  colnames(mat) <- mat[1,]
  mat <- mat[-1,order(as.integer(colnames(mat)))]
  colnames(mat) <- my_xlabels[as.integer(colnames(mat))]
  
  
  barplot(mat, col=bar_col, legend = c("E. coli", "P. putida"), xlim=c(0, ncol(mat) + 5),
          main=ex_id_vec[i], las = 2)
  
}




####### FUNCTIONS
bin_sta_plot <- function(temp_table, pdf_file_name, rep_num){
  pdf(pdf_file_name)
  my_bins <- seq(from=1.70, to=1.76, by=0.005)
  temp_table$exper_group <- paste(temp_table$Experiment, temp_table$Replicate, sep="_")
  temp_table$bins <- cut(temp_table$Density, c(-Inf, my_bins, Inf), 
                         include.lowest = T, labels=1:14)
  temp_table <- temp_table[,-c(1:3)]
  colnames(temp_table) <- c("basecov", "group", "bins")
  
  aver_bins <- ddply(temp_table,.(group,bins),summarise,basecov = mean(basecov))
  aver_bins$exp_id <- gsub("_.*","", aver_bins$group)
  ex_id_vec <- unique(aver_bins$exp_id)
  graph_table <- function(aver_bins, ex_id_vec_i){
    
    tab <- table(aver_bins[aver_bins$exp_id == ex_id_vec_i,"bins"])
    rep_bins <- as.numeric(names(tab[tab > (rep_num-1)]))
    exp_table <- aver_bins[aver_bins$exp_id == ex_id_vec_i & aver_bins$bins %in% rep_bins,]
    exp_table_m <- ddply(exp_table,.(bins),summarise,basecov = mean(basecov))
    exp_table_sd <- ddply(exp_table,.(bins),summarise,basecov = sd(basecov))
    exp_table_m$sd <- exp_table_sd$basecov
    return(exp_table_m)
  }
  
  l_col <- c("orange", "red", "pink", "black", "blue", "grey")
  my_xlabels <- c("0-1.700", 
                  "1.700-1.705", 
                  "1.705-1.710",
                  "1.710-1.715",
                  "1.715-1.720",
                  "1.720-1.725",
                  "1.725-1.730",
                  "1.730-1.735",
                  "1.735-1.740",
                  "1.740-1.745",
                  "1.745-1.750",
                  "1.750-1.755",
                  "1.755-1.760",
                  "1.760-inf")
  
  
  t <- graph_table(aver_bins, ex_id_vec[1])
  plot(as.numeric(t$bins), t$basecov, type="l", 
       ylim = c(0, max(temp_table$basecov)),
       xlim = c(1, length(my_xlabels)),
       col=l_col[1],xaxt='n', ylab = "Base coverage, mean",
       xlab = "")
  arrows(as.numeric(t$bins),t$basecov-t$sd,
         as.numeric(t$bins),t$basecov+t$sd, 
         code=3, length=0.02, angle = 90, col=l_col[1])
  for(i in 2:length(ex_id_vec)){
    t2 <- graph_table(aver_bins, ex_id_vec[i])
    lines(as.numeric(t2$bins), t2$basecov, type="l", 
          ylim = c(0, (max(t2$basecov+t2$sd))),
          col=l_col[i])
    arrows(as.numeric(t2$bins),t2$basecov-t2$sd,
           as.numeric(t2$bins),t2$basecov+t2$sd, 
           code=3, length=0.02, angle = 90, col=l_col[i])
  }
  legend("topright", legend=ex_id_vec, col=l_col, lty=1, cex=0.8)
  axis(1, at=1:14, labels = my_xlabels,las=2)
  dev.off()
}

bin_sta_plot(norm_gc[,c(2,3,4,19)], "test_normGC_ecoli.pdf", 3)
bin_sta_plot(norm_gc[,c(2,3,4,18)], "test_normGC_pputida.pdf", 3)
