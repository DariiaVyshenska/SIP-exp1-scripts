library(data.table)
require(plyr)


####### FUNCTIONS
allex_graph <- function(test_table, file_name){
  
  my_bins <- seq(from=1.70, to=1.76, by=0.005)
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
  
  test_table$exper_group <- paste(test_table$Experiment, test_table$Replicate, sep="_")
  test_table$bins <- cut(test_table$Density, c(-Inf, my_bins, Inf), include.lowest = T, labels=1:14)
  # write.csv(test_table, "ecoli_bins.csv")
  
  ex_group_ids <- unique(test_table$exper_group)
  l_col <- c("orange", "red", "pink", "black", "blue", "grey")
  
  get_values <- function(test_table, i){
    test_subtable <- data.table(test_table[test_table$exper_group == ex_group_ids[i],-c(1,2,3)])
    colnames(test_subtable) <- c("basecov", "eper_group", "bins")
    subtable_plot <- test_subtable[, mean(basecov), by=bins]
    return(list(x = as.vector(subtable_plot$bins), y = subtable_plot$V1))
  }
  
  pdf(file_name)
  first_curve <- get_values(test_table, 1)
  plot(first_curve$x, first_curve$y, type="l", ylab = "Average coverage in the bin", xlab = "",
       ylim = c(0,max(test_table[,4])),
       xlim = c(1,14), xaxt='n', col=l_col[1])
  axis(1, at=1:14, labels = my_xlabels,las=2)
  for(i in c(2:6)){
    test_curve <- get_values(test_table, i)
    lines(test_curve$x, test_curve$y, 
          col=l_col[i])
    
  }
  legend("topright", legend=ex_group_ids, col=l_col, lty=1, cex=0.8)
  dev.off(which = dev.cur())
  
}
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
##################
###
# importing data
rawtable_median <- read.csv("./summary_files/summary_ep_medianBasecov.csv", 
                            stringsAsFactors = F)
rawtable_mean <- read.csv("./summary_files/summary_ep_meanBasecov.csv", 
                          stringsAsFactors = F)

###
# plotting mean&median base coverage per bin for all experiments
allex_graph(rawtable_mean[,c(4,5,6, 23)], "ecoli_meanBins_allex.pdf")
allex_graph(rawtable_mean[,c(4,5,6, 24)], "pputida_meanBins_allex.pdf")


allex_graph(rawtable_median[,c(4,5,6, 23)], "ecoli_medianBins_allex.pdf")
allex_graph(rawtable_median[,c(4,5,6, 24)], "pputida_medianBins_allex.pdf")

#################
###
# now, find the mean for each bin across replicates, plot the mean
#ecoli,within bin median, between replicates - mean
bin_sta_plot(temp_table = rawtable_median[,c(4,5,6,23)], 
            pdf_file_name = "ecoli_medianBins_mean.pdf",rep_num = 3)
bin_sta_plot(rawtable_median[,c(4,5,6,24)],"pputida_medianBins_mean.pdf",3)


###
# next, plot the same, but without exp 2 repl 1
rawtable_mean_oo <- rawtable_mean[!(rawtable_mean$Experiment == "40Pputida" & rawtable_mean$Replicate == "A"),]
rawtable_median_oo <- rawtable_median[!(rawtable_median$Experiment == "40Pputida" & rawtable_median$Replicate == "A"),]

#ecoli,within bin median, between replicates - mean
bin_sta_plot(rawtable_median_oo[,c(4,5,6, 23)], "ecoli_medianBins_mean_oo.pdf",2)
bin_sta_plot(rawtable_median_oo[,c(4,5,6, 24)], "pputida_medianBins_mean_oo.pdf",2)

