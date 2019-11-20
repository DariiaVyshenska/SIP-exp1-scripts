setwd("/Users/DVyshenska/python_wd")

plot_bac_separately <- function(in_table, stat_name, microbe, yaxis_name){
  e <- c("0Pputida", "40Pputida")
  r <- c("A", "B", "C")
  legend_e <- vector()
  for (e_i in e){legend_e <- append(legend_e, (paste(e_i, r, sep = "")))}
  legend_e <- paste(microbe, legend_e, sep = " ")
  
  
  pdf_file_name <- paste(microbe, "_only_cov_", stat_name, ".pdf", sep = "")
  bac_col <- paste(microbe, "_basecov", sep = "")
  l_col <- c("orange", "red", "pink", "black", "blue", "grey")
  
  pdf(pdf_file_name)
  c <- 1
  plot(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[1])], 
       in_table[(in_table$Experiment == e[1] & in_table$Replicate == r[1]),bac_col], 
       type="l", col=l_col[c], xlab = "Density", ylab = yaxis_name,
       xlim = c(1.7026, 1.7533),
       ylim = c(0,max(in_table[,bac_col])))
  c <- c+1
  for(i in c(2,3)){
    lines(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
          in_table[(in_table$Experiment == e[1] & in_table$Replicate == r[i]), bac_col], 
          col=l_col[c], xlim = c(1.7026, 1.7533),
          ylim = c(0,max(in_table[, bac_col])))
    c <- c+1
  }
  for(i in c(1,2,3)){
    lines(in_table$Density[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
          in_table[(in_table$Experiment == e[2] & in_table$Replicate == r[i]), bac_col], 
          col=l_col[c], xlim = c(1.7026, 1.7533),
          ylim = c(0,max(in_table[, bac_col])))
    c <- c+1
  }
  legend("topright", legend=legend_e, col=l_col, lty=1, cex=0.8)
  dev.off()
  
}


# mean
rawtable <- read.csv("./summary_files/summary_out_bc_mean.csv")
rawtable[is.na(rawtable)] <- 0


plot_bac_separately(rawtable, "mean", "ecoli", "Base coverage")
plot_bac_separately(rawtable, "mean", "pputida", "Base coverage")


# median
rawtable2 <- read.csv("./summary_files/summary_out_bc_median.csv")
rawtable2[is.na(rawtable2)] <- 0

plot_bac_separately(rawtable2, "median", "ecoli", "Base coverage")
plot_bac_separately(rawtable2, "median", "pputida", "Base coverage")
