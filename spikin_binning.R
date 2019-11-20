library(data.table)
require(plyr)
####
# FUNCTIONS

table_reorder <- function(table){
  spikin <- colnames(table[-c(1:7)])
  spikin_leg <- gsub("_", " & ", gsub('Batch146B_p[0-9]+_|pSET152_|pW5Y.AprR_', 
                                      "", spikin))
  
  clable <- as.numeric(gsub(".*label.", "",spikin_leg))
  table[-c(1:7)] <- table[-c(1:7)][order(clable)]
  colnames(table) <- c(colnames(table[c(1:7)]), 
                       colnames(table[-c(1:7)][order(clable)]))
  spikin <- spikin[order(clable)]
  spikin_leg <- spikin_leg[order(clable)]
  return(list(spikin=spikin, spikin_leg=spikin_leg, spikin_table=table))
}
get_table <- function(table, exp_group,in_list, i, a){
  subtable <- table[table$exper_group == exp_group[i],c("bins", in_list$spikin[a])]
  colnames(subtable) <- c("bins", "spikin")
  single_binT <- ddply(subtable,.(bins),summarise,spikin = dispatch2(spikin,dis_fun))
  return(single_binT)
}
plot_bin_spikin <- function(table, dis_fun, pdf_file_name){
  dispatch2 <- function(vec, dis_fun){
    if(dis_fun=="mean"){
      return(mean(vec))
    } 
    else if(dis_fun=="median"){
      return(median(vec))}
  }
  my_bins <- seq(from=1.70, to=1.76, by=0.005)
  
  table$exper_group <- paste(table$Experiment, table$Replicate, sep="_")
  table <- table[,-c(1:3,5:7)]
  table$bins <- cut(table$Density, c(-Inf, my_bins, Inf), 
                    include.lowest = T)
  exp_group <- unique(table$exper_group)
  
  pdf(pdf_file_name)
  par(mar=c(5.1, 4.1, 4.1, 12), xpd=TRUE, mfrow=c(3,1))
  l_col <- c("orange", "red", "pink", "black", "blue", "grey")
  my_xlabels <- as.character(levels(table$bins))
  
  for(i in 1:length(exp_group)){
    single_binT <- get_table(table,exp_group,in_list, i,1)
    plot(as.numeric(single_binT$bins), single_binT$spikin, type="l", 
         ylim = c(0, max(table[,-c(1,max(ncol(table)-1), ncol(table))])) ,
         xlim = c(1, length(my_xlabels)), main = exp_group[i],
         col=l_col[1],xaxt='n', ylab = paste("Base coverage within bin, ", 
                                             dis_fun, sep = ""), xlab = "")
    
    for(a in 2:length(in_list$spikin)){
      t2 <- get_table(table,exp_group, in_list, i,a)
      lines(as.numeric(t2$bins), t2$spikin, type="l", col=l_col[a])}
    legend(x=length(levels(table$bins))+1, 
           y=max(table[,-c(1,max(ncol(table)-1), ncol(table))]), 
           legend=in_list$spikin_leg, col=l_col, lty=1, cex=0.8)
    axis(1, at=1:14, labels = my_xlabels,las=2)
  }
  dev.off()
}
####
# import raw data
rawtable <- read.csv("./summary_files/summary_out_wd.csv",stringsAsFactors = F)
sptable <- rawtable[-(c(2,3,7,22,21))]

# getting reordered spikin abundancy table
sptable_fr <- sptable
sptable_fr[-c(1:7)] <- sptable[-c(1:7)]/sptable$total_spkn_reads

in_list <- table_reorder(sptable_fr)
######

dis_fun <- "mean"
table <- in_list$spikin_table
pdf_file_name

plot_bin_spikin(in_list$spikin_table, "mean", "spikin_binned_allexp.pdf")
