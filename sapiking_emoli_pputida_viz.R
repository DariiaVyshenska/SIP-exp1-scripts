# normalizing all the table together was not a wise idea -
# it brought all frequencies of bacteria to the same
# number over all the samples.


setwd("/Users/DVyshenska/python_wd")
#library("ggplot2")
#library(phyloseq)
library(preprocessCore)

# importing and wrangling the raw reads data
rawtable <- read.csv("summary_out_wd.csv")
sptable <- rawtable[-(c(2,3,7,22,21))]

spikin <- colnames(sptable[-c(1:7)])
clable <- as.numeric(gsub(".*label.", "",spikin))
sptable[-c(1:7)] <- sptable[-c(1:7)][order(clable)]
colnames(sptable) <- c(colnames(sptable[c(1:7)]), 
                          colnames(sptable[-c(1:7)][order(clable)]))
raw_reads_table <- cbind(sptable[-c(1:7)], rawtable[,c(21,22)])
rownames(raw_reads_table) <- rawtable$Library_Name
raw_reads_table <- as.data.frame(t(raw_reads_table))

# getting frequencies and quantile normalizing the data
# relativizing the data
rel_reads_table <- raw_reads_table
rel_reads_table <- t(t(raw_reads_table) /colSums(raw_reads_table))

# quantile normalizing
data_mat <- data.matrix(rel_reads_table) 
q_reads_table <- normalize.quantiles(data_mat, copy = TRUE)
rownames(q_reads_table) <- rownames(data_mat)
colnames(q_reads_table) <- colnames(data_mat)
q_reads_table <- as.data.frame(q_reads_table)

# adding metadata to normalized reads
norm_table <- as.data.frame(t(q_reads_table))
meta_table <- rawtable
rownames(meta_table) <- meta_table$Library_Name
meta_table <- meta_table[,c(4:6)]
graph_in_table <- merge(meta_table, norm_table, by=0)

write.csv(graph_in_table, "test.csv")

############


pdf_file_name <- "pputida_only_reads.pdf"
in_table <- graph_in_table


e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
legend_e <- c("pputida A-pputida0", "pputida B-pputida0", "pputida C-pputida0", "pputida A-pputida40", "pputida B-pputida40", "pputida C-pputida40")

l_col <- c("orange", "red", "pink", "black", "blue", "grey")

pdf(pdf_file_name)

c <- 1
plot(in_table$Density[(in_table$Experiment == e[2] & in_table$Replicate == r[1])], 
     in_table$pputida_assigned_reads[(in_table$Experiment == e[2] & in_table$Replicate == r[1])], 
     type="l", col=l_col[c], xlab = "Density", ylab = "Bacteria assigned reads", 
     xlim = c(1.7026, 1.7533),
     ylim = c(0,max(in_table$pputida_assigned_reads)))
c <- c+1
for(i in c(2,3)){
  lines(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
        in_table$pputida_assigned_reads[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
for(i in c(1,2,3)){
  lines(in_table$Density[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        in_table$pputida_assigned_reads[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
legend("topright", legend=legend_e, col=l_col, lty=1, cex=0.8)

dev.off()



# yaxis_lab <- "Relative abundance"
# fileO_name <- "norm_fr_spikin.pdf"
# # plotting spikin graphs
# #all_sp_graphs <- function(sptable_fr, fileO_name, yaxis_lab){
#   spikin <- colnames(graph_in_table[-c(1,2,3,15,14)])
#   spikin_leg <- gsub("_", " & ", gsub('Batch146B_p[0-9]+_|pSET152_|pW5Y.AprR_', 
#                                       "", spikin))
# 
#   
#   e <- c("0Pputida", "40Pputida")
#   r <- c("A", "B", "C")
#   line_col <- c(
#     "#56B4E9", 
#     "#009E73", 
#     "#999999", 
#     "#E69F00", 
#     "#F0E442", 
#     "#0072B2", 
#     "#D55E00", 
#     "#CC79A7",
#     "red",
#     "#000000")
#   
#   pdf(fileO_name)
#   par(mar=c(5.1, 4.1, 4.1, 12), xpd=TRUE, mfrow=c(3,1))
#   for(e_i in c(1,2)){
#     for(r_i in c(1,2,3)){
#       e_i <- 1
#       r_i <- 1
#       l_vec <- graph_in_table$Experiment == e[e_i] & graph_in_table$Replicate == r[r_i]
#       
#       plot(graph_in_table[l_vec, "Density"], graph_in_table[l_vec, spikin[1]], 
#            type = "l", col=line_col[1], ylim = c(0, max(graph_in_table[-c(1,2,3,15,14)])), xlim = c(1.7026, 1.7533),
#            xlab = "Density", ylab = yaxis_lab, 
#            main = paste("Experiment: ", e[e_i], ", Replicate ", r[r_i], sep = ""))
#       for(i_t in 2:length(spikin)){
#         lines(graph_in_table[l_vec, "Density"], graph_in_table[l_vec, spikin[i_t]], 
#               type = "l", col=line_col[i_t])
#       }
#       legend(x=1.7555, y=max(graph_in_table[-c(1:7)]), inset=c(-0.2,0), legend=spikin_leg[order(clable)], 
#              col=line_col, lty=1, cex=0.8)
#     }
#   }
#   dev.off()
# #}
# 
# 
# all_sp_graphs(sptable_fr, "ord_spik_per_rep_fr.pdf", "Relative abundance")
# all_sp_graphs(sptable, "ord_spik_per_rep.pdf", "Reads number")
# 
# 
# 
# ######
# #spikin_leg <- gsub("_", " & ", gsub('Batch146B_p[0-9]+_|pSET152_|pW5Y.AprR_', 
# #                                    "", spikin))