setwd("/Users/DVyshenska/python_wd")

rawtable <- read.csv("RIGHTsummary_out_wd.csv")
speptable <- rawtable[-(c(2,3,7))]
speptable$bac_read_sum <- speptable$ecoli_assigned_reads + speptable$pputida_assigned_reads


speptable_fr <- speptable
speptable_fr$pputida_assigned_reads <- speptable_fr$pputida_assigned_reads/speptable_fr$bac_read_sum
speptable_fr$ecoli_assigned_reads <- speptable_fr$ecoli_assigned_reads/speptable_fr$bac_read_sum



#pdf_file_name <- "pputida_only_reads.pdf"
in_table <- speptable_fr


e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
legend_e <- c("pputida A-pputida0", "pputida B-pputida0", "pputida C-pputida0", "pputida A-pputida40", "pputida B-pputida40", "pputida C-pputida40")

l_col <- c("orange", "red", "pink", "black", "blue", "grey")

#pdf(pdf_file_name)

c <- 1
plot(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[1])], 
     in_table$pputida_assigned_reads[(in_table$Experiment == e[1] & in_table$Replicate == r[1])], 
     type="l", col=l_col[c], xlab = "Density", ylab = "Bacteria assigned reads", 
     xlim = c(1.7026, 1.7533),
     ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
c <- c+1
for(i in c(2,3)){
  lines(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
        in_table$pputida_assigned_reads[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
for(i in c(1,2,3)){
  #i <- 1
  lines(in_table$Density[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        in_table$pputida_assigned_reads[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
legend("topright", legend=legend_e, col=l_col, lty=1, cex=0.8)





e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
legend_e <- c("ecoli A-pputida0", "ecoli B-pputida0", "ecoli C-pputida0", "ecoli A-pputida40", "ecoli B-pputida40", "ecoli C-pputida40")




l_col <- c("orange", "red", "pink", "black", "blue", "grey")
c <- 1
plot(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[1])], 
     in_table$ecoli_assigned_reads[(in_table$Experiment == e[1] & in_table$Replicate == r[1])], 
     type="l", col=l_col[c], xlab = "Density", ylab = "Bacteria assigned reads",
     xlim = c(1.7026, 1.7533),
     ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
c <- c+1
for(i in c(2,3)){
  lines(in_table$Density[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
        in_table$ecoli_assigned_reads[(in_table$Experiment == e[1] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
for(i in c(1,2,3)){
  lines(in_table$Density[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        in_table$ecoli_assigned_reads[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
legend("topright", legend=legend_e, col=l_col, lty=1, cex=0.8)


