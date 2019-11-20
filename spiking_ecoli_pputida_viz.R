setwd("/Users/DVyshenska/python_wd")

rawtable <- read.csv("summary_out_wd.csv")
sptable <- rawtable[-(c(2,3,7,22,21))]

###### Plotting density vs DNA concentrations
e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
l_col <- c("orange", "blue", "purple")
plot_dnaconc <- function(e,r,l_col,in_table){
plot(in_table$Density[(in_table$Experiment == e & in_table$Replicate == r[1])], 
     in_table$DNAconc_NgperUL[(in_table$Experiment == e & in_table$Replicate == r[1])], 
     type="l", col=l_col[1], xlab = "Density", ylab = "DNA conc (ng/ul)", main = e,
     ylim = c(0,5))
for(i in c(2,3)){
  lines(in_table$Density[(in_table$Experiment == e & in_table$Replicate == r[i])], 
        in_table$DNAconc_NgperUL[(in_table$Experiment == e & in_table$Replicate == r[i])], 
        col=l_col[i], ylim=c(0,5))
  }
}


pdf("densityVSdnaconc.pdf")
plot_dnaconc(e[1],r,l_col,in_table = sptable)
legend(x=1.73, y= 5, legend=r, col=l_col, lty=1, cex=0.8)
plot_dnaconc(e[2],r,l_col,in_table = sptable)
legend(x=1.74, y= 5, legend=r, col=l_col, lty=1, cex=0.8)
dev.off()

############
# plotting spikin only
all_sp_graphs <- function(sptable_fr, fileO_name, yaxis_lab){
  spikin <- colnames(sptable_fr[-c(1:7)])
  spikin_leg <- gsub("_", " & ", gsub('Batch146B_p[0-9]+_|pSET152_|pW5Y.AprR_', 
                                      "", spikin))
  
  clable <- as.numeric(gsub(".*label.", "",spikin_leg))
  sptable_fr[-c(1:7)] <- sptable_fr[-c(1:7)][order(clable)]
  colnames(sptable_fr) <- c(colnames(sptable_fr[c(1:7)]), 
                            colnames(sptable_fr[-c(1:7)][order(clable)]))
  spikin_o <- spikin[order(clable)]
  
  e <- c("0Pputida", "40Pputida")
  r <- c("A", "B", "C")
  line_col <- c(
    "#56B4E9", 
    "#009E73", 
    "#999999", 
    "#E69F00", 
    "#F0E442", 
    "#0072B2", 
    "#D55E00", 
    "#CC79A7",
    "red",
    "#000000")
  
  pdf(fileO_name)
  par(mar=c(5.1, 4.1, 4.1, 12), xpd=TRUE, mfrow=c(3,1))
  for(e_i in c(1,2)){
    for(r_i in c(1,2,3)){
      l_vec <- sptable_fr$Experiment == e[e_i] & sptable_fr$Replicate == r[r_i]
      
      plot(sptable_fr[l_vec, "Density"], sptable_fr[l_vec, spikin_o[1]], 
           type = "l", col=line_col[1], ylim = c(0, max(sptable_fr[-c(1:7)])), xlim = c(1.7026, 1.7533),
           xlab = "Density", ylab = yaxis_lab, 
           main = paste("Experiment: ", e[e_i], ", Replicate ", r[r_i], sep = ""))
      for(i_t in 2:length(spikin)){
        lines(sptable_fr[l_vec, "Density"], sptable_fr[l_vec, spikin_o[i_t]], 
              type = "l", col=line_col[i_t])
      }
      legend(x=1.7555, y=max(sptable_fr[-c(1:7)]), inset=c(-0.2,0), legend=spikin_leg[order(clable)], 
             col=line_col, lty=1, cex=0.8)
    }
  }
  dev.off()
}

sptable_fr <- sptable
sptable_fr[-c(1:7)] <- sptable[-c(1:7)]/sptable$total_spkn_reads

all_sp_graphs(sptable_fr, "ord_spik_per_rep_fr.pdf", "Relative abundance")
all_sp_graphs(sptable, "ord_spik_per_rep.pdf", "Reads number")

####### Plotting ecoli vs pputida
par(mar=c(5.1, 4.1, 4.1, 2), xpd=TRUE, mfrow=c(1,1))
rawtable <- read.csv("RIGHTsummary_out_wd.csv")
speptable <- rawtable[-(c(2,3,7))]
speptable$bac_read_sum <- speptable$ecoli_assigned_reads + speptable$pputida_assigned_reads

plot_ecpp <- function(e,r,l_col,in_table){
  c <- 1
  plot(in_table$Density[(in_table$Experiment == e & in_table$Replicate == r[1])], 
       in_table$ecoli_assigned_reads[(in_table$Experiment == e & in_table$Replicate == r[1])], 
       type="l", col=l_col[c], xlab = "Density", ylab = "Bacteria assigned reads", main = e,
       xlim = c(1.7026, 1.7533),
       ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
  lines(in_table$Density[(in_table$Experiment == e & in_table$Replicate == r[1])], 
        in_table$pputida_assigned_reads[(in_table$Experiment == e & in_table$Replicate == r[1])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads))) 
  for(i in c(2,3)){
    c <- c+1
    lines(in_table$Density[(in_table$Experiment == e & in_table$Replicate == r[i])], 
          in_table$ecoli_assigned_reads[(in_table$Experiment == e & in_table$Replicate == r[i])], 
          col=l_col[c], xlim = c(1.7026, 1.7533),
          ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
    c <- c+1
    lines(in_table$Density[(in_table$Experiment == e & in_table$Replicate == r[i])], 
          in_table$pputida_assigned_reads[(in_table$Experiment == e & in_table$Replicate == r[i])], 
          col=l_col[c], xlim = c(1.7026, 1.7533),
          ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  }
}


e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
l_col <- c("orange", "blue", "purple", "green", "red", "grey")
legend_ep <- c("ecoli A", "pputida A", "ecoli B", "pputida B", "ecoli C", "pputida C")


plot_ecpp(e[1],r,l_col,speptable)

plot_ecpp(e[2],r,l_col,speptable)



speptable_fr <- speptable
speptable_fr$pputida_assigned_reads <- speptable_fr$pputida_assigned_reads/speptable_fr$bac_read_sum
speptable_fr$ecoli_assigned_reads <- speptable_fr$ecoli_assigned_reads/speptable_fr$bac_read_sum

pdf("freq_ep.pdf")
plot_ecpp(e[1],r,l_col,speptable_fr)
legend(x=1.745, y= 1, legend=legend_ep, col=l_col, lty=1, cex=0.8)
plot_ecpp(e[2],r,l_col,speptable_fr)
legend(x=1.745, y= 1, legend=legend_ep, col=l_col, lty=1, cex=0.8)
dev.off()

pdf("reads_ep.pdf")
plot_ecpp(e[1],r,l_col,speptable)
legend(x=1.740, y= max(speptable$ecoli_assigned_reads, speptable$pputida_assigned_reads), legend=legend_ep, col=l_col, lty=1, cex=0.8)
plot_ecpp(e[2],r,l_col,speptable)
legend(x=1.740, y=max(speptable$ecoli_assigned_reads, speptable$pputida_assigned_reads), legend=legend_ep, col=l_col, lty=1, cex=0.8)
dev.off()

###############
# ecoli only
speptable
speptable_fr
pdf_file_name <- "ecoli_only_reads.pdf"
in_table <- speptable




e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
legend_e <- c("ecoli A-pputida0", "ecoli B-pputida0", "ecoli C-pputida0", "ecoli A-pputida40", "ecoli B-pputida40", "ecoli C-pputida40")

pdf(pdf_file_name)


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
dev.off()


# pputida only
speptable
speptable_fr

pdf_file_name <- "pputida_only_reads.pdf"
in_table <- speptable


e <- c("0Pputida", "40Pputida")
r <- c("A", "B", "C")
legend_e <- c("pputida A-pputida0", "pputida B-pputida0", "pputida C-pputida0", "pputida A-pputida40", "pputida B-pputida40", "pputida C-pputida40")

l_col <- c("orange", "red", "pink", "black", "blue", "grey")

pdf(pdf_file_name)

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
  lines(in_table$Density[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        in_table$pputida_assigned_reads[(in_table$Experiment == e[2] & in_table$Replicate == r[i])], 
        col=l_col[c], xlim = c(1.7026, 1.7533),
        ylim = c(0,max(in_table$ecoli_assigned_reads, in_table$pputida_assigned_reads)))
  c <- c+1
}
legend("topright", legend=legend_e, col=l_col, lty=1, cex=0.8)

dev.off()

##########
# barplot for spikings
# ATTENTION! before reusing this script check if your table has correctly 
# ordered columns


spikin <- colnames(sptable_fr[-c(1:7)])
spikin_leg <- gsub("_", " & ", gsub('Batch146B_p[0-9]+_|pSET152_|pW5Y.AprR_', 
                                    "", spikin))

clable <- as.numeric(gsub(".*label.", "",spikin_leg))
sptable_fr[-c(1:7)] <- sptable_fr[-c(1:7)][order(clable)]
colnames(sptable_fr) <- c(colnames(sptable_fr[c(1:7)]), 
                          colnames(sptable_fr[-c(1:7)][order(clable)]))
spikin_o <- spikin[order(clable)]


full_data <- sptable_fr
bar_col <- c(
  "#56B4E9", 
  "#009E73", 
  "#999999", 
  "#E69F00", 
  "#F0E442", 
  "#0072B2", 
  "#D55E00", 
  "#CC79A7",
  "red",
  "magenta")


pdf("barplot_spikin_l.pdf")
for(e_i in e){
    for(r_i in r){
    temp_data <- full_data[full_data$Experiment == e_i & full_data$Replicate == r_i,c(4,8:17)]
    mat1 <- as.matrix(temp_data[,-1])
    mat <- t(mat1[ nrow(mat1):1,])
    colnames(mat) <- rev(temp_data$Density)
    
    barplot(as.matrix(mat), col=bar_col, legend = spikin_leg, xlim=c(0, nrow(t_temp_data) + 40),
            main=paste(e_i, " & ", r_i, sep =""), las = 2)
    #legend()
    }
}
dev.off()

mat1 <- as.matrix(temp_data[,-1])
mat <- t(mat1[ nrow(mat1):1,])
colnames(mat) <- rev(temp_data$Density)

e_i <- "0Pputida"
r_i <- "A"



