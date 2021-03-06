d <- read.table("benchmark_no_stm.csv", sep=";", header=T)
d_t2 <- read.table("benchmark_targetrna2_no_stm.csv", sep=";", header=T)

# CopraRNA 1
c1 <- d$c1_rank
c1_cumulated <- c()

# CopraRNA 2 ooi
c2_ooi <- d$c2_ooi
c2_ooi_cumulated <- c()

# CopraRNA2 ooi cons
c2_ooi_cons <- d$c2_ooi_cons
c2_ooi_cons_cumulated <- c()

# CopraRNA2 ooi ooicons
c2_ooi_ooicons <- d$c2_ooi_ooicons
c2_ooi_ooicons_cumulated <- c()

# CopraRNA2 balanced
c2_bal <- d$c2_bal
c2_bal_cumulated <- c()

# CopraRNA2 balanced cons
c2_bal_cons <- d$c2_bal_cons
c2_bal_cons_cumulated <- c()

# CopraRNA 2 evo
c2_evo <- d$c2_evo
c2_evo_cumulated <- c()

# IntaRNA
intarna <- d$intarna
intarna_cumulated <- c()

# interactome
interactome <- d$interactome
interactome_cumulated <- c()

# TargetRNA2
targetrna <- d_t2$targetrna2
targetrna_cumulated <- c()

# ooi_p-filtered
ooi_p_filtered <- d$ooi_p_filtered
ooi_p_filtered_cumulated <- c()

# ooi_p-cons-fltered
ooi_p_cons_fltered <- d$ooi_p_cons_fltered
ooi_p_cons_fltered_cumulated <- c()

# ooi_ooi_cons_p-filtered
ooi_ooi_cons_p_filtered <- d$ooi_ooi_cons_p_filtered
ooi_ooi_cons_p_filtered_cumulated <- c()

# ooi_cons_p-filtered
ooi_cons_p_filtered <- d$ooi_cons_p-filtered
ooi_cons_p_filtered_cumulated <- c()

for (i in 1:200) {
    # CopraRNA 1
    c1_cumulated <- c(c1_cumulated, length(which(c1<=i)))
    # CopraRNA 2 ooi
    c2_ooi_cumulated <- c(c2_ooi_cumulated, length(which(c2_ooi<=i)))
    # CopraRNA2 ooi cons
    c2_ooi_cons_cumulated <- c(c2_ooi_cons_cumulated,length(which(c2_ooi_cons<=i)))
    # CopraRNA2 ooi ooicons
    c2_ooi_ooicons_cumulated <- c(c2_ooi_ooicons_cumulated, length(which(c2_ooi_ooicons<=i)))
    # CopraRNA2 balanced
    c2_bal_cumulated <- c(c2_bal_cumulated, length(which(c2_bal<=i)))
    # CopraRNA2 balanced cons
    c2_bal_cons_cumulated <- c(c2_bal_cons_cumulated, length(which(c2_bal_cons<=i)))
    # CopraRNA 2 evo 
    c2_evo_cumulated <- c(c2_evo_cumulated, length(which(c2_evo<=i)))
    # IntaRNA
    intarna_cumulated <- c(intarna_cumulated, length(which(intarna<=i)))
    # interactome
    interactome_cumulated <- c(interactome_cumulated, length(which(interactome<=i)))
    # TargetRNA2
    targetrna_cumulated <- c(targetrna_cumulated, length(which(targetrna<=i)))
    # ooi_p-filtered
    ooi_p_filtered_cumulated <- c(ooi_p_filtered_cumulated, length(which(ooi_p_filtered<=i)))
    # ooi_p-cons-fltered
    ooi_p_cons_fltered_cumulated <- c(ooi_p_cons_fltered_cumulated, length(which(ooi_p_cons_fltered<=i)))
    # ooi_ooi_cons_p-filtered
    ooi_ooi_cons_p_filtered_cumulated <- c(ooi_ooi_cons_p_filtered_cumulated, length(which(ooi_ooi_cons_p_filtered<=i)))
    # ooi_cons_p-filtered
    ooi_cons_p_filtered_cumulated <- c(ooi_cons_p_filtered_cumulated, length(which(ooi_cons_p_filtered<=i)))
}

pdf("CopraRNA2_benchmark_no_stm.pdf")
    # CopraRNA 1
    plot(c1_cumulated, type="l", main="", ylab="# True Positives", xlab="# Target predictions per sRNA", ylim=c(0,70), lwd=2, col="#a6cee3")
    # IntaRNA
    lines(intarna_cumulated, col="black", lwd=2)
    # CopraRNA 2 ooi
    lines(c2_ooi_cumulated, col="#b2df8a", lwd=2)
    # CopraRNA2 ooi cons
    lines(c2_ooi_cons_cumulated, col="#33a02c", lwd=2)
    # CopraRNA2 ooi ooicons
    lines(c2_ooi_ooicons_cumulated, col="#fb9a99", lwd=2)
    # CopraRNA2 balanced
    lines(c2_bal_cumulated, col="#e31a1c", lwd=2)
    # CopraRNA2 balanced cons
    lines(c2_bal_cons_cumulated, col="#fdbf6f", lwd=2)
    # CopraRNA 2 evo
    lines(c2_evo_cumulated, col="#ff7f00", lwd=2)
    # interactome
    lines(interactome_cumulated, col="#cab2d6", lwd=2)
    # TargetRNA 2
    lines(targetrna_cumulated, col="#6a3d9a", lwd=2)
    # ooi_p-filtered
    lines(ooi_p_filtered_cumulated, col="#b2df8a", lwd=2, lty=2)
    # ooi_p-cons-fltered
    lines(ooi_p_cons_fltered_cumulated, col="#b2df8a", lwd=2, lty=3)
    # ooi_ooi_cons_p-filtered
    lines(ooi_ooi_cons_p_filtered_cumulated, col="#fb9a99", lwd=2, lty=2)
    # ooi_cons_p-filtered
    lines(ooi_cons_p_filtered_cumulated, col="#33a02c", lwd=2, lty=2)
     
    legend("bottomright", c("CopraRNA 1", "IntaRNA", "CopraRNA 2 ooi", "CopraRNA 2 ooi cons", "CopraRNA 2 ooi ooicons", "CopraRNA 2 balanced", "CopraRNA 2 balanced cons", "CopraRNA 2 evo", "Interactome", "TargetRNA 2", "ooi_p-filtered", "ooi_p-cons-filtered","ooi_ooi_cons_p-filtered","ooi_cons_p-filtered"), col=c("#a6cee3","black","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#b2df8a","#b2df8a","#fb9a99","#33a02c"), lty=c(rep(1,10),2,3,2,2),lwd=2, bty="n", cex=0.6)
dev.off()

pdf("CopraRNA2_benchmark_no_stm_trunc.pdf")
    # CopraRNA 1
    plot(c1_cumulated, type="l", main="", ylab="# True Positives", xlab="# Target predictions per sRNA", ylim=c(0,70), lwd=2, col="#a6cee3")
    # IntaRNA
    lines(intarna_cumulated, col="black", lwd=2)
    # CopraRNA 2 ooi
    lines(c2_ooi_cumulated, col="#b2df8a", lwd=2)
    # CopraRNA2 ooi ooicons
    lines(c2_ooi_ooicons_cumulated, col="#fb9a99", lwd=2)
    # interactome
    lines(interactome_cumulated, col="#cab2d6", lwd=2)
    # TargetRNA 2
    lines(targetrna_cumulated, col="#6a3d9a", lwd=2)

    legend("bottomright", c("CopraRNA 1", "IntaRNA", "CopraRNA 2 ooi", "CopraRNA 2 ooi ooicons", "Interactome (Melamed et al., 2016)", "TargetRNA 2"), col=c("#a6cee3","black","#b2df8a","#fb9a99","#cab2d6","#6a3d9a"), lwd=2, bty="n", cex=0.6)
dev.off()




