setwd("/home/samantha-capel/Documents/Updated Fst and subsequent results/")

library(SNPRelate)
library(gdsfmt)

snpgdsVCF2GDS(vcf.fn = "gatk.snp.qual_hard_filtered_autosomes_thin.vcf.gz",
              out.fn = "snp_thin.gds")
gds <- openfn.gds("snp_thin.gds", allow.duplicate = T)
snpgdsSummary("snp_thin.gds")
genofile <- snpgdsOpen("snp_thin.gds", allow.duplicate = T)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snpset.id <- c(snpgdsSNPList(genofile)[,1])
snp.id <- sample(snpset.id, round(length(snpset.id)/10))
set.seed(100)
ibd <- snpgdsIBDMLE(genofile, sample.id=sample.id, snp.id = snp.id, 
                    autosome.only = F, maf = NaN, missing.rate = NaN, 
                    num.thread=2)
snpgdsClose(genofile)
sample_ids <- ibd$sample.id
ibd.coeff <- snpgdsIBDSelection(ibd)
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=2)
kin_mat <- matrix(NA, nrow = length(sample_ids), ncol = length(sample_ids),
                  dimnames = list(sample_ids, sample_ids))
for (i in 1:nrow(ibd.coeff)) {
  id1 <- ibd.coeff$ID1[i]
  id2 <- ibd.coeff$ID2[i]
  kin_mat[id1, id2] <- ibd.coeff$kinship[i]
  kin_mat[id2, id1] <- ibd.coeff$kinship[i]
}

library(pheatmap)

pops <- read.csv("popmap_vcf_order.tsv", sep = "\t", header = F)
pops$V2 <- factor(pops$V2, levels = c("NYPRE", "NYPOST", "PAPRE", "PAPOST"))
pops_ordered <- pops[order(pops$V2), ]
kin_mat_ordered <- kin_mat[pops_ordered$V1, pops_ordered$V1]
pdf("Figure S6 SNPRelate Kinship.pdf", height = 6.7, width = 7.35)
pheatmap(kin_mat_ordered, cluster_rows = F, cluster_cols = F, fontsize = 8)
dev.off()


