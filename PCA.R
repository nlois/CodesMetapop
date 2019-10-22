## PCA Rockhoppers regional

rm(list=ls())
graphics.off()

library(RColorBrewer)
library(SNPRelate)
library(gdsfmt)
library(ggplot2)
#
setwd("/Users/user/Documents/Rockhopper")
#
snpgdsVCF2GDS(vcf.fn="./Sep2018/NORTH/batch_1_allSNPs.vcf", out.fn="NORTH.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
#4572 variants
snpgdsVCF2GDS(vcf.fn="./Sep2018/SOUTH/batch_1_allSNPs.vcf", out.fn="SOUTH.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
#4918 variants
snpgdsVCF2GDS(vcf.fn="./Sep2018/popGen_dataSET/batch_1_allSNPs.vcf", out.fn="ALL.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
#4975 variants
snpgdsVCF2GDS(vcf.fn="./Sep2018/popGen_dataSET_andDEAD/batch_1_allSNPs.vcf", out.fn="andDEAD.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
#4705 variants
snpgdsVCF2GDS(vcf.fn="./Sep2018/NORTH_andDEAD/batch_1_allSNPs.vcf", out.fn="NORTHandDEAD.gds",method = c("biallelic.only"),compress.annotation="ZIP.max", snpfirstdim=FALSE, verbose=TRUE)
#4140 variants
snpgdsSummary("./NORTH.gds")
snpgdsSummary("./SOUTH.gds")
snpgdsSummary("./ALL.gds")
snpgdsSummary("./andDEAD.gds")
snpgdsSummary("./NORTHandDEAD.gds")

genofile1 <- snpgdsOpen("./NORTH.gds")
genofile2 <- snpgdsOpen("./SOUTH.gds")
genofile3 <- snpgdsOpen("./ALL.gds")
genofile4 <- snpgdsOpen("./andDEAD.gds")
genofile5 <- snpgdsOpen("./NORTHandDEAD.gds")

miss1 <- snpgdsSampMissRate(genofile1, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss2 <- snpgdsSampMissRate(genofile2, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss3 <- snpgdsSampMissRate(genofile3, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss4 <- snpgdsSampMissRate(genofile4, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss5 <- snpgdsSampMissRate(genofile5, sample.id=NULL, snp.id=NULL, with.id=TRUE)

pca1 <- snpgdsPCA(gdsobj = genofile1,autosome.only=FALSE)
pca2 <- snpgdsPCA(gdsobj = genofile2,autosome.only=FALSE)
pca3 <- snpgdsPCA(gdsobj = genofile3,autosome.only=FALSE)
pca4 <- snpgdsPCA(gdsobj = genofile4,autosome.only=FALSE)
pca5 <- snpgdsPCA(gdsobj = genofile5,autosome.only=FALSE)

pc.percent1 <- pca1$varprop*100
head(round(pc.percent1, 2))
#5.33 4.83 4.71 4.64 4.58 4.50

pc.percent2 <- pca2$varprop*100
head(round(pc.percent2, 2))
#2.56 1.77 1.74 1.67 1.64 1.63

pc.percent3 <- pca3$varprop*100
head(round(pc.percent3, 2))
#2.65 1.94 1.41 1.36 1.31 1.29

pc.percent4 <- pca4$varprop*100
head(round(pc.percent4, 2))
#2.92 1.66 1.21 1.13 1.10 1.09

pc.percent5 <- pca5$varprop*100
head(round(pc.percent5, 2))
#2.86 2.47 2.47 2.43 2.38 2.36

tab1 <- data.frame(sample.id = pca1$sample.id,
                  EV1 = pca1$eigenvect[,1],    # the first eigenvector
                  EV2 = pca1$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

tab2 <- data.frame(sample.id = pca2$sample.id,
                   EV1 = pca2$eigenvect[,1],    # the first eigenvector
                   EV2 = pca2$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)

tab3 <- data.frame(sample.id = pca3$sample.id,
                   EV1 = pca3$eigenvect[,1],    # the first eigenvector
                   EV2 = pca3$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)

tab4 <- data.frame(sample.id = pca4$sample.id,
                   EV1 = pca4$eigenvect[,1],    # the first eigenvector
                   EV2 = pca4$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)

tab5 <- data.frame(sample.id = pca5$sample.id,
                   EV1 = pca5$eigenvect[,1],    # the first eigenvector
                   EV2 = pca5$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)


#Order in the vcf
vec1 <- c("Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas RV","Malvinas RV","Malvinas RV","Malvinas RV","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado")
vec2 <- c("Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","Franklin","Franklin")
vec3 <- c("Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas RV","Malvinas RV","Malvinas RV","Malvinas RV","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","Franklin","Franklin")
vec4 <- c("Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas RV","Malvinas RV","Malvinas RV","Malvinas RV","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Franklin","Dead TdF","Dead TdF","Dead TdF","Dead TdF","Dead TdF","Dead TdF","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","Chile","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","San Juan","Franklin","Franklin")
vec5 <- c("Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas GJ","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas SLI","Malvinas RV","Malvinas RV","Malvinas RV","Malvinas RV","Dead TdF","Dead TdF","Dead TdF","Dead TdF","Dead TdF","Dead TdF","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Dead SC","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado","Deseado")

tab1f <- cbind(tab1,vec1)
tab2f <- cbind(tab2,vec2)
tab3f <- cbind(tab3,vec3)
tab4f <- cbind(tab4,vec4)
tab5f <- cbind(tab5,vec5)


colorlist<-brewer.pal(9,"Set1")
pdf("rockhopperPCAs_DeadSNPs.pdf",width=10,height=10)
       
plot(tab4f$EV1, tab4f$EV2, ymin=-0.1, ymax=0.1, ylab="PC2 1.66%", xlab="PC1 2.92%",col=colorlist[as.integer(tab4f$vec4)],pch=16,cex=1,main="PCA 4705 SNPs - all populations and Dead")
legend("topleft", legend=levels(tab4f$vec4),cex=1, pch=16,col=colorlist)
#
dev.off() 
#
#
#
