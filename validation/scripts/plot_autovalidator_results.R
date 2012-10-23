#read in data
data = read.csv("../autovalidate_output/mean.csv", head=TRUE, sep="\t", row.names=NULL)
attach(data)

#load scatterplot3d library
library(scatterplot3d)

#plot false positives
pdf("../plots/falsepositive.pdf")
s3d = scatterplot3d(number.of.reads, assemblage.taxa, false.positive, pch=16, highlight.3d=TRUE, type="h", main="False Positive", xlab = "Number of reads", ylab = "Assemblage / All Taxa", zlab = "False Positive Rate")
fit = lm(false.positive ~ number.of.reads+assemblage.taxa)
s3d$plane3d(fit, lty.box = "solid")
dev.off()

#plot false negatives
pdf("../plots/falsenegative.pdf")
s3d = scatterplot3d(number.of.reads, assemblage.taxa, false.negative, pch=16, highlight.3d=TRUE, type="h", main="False Negative", xlab = "Number of reads", ylab = "Assemblage / All Taxa", zlab = "False Negative Rate")
fit = lm(false.negative ~ number.of.reads+assemblage.taxa)
s3d$plane3d(fit, lty.box = "solid")
dev.off()

#plot minspec false negatives
pdf("../plots/minspecfalsenegative.pdf")
s3d = scatterplot3d(number.of.reads, assemblage.taxa, minspec.false.negative, pch=16, highlight.3d=TRUE, type="h", main="False Negative (minspec)", xlab = "Number of reads", ylab = "Assemblage / All Taxa", zlab = "False Negative Rate (minspec)")
fit = lm(minspec.false.negative ~ number.of.reads+assemblage.taxa)
s3d$plane3d(fit, lty.box = "solid")
dev.off()

#plot false taxa removed
pdf("../plots/falsetaxaremoved.pdf")
s3d = scatterplot3d(number.of.reads, assemblage.taxa, false.taxa.removed, pch=16, highlight.3d=TRUE, type="h", main="False Taxa Removed", xlab = "Number of reads", ylab = "Assemblage / All Taxa", zlab = "False Taxa Removed")
fit = lm(false.taxa.removed ~ number.of.reads+assemblage.taxa)
s3d$plane3d(fit, lty.box = "solid")
dev.off()
