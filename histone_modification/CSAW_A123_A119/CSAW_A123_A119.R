require(csaw)
require(edgeR)

#Read bam file
param <- readParam(minq=50, max.frag=400, pe="both")
bam.files = c("A119R1_K4.bam", "A119R2_K4.bam", "A123R1_K4.bam", "A123R2_K4.bam")
data <- windowCounts(bam.files, ext=100, spacing=50, width=150, param=param)
#keep <- aveLogCPM(asDGEList(data)) >= -1
#data <- data[keep,]

##Filtering and Normalization
binned <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)
filter.stat <- filterWindows(data, background=binned, type="global")
# a fold change of 3 is necessary for a window to be considered as containing a binding site
keep <- filter.stat$filter > log2(2)
data <- data[keep,]
sum(keep)
#Normalization
normfacs <- normOffsets(binned)
y <- asDGEList(data, norm.factors=normfacs)

#Test differences
design <- model.matrix(~factor(c('A119', 'A119', 'A123', 'A123')))
colnames(design) <- c("intercept", "Strain")
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit, contrast=c(0, 1))

##merge window and output to file
#tol is minimum distance between two binding events are treated as separate sites
merged <- mergeWindows(rowRanges(data), tol=1000L)
tabcom <- combineTests(merged$id, results$table)
write.table(tabcom,file="CSAW.A119_A123.table",sep="\t")
