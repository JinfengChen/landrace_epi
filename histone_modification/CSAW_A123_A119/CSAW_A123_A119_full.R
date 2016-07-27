require(csaw)
require(edgeR)
pdf("CSAW_A123_A119_full.pdf")


##Read bam file
param <- readParam(minq=50, max.frag=400, pe="both")
bam.files = c("A119R1_K4.bam", "A119R2_K4.bam", "A123R1_K4.bam", "A123R2_K4.bam")
data <- windowCounts(bam.files, ext=100, spacing=50, width=150, param=param)
#keep <- aveLogCPM(asDGEList(data)) >= -1
#data <- data[keep,]

##Fragment size distribution
par(mfrow=c(2, 2), mar=c(5, 4, 2, 1.5))
for (i in 1:(length(bam.files))) {
    out <- getPESizes(bam.files[i])
    frag.sizes <- out$sizes[out$sizes<=800]
    hist(frag.sizes, breaks=50, xlab="Fragment sizes (bp)", ylab="Frequency", main=bam.files[i])
    abline(v=400, col="red")
}

##Estimating the average fragment length
max.delay <- 500
dedup.on <- reform(param, dedup=TRUE)
x <- correlateReads(bam.files, max.delay, param=dedup.on)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

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

##Visualze normalizaiton
par(mfrow=c(2, 2), mar=c(5, 4, 2, 1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
for (i in 1:(length(bam.files)-1)) {
    cur.x <- adj.counts[,1]
    cur.y <- adj.counts[,1+i]
    smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
        xlab="Mean of count", ylab="Log2FoldChanges", main=paste("1 vs", i+1))
    all.dist <- diff(log2(normfacs[c(i+1, 1)]))
    abline(h=all.dist, col="red")
}

##Test differences
design <- model.matrix(~factor(c('A119', 'A119', 'A123', 'A123')))
colnames(design) <- c("intercept", "Strain")
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit, contrast=c(0, 1))

##The effect of empirical Bayes stabilisation
par(mfrow=c(2,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
   ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
   ylab=("Biological coefficient of variation"))
plotQLDisp(fit)

#raw and filtered
relevant <- rowSums(assay(data)) >= 20 # some filtering; otherwise, it takes too long.
yo <- asDGEList(data[relevant], norm.factors=normfacs)
yo <- estimateDisp(yo, design)
oo <- order(yo$AveLogCPM)
plot(yo$AveLogCPM[oo], sqrt(yo$trended.dispersion[oo]), type="l", lwd=2,
    ylim=c(0, max(sqrt(yo$trended))), xlab=expression("Ave."~Log[2]~"CPM"),
    ylab=("Biological coefficient of variation"), col="black")
lines(y$AveLogCPM[o], sqrt(y$trended[o]), lwd=2, col="grey")
legend("topright", c("raw", "filtered"), col=c("black", "grey"), lwd=2)

##Replicate similarity
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts <- cpm(y, log=TRUE)
for (top in c(100, 500, 1000, 5000)) {
     plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"), labels=c("sample1.1", "sample1.2", "sample2.1", "sample2.2"), top=top)
}


##merge window and output to file
#tol is minimum distance between two binding events are treated as separate sites
merged <- mergeWindows(rowRanges(data), tol=1000L)
tabcom <- combineTests(merged$id, results$table)
#merged results
regions <- data.frame(seqnames=seqnames(merged$region), start=start(merged$region)-1, ends=end(merged$region))
output <- cbind(regions, tabcom)
write.table(output,file="CSAW_A119_A123.merged.table",sep="\t", row.names = FALSE, quote = FALSE)
#windows based results
windows <- data.frame(seqnames=seqnames(rowRanges(data)), start=start(rowRanges(data))-1, ends=end(rowRanges(data)))
win_output <- cbind(windows, results$table)
write.table(win_output,file="CSAW_A119_A123.windows.table",sep="\t", row.names = FALSE, quote = FALSE)

dev.off()
