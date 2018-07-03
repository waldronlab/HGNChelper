raw <- read.csv("gpls_already_tested.csv", as.is=TRUE)
raw$valid.after.hgnchelper.frac <- as.numeric(raw$valid.after.hgnchelper.frac)
##raw$valid.frac <- as.numeric(raw$valid.frac)
raw <- raw[!is.na(raw$valid.frac), ]
raw <- raw[raw$valid.frac > 0.2, ]
raw <- raw[!grepl("2014", raw$submission_date), ]
raw$year <- as.numeric(gsub(".+[ ]+", "", raw$submission_date))


pdf("GPL_beforeafter.pdf", width=5, height=5)
n.years <- length(unique(raw$year))
boxplot(valid.frac ~ year, data=raw, boxwex=0.3, xlab="Year", ylab="Fraction of gene symbols that are valid", names=sub("20", "'", sort(unique(raw$year))))
boxplot(valid.after.hgnchelper.frac~year, data=raw, boxwex=0.3, xaxt='n',
        at=(1:n.years)+0.3,
            add=TRUE, col="grey")
legend("bottomleft", legend=c("Before", "After"), pch=c(0, 15), col=c("black", "grey"), lty=-1, bty='n', cex=1)
dev.off()

##system("open GPL_beforeafter.pdf &")
