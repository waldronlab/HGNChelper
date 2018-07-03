x <- read.csv("gpls_already_tested.csv", as.is=TRUE)
x$valid.after.hgnchelper.frac <- as.numeric(x$valid.after.hgnchelper.frac)
##x$valid.frac <- as.numeric(x$valid.frac)
x <- x[!is.na(x$valid.frac), ]
x <- x[x$valid.frac > 0, ]

## plot(100*x$valid.frac, 100*(x$valid.after.hgnchelper.frac - x$valid.frac))
## abline(a=0, b=1)

## plot(100*x$valid.frac, (x$valid.after.hgnchelper.frac - x$valid.frac)/(1-x$valid.frac),
##      xlab="% Valid Before HGNChelper",
##      ylab="Fraction of invalid gene symbols fixed")

## hist((x$valid.after.hgnchelper.frac - x$valid.frac)/(1-x$valid.frac),
##      main="", xlab="Fraction of invalid gene symbols corrected by HGNChelper")

df <- data.frame(before=cut(100*x$valid.frac, breaks=seq(0, 100, by=20)),
                 fixed=(x$valid.after.hgnchelper.frac - x$valid.frac)/(1-x$valid.frac))
pdf("GEOanalysis.pdf", width=6.5, height=3)
par(mar=c(4,5,2,0.5))
boxplot(fixed ~ before, data=df,
        main=paste("Correcting the annotations of", nrow(df), "GEO platforms"),
     xlab="% Valid Before HGNChelper",
     ylab="Fraction of invalid \n gene symbols fixed",
        col="grey", boxwex=1.1, varwidth=TRUE)
dev.off()

as.Date(head(x$submission_date), format="%b %d %Y")

