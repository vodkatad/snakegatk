loweraf <- as.numeric(snakemake@wildcards[["loweraf"]])
higheraf <- as.numeric(snakemake@wildcards[["higheraf"]])
data <- snakemake@input[["fastamuts"]]
outfile <- snakemake@output[['phyper']]

save.image('pippo.Rdata')

id <- ''
connection <- file(data, open='r')
fileended <- FALSE
samples <- c()
df <- NULL
onefasta <- c()
while (!fileended) {
    line <- readLines(connection, n = 1)
    if (identical(line, character(0))) {
        fileended <- TRUE
    } else {
        if (grepl('>', line)) {
            if (length(onefasta) != 0) {
                if (!is.null(df)) {
                    df <- rbind(df, onefasta)
                } else {
                    df <- onefasta
                }
            }
            id <- gsub('>','', line, fixed=TRUE)
            samples <- c(samples, id)
            onefasta <- c()
        } else {
           if (id == '') {
               stop('Malformed fasta!')
           }
           muts <- as.numeric(unlist(strsplit(line, '\t'))[2])
           names(muts) <- unlist(strsplit(line, '\t'))[1]
           onefasta <- c(onefasta, muts)
        }
    }
}
if (length(onefasta) != 0) {
    df <- rbind(df, onefasta)
}
rownames(df) <- samples


get_phy <- function(x) {
    phy_cl <- phyper(x[4]-1, x[2], x[1]-x[2], x[3], lower.tail=TRUE)
    phy_sc <- phyper(x[7]-1, x[5], x[1]-x[5], x[6], lower.tail=TRUE)
    return(c(phy_cl, phy_sc))
}


apply(df, 1, get_phy)
close(connection)

