# Adapted plotAlphaDistributionOneChr function from the MethylSeekR R package
# Minor modifications were made to calculate and extract posterior mean of alpha values from multiple chromosomes at a time
# The rest of the function works the same way as originally designed, see DOI: 10.18129/B9.bioc.MethylSeekR
# Reason for the adaptation:
# The original function only works for one chromosome at a time and only outputs a graph, not the values
# I needed the actual values to conduct other analyses and to make my own figures

# PLEASE CITE ORIGINAL DEVELOPERS - DOI: 10.18129/B9.bioc.MethylSeekR

# Function

plotAlphaDistributionOneChr_EE <- function (m, chr.sel, pdfFilename = NULL, num.cores = 1, nCGbin = 101) {
  results <- numeric()
  for (chr in chr.sel) {
    message("Determining alpha distribution for chromosome: ", 
            chr)
    indx <- as.character(seqnames(m)) == chr
    if (sum(indx) < nCGbin) {
      warning(sprintf("Skipping chromosome %s: less than %d covered CpGs", chr, nCGbin))
      next
    }
    T <- as.numeric(values(m[indx])[, 1])
    M <- as.numeric(values(m[indx])[, 2])
    score <- MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores)
    results <- c(results, score)
    rm(indx, T, M, score)
    gc()
  }
  if (!is.null(pdfFilename)) {
    pdf(pdfFilename, width = 5, height = 5)
  }
  hist(results, probability = TRUE, breaks = 30, xlab = "Posterior mean of alpha", main = "")
  if (!is.null(pdfFilename)) 
    dev.off()
  return(results)
}

# Example

library(MethylSeekR)

# meth.gr : genomic ranges object containing methylated and total counts per CpG (same structure as original function)
# chr.sel : vector of characters corresponding to the chromosomes e.g. c("chr1, "chr2", "chr3"...)
# pmoa_distribution : vector containing posterior mean of alpha values

pmoa_distribution = plotAlphaDistributionOneChr_EE(m=meth.gr, chr.sel=chromosomes, pdfFilename = NULL, num.cores = 1, nCGbin = 101)
