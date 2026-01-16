# Adapted getBSseqIndex and DMLtest functions from the DSS R package
# getBSseqIndex is a helper function inside the DMLtest function
# Minor modifications were made to filter CpGs in the DMLtest function based on more than 2 groups
# The rest of the functions work the same way as originally designed, see DOI: 10.18129/B9.bioc.DSS
# getBSseqIndex_EE_3 and DMLtest_EE_3 for 3 groups, getBSseqIndex_EE_4 and DMLtest_EE_4 for 4 groups
# Reason for the adaptation:
# I needed smoothed methylation levels of all groups to be comparable even if DMLs / DMRs are only identified between 2 groups at a time

# PLEASE CITE ORIGINAL DEVELOPERS - DOI: 10.18129/B9.bioc.DSS

# 3 groups

getBSseqIndex_EE_3 <- function (sName, group1, group2, group3) 
{
  check <- function(group, id) {
    thisGrp = paste("group", id, sep = "")
    if (is.character(group)) {
      if (!all(group %in% sName)) 
        stop("Some sample names not found in", thisGrp)
      group <- match(group, sName)
    }
    if (is.numeric(group)) {
      if (min(group) < 1 | max(group) > length(sName)) 
        stop("Some group indices are wrong in", thisGrp)
    }
    else stop(paste("problems with argument", thisGrp))
    group
  }
  group1 <- check(group1, 1)
  group2 <- check(group2, 2)
  group3 <- check(group3, 2)
  if (length(group1) <= 0) 
    stop("group1 not found.")
  if (length(group2) <= 0) 
    stop("group2 not found.")
  if (length(group3) <= 0) 
    stop("group3 not found.")
  list(group1 = group1, group2 = group2, group3 = group3)
}

DMLtest_EE_3 <- function (BSobj, group1, group2, equal.disp = FALSE, smoothing = FALSE, 
                        smoothing.span = 500, ncores, group3) 
{
  if (missing(ncores)) {
    if (.Platform$OS.type == "windows" | Sys.info()["sysname"] == 
        "Windows") 
      ncores = 1
    else ncores = max(detectCores() - 3, 1)
  }
  if (ncores > detectCores()) 
    stop("ncores exceeds the number of CPU cores on the system.")
  tmp <- getBSseqIndex_EE_3(sampleNames(BSobj), group1, group2, group3)
  BS1 <- BSobj[, tmp$group1]
  BS2 <- BSobj[, tmp$group2]
  BS3 <- BSobj[, tmp$group3]
  n1 <- as.array(getBSseq(BS1, "Cov"))
  n2 <- as.array(getBSseq(BS2, "Cov"))
  n3 <- as.array(getBSseq(BS3, "Cov"))
  allpos <- start(BSobj)
  ix1 <- DSS:::hasCoverage(n1, allpos)
  ix2 <- DSS:::hasCoverage(n2, allpos)
  ix3 <- DSS:::hasCoverage(n3, allpos)
  ix <- ix1 & ix2 & ix3
  BS1 <- BS1[ix]
  BS2 <- BS2[ix]
  nreps1 <- dim(BS1)[2]
  nreps2 <- dim(BS2)[2]
  if ((nreps1 == 1 | nreps2 == 1) & !equal.disp) {
    if (!smoothing) 
      stop("There is no biological replicates in at least one condition. Please set smoothing=TRUE or equal.disp=TRUE and retry.")
  }
  if (!smoothing) {
    dmls <- DSS:::DMLtest.noSmooth(BS1, BS2, equal.disp, ncores)
  }
  else {
    dmls <- DSS:::DMLtest.Smooth(BS1, BS2, equal.disp, smoothing.span, 
                                 ncores)
  }
  class(dmls) = c("DMLtest", class(dmls))
  invisible(dmls)
}

# 4 groups

getBSseqIndex_EE_4 <- function (sName, group1, group2, group3, group4) 
{
  check <- function(group, id) {
    thisGrp = paste("group", id, sep = "")
    if (is.character(group)) {
      if (!all(group %in% sName)) 
        stop("Some sample names not found in", thisGrp)
      group <- match(group, sName)
    }
    if (is.numeric(group)) {
      if (min(group) < 1 | max(group) > length(sName)) 
        stop("Some group indices are wrong in", thisGrp)
    }
    else stop(paste("problems with argument", thisGrp))
    group
  }
  group1 <- check(group1, 1)
  group2 <- check(group2, 2)
  group3 <- check(group3, 2)
  group4 <- check(group4, 2)
  if (length(group1) <= 0) 
    stop("group1 not found.")
  if (length(group2) <= 0) 
    stop("group2 not found.")
  if (length(group3) <= 0) 
    stop("group3 not found.")
  if (length(group4) <= 0) 
    stop("group3 not found.")
  list(group1 = group1, group2 = group2, group3 = group3, group4 = group4)
}

DMLtest_EE_4 <- function (BSobj, group1, group2, equal.disp = FALSE, smoothing = FALSE, 
                         smoothing.span = 500, ncores, group3, group4) 
{
  if (missing(ncores)) {
    if (.Platform$OS.type == "windows" | Sys.info()["sysname"] == 
        "Windows") 
      ncores = 1
    else ncores = max(detectCores() - 3, 1)
  }
  if (ncores > detectCores()) 
    stop("ncores exceeds the number of CPU cores on the system.")
  tmp <- getBSseqIndex_EE_4(sampleNames(BSobj), group1, group2, group3, group4)
  BS1 <- BSobj[, tmp$group1]
  BS2 <- BSobj[, tmp$group2]
  BS3 <- BSobj[, tmp$group3]
  BS4 <- BSobj[, tmp$group4]
  n1 <- as.array(getBSseq(BS1, "Cov"))
  n2 <- as.array(getBSseq(BS2, "Cov"))
  n3 <- as.array(getBSseq(BS3, "Cov"))
  n4 <- as.array(getBSseq(BS4, "Cov"))
  allpos <- start(BSobj)
  ix1 <- DSS:::hasCoverage(n1, allpos)
  ix2 <- DSS:::hasCoverage(n2, allpos)
  ix3 <- DSS:::hasCoverage(n3, allpos)
  ix4 <- DSS:::hasCoverage(n4, allpos)
  ix <- ix1 & ix2 & ix3 & ix4
  BS1 <- BS1[ix]
  BS2 <- BS2[ix]
  nreps1 <- dim(BS1)[2]
  nreps2 <- dim(BS2)[2]
  if ((nreps1 == 1 | nreps2 == 1) & !equal.disp) {
    if (!smoothing) 
      stop("There is no biological replicates in at least one condition. Please set smoothing=TRUE or equal.disp=TRUE and retry.")
  }
  if (!smoothing) {
    dmls <- DSS:::DMLtest.noSmooth(BS1, BS2, equal.disp, ncores)
  }
  else {
    dmls <- DSS:::DMLtest.Smooth(BS1, BS2, equal.disp, smoothing.span, 
                                 ncores)
  }
  class(dmls) = c("DMLtest", class(dmls))
  invisible(dmls)
}
