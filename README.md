Hello, DNA methylation enthusiasts! I'm Elizabeth, a scientific researcher based in Montreal, Quebec, Canada working in the fields of epigenetics and bioinformatics.

This repository will collect miscellaneous functions, scripts, and workflow snippets for analyzing DNA methylation data that I created or adapted during various research projects.

PLEASE CITE ORIGINAL DEVELOPERS OF MY ADAPTED FUNCTIONS.

Below, you’ll find a list of the repository’s contents along with brief descriptions, which will be updated as new items are added.

- getBSseqIndex_DMLtest_EE.R : adapted functions from the DSS R package (DOI: 10.18129/B9.bioc.DSS)
  - getBSseqIndex_EE : adapted getBSseqIndex function to include more than 2 experimental conditions in the CpG filtering steps of the DMLtest function
  - DMLtest_EE : adapted DMLtest function to include more than 2 experimental conditions in the CpG filtering steps

- plotAlphaDistributionOneChr_EE.R : adapted function from the MethylSeekR R package (DOI: 10.18129/B9.bioc.MethylSeekR)
  - Adapted plotAlphaDistributionOneChr function to calculate and extract posterior mean of alpha values from multiple chromosomes at a time
