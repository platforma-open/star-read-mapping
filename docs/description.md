# Overview

Processes bulk RNA sequencing raw data by aligning sequencing reads to a reference genome using STAR (Spliced Transcripts Alignment to a Reference) and quantifying gene expression with featureCounts from the Subread package. The block generates both raw and normalized count matrices suitable for downstream differential expression analysis and other statistical comparisons.

The block performs comprehensive quality control by calculating Principal Component Analysis (PCA) components and sample distance matrices to visually evaluate sample quality, identify potential batch effects, and detect outliers. These QC metrics help assess data quality before proceeding with differential expression analysis or other downstream analyses. The normalized count matrices account for library size differences between samples, enabling accurate comparisons of gene expression levels across experimental conditions.

The block uses STAR for read alignment and featureCounts for gene counting. When using this block in your research, cite the STAR publication (Dobin et al. 2013) and the Subread/featureCounts publication (Liao et al. 2014) listed below.

The following publications describe the methodologies used:

> Dobin, A., Davis, C. A., Schlesinger, F. et al. (2013). STAR: ultrafast universal RNA-seq aligner. _Bioinformatics_ **29**, 15–21 (2013). [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

> Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. _Bioinformatics_ **30**, 923–930 (2014). [https://doi.org/10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)