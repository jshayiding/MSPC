## MSPC : Bioconductor Package for Multiple Sample Peak Calling

The analysis of ChIP-seq samples outputs a number of enriched regions, each indicating a protein-DNA interaction or a specific chromatin modification. Enriched regions (commonly known as "peaks") are called when the read distribution is significantly different from the background and its corresponding significance measure (p-value) is below a user-defined threshold.

When replicate samples are analysed, overlapping enriched regions are expected. This repeated evidence can therefore be used to locally lower the minimum significance required to accept a peak. Here, we propose a method for joint analysis of weak peaks.

Given a set of peaks from (biological or technical) replicates, the method combines the p-values of overlapping enriched regions: users can choose a threshold on the combined significance of overlapping peaks and set a minimum number of replicates where the overlapping peaks should be present. The method allows the "rescue" of weak peaks occuring in more than one replicate and outputs a new set of enriched regions for each replicate.


MSPC packages implements set of functions to retrieve overlapped enriched regions across multiple replicates in parallel (a.k.a, pair-wise), statistical method for overlapped regions. simultaneous presence of an enriched regions in replicates experiment would justify a local decrease of the stringency criterion, leveraging on the principal that repeated evidence is compensating for weak evidence. This packages jointly analyzes the enriched regions of multiple replicates, distinguishing between biological and technical replicates, and accepting user defined parameters.
