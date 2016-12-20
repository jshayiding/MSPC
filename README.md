## MSPC : Bioconductor Package for Multiple Sample Peak Calling

ChIP-seq technology is nowadays routinely used to identify DNA-protein 
interaction and chromatin modifications, while DNase-seq is one of the
most prominent methods to identify open chromatin regions. The output 
of both these techniques is a list of enriched regions, whose statistical 
significance is usually quantified with a score (p-value). Typically, only 
enriched regions whose significance is above a user-defined threshold 
(e.g. p-value < 10-8) are considered. 

The analysis of ChIP-seq samples outputs a number of enriched regions, each 
indicating a protein-DNA interaction or a specific chromatin modification. 
Enriched regions (commonly known as "peaks") are called when the read 
distribution is significantly different from the background and its 
corresponding significance measure (p-value) is below a user-defined threshold. 
When replicate samples are analysed, overlapping enriched regions are expected. 
This repeated evidence can therefore be used to locally lower the minimum 
significance required to accept a peak. Here, we propose a method for joint 
analysis of weak peaks.

Given a set of peaks from (biological or technical) replicates, the method 
combines the p-values of overlapping enriched regions: users can choose a 
threshold on the combined significance of overlapping peaks and set a 
minimum number of replicates where the overlapping peaks should be present. 
The method allows the "rescue" of weak peaks occuring in more than one replicate 
and outputs a new set of enriched regions for each replicate.

Original method is presented in [@Vahid_Jalili_MSPC_2015]. Vahid Jalili, Matteo Matteucci, 
Marco Masseroli,and Marco J. Morelli : Using combined evidence from replicates to evaluate 
ChIP-seq peaks. Bioinformatics 2015, 31(17):2761-2769.doi:[10.1093/bioinformatics/btv293]
(http://bioinformatics.oxfordjournals.org/content/31/17/2761.full). MSPC method was 
implemented in C# programming language and publicly available at https://mspc.codeplex.com. 
Graphical version of MSPC software can be found and available at project site http://musera.codeplex.com.

MSPC packages implements set of functions to retrieve overlapped enriched regions across 
multiple replicates in parallel (a.k.a, pair-wise), statistical method for overlapped regions.
simultaneous presence of an enriched regions in replicates experiment would justify a local 
decrease of the stringency criterion, leveraging on the principal that repeated evidence is
compensating for weak evidence. This packages jointly analyzes the enriched regions of multiple
replicates, distinguishing between biological and technical replicates, and accepting user defined
parameters. Goal of developing R/Bioconductor package of Multiple Sample Peak Calling, to implement 
our algorithm in R and make sure R community get benefit from our method to solve issue raised by high-throughput genomic data. 
