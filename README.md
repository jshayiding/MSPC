## MSPC : Bioconductor Package for Multiple Sample Peak Calling

MSPC packages implements set of functions to retrieve overlapped enriched regions across 
multiple replicates in parallel (a.k.a, pair-wise), statistical method for overlapped regions.
simultaneous presence of an enriched regions in replicates experiment would justify a local 
decrease of the stringency criterion, leveraging on the principal that repeated evidence is
compensating for weak evidence. This packages jointly analyzes the enriched regions of multiple
replicates, distinguishing between biological and technical replicates, and accepting user defined
parameters. Goal of developing R/Bioconductor package of Multiple Sample Peak Calling, to implement 
our algorithm in R and make sure R community get benefit from our method to solve issue raised by
high-throughput genomic data. 

Original method is presented in [@Vahid_Jalili_MSPC_2015]. Vahid Jalili, Matteo Matteucci, 
Marco Masseroli,and Marco J. Morelli : Using combined evidence from replicates to evaluate 
ChIP-seq peaks. Bioinformatics 2015, 31(17):2761-2769.doi:[10.1093/bioinformatics/btv293]
(http://bioinformatics.oxfordjournals.org/content/31/17/2761.full). MSPC method was 
implemented in C# programming language and publicly available at https://mspc.codeplex.com. 
Graphical version of MSPC software can be found and available at project site http://musera.codeplex.com.

### RUN MSPC Package

```
Depends:
    R(>= 3.6)
Imports:
    GenomicRanges,
    rtracklayer,
    S4Vectors,
    IRanges,
    dplyr,
    magrittr,
    purrr,
    tidyr,
    ggplot2
```

Local files:

```
    MSPC_Package.R
    readPeakFile.R
    denoise_ER.R
    exportER.R
    getPlot.R
    runMSPC.R
    FDR_stats.R
    utilities.R
    zzz.R

```
### MSPC implementation pipeline
Here is final stage of full-stack implementation of MSPC package:

![screenshot1](output/package_pipeline.PNG)

## Result of Enrichment Analysis

![screenshot1](app_images/crime&realEstate.PNG)

![screenshot2](app_images/all_crime_viz.PNG)
