# External Data Parsing for the Yeast CC Database

This R project stores the data and code used to parse the data listed below
into a format suitable for upload to the yeast CC database. See the data
subdirectories for more details on each.

## Harbison ChIP

[http://younglab.wi.mit.edu/regulatory_code/GWLD.html](http://younglab.wi.mit.edu/regulatory_code/GWLD.html)

> Harbison CT, Gordon DB, Lee TI, Rinaldi NJ, Macisaac KD, Danford TW, Hannett NM, Tagne JB, Reynolds DB, Yoo J, et al. 2004. Transcriptional regulatory code of a eukaryotic genome. Nature 431: 99–104.doi:10.1038/nature02800

## McIsaac Over Expression

[https://idea.research.calicolabs.com/data](https://idea.research.calicolabs.com/data)

> Hackett SR, Baltz EA, Coram M, Wranik BJ, Kim G, Baker A, Fan M, Hendrickson DG, Berndl M, McIsaac RS. 2020. Learning causal networks using inducible transcription factors and transcriptome-wide time series. Mol Syst Biol (in press).

## Kemmeren TFKO

[https://deleteome.holstegelab.nl/](https://deleteome.holstegelab.nl/)

Note: I joined the M values adjusted by removing the first principle component,
which the authors published a paper on showing that this can be confidently
associated with a 'slow growth' signature in a population of the cells.
Additionally, there are some genes for which there are multiple probes, and 
some probes which are noted in the original data to different loci which have
since been merged in the annotation set. In both cases, I deduplicated the data
before uploading by selecting, for a locus with multiple probes, the maximum M,
M-adjusted (Madj), A and minimum p-value.

> Kemmeren P, Sameith K, van de Pasch LA, Benschop JJ, Lenstra TL, Margaritis T, O'Duibhir E, Apweiler E, van Wageningen S, Ko CW, et al. 2014. Large-scale genetic perturbations reveal regulatory networks and an abundance of gene-specific repressors. Cell 157: 740–752.doi:10.1016/j.cell.2014.02.054

## Hu et al TFKO re-analyzed by Reimand et al

[https://pubmed.ncbi.nlm.nih.gov/20385592/](https://pubmed.ncbi.nlm.nih.gov/20385592/)

> Reimand J, Vaquerizas JM, Todd AE, Vilo J, Luscombe NM. Comprehensive reanalysis of transcription factor knockout expression data in Saccharomyces cerevisiae reveals many new targets. Nucleic Acids Res. 2010 Aug;38(14):4768-77. doi: 10.1093/nar/gkq232. Epub 2010 Apr 12. PMID: 20385592; PMCID: PMC2919724.
