# External Data Parsing for the Yeast CC Database

This R project stores the data and code used to parse the data listed below
into a format suitable for upload to the yeast CC database. See the data
subdirectories for more details on each.

## Harbison ChIP

see `data/harbison`

[http://younglab.wi.mit.edu/regulatory_code/GWLD.html](http://younglab.wi.mit.edu/regulatory_code/GWLD.html)

> Harbison CT, Gordon DB, Lee TI, Rinaldi NJ, Macisaac KD, Danford TW, Hannett NM, Tagne JB, Reynolds DB, Yoo J, et al. 2004. Transcriptional regulatory code of a eukaryotic genome. Nature 431: 99–104.doi:10.1038/nature02800

I pulled the 'Excel' versions of both the P values 
(the downloaded file is called files_for_paper_abbr.zip)
and Binding ratios (Ratio_forpaper_abbr.zip). Upon
extracting the files from the zip, I retained the pvalbygene
and ratiobygene files.

## McIsaac Over Expression

[https://idea.research.calicolabs.com/data](https://idea.research.calicolabs.com/data)

> Hackett SR, Baltz EA, Coram M, Wranik BJ, Kim G, Baker A, Fan M, Hendrickson DG, Berndl M, McIsaac RS. 2020. Learning causal networks using inducible transcription factors and transcriptome-wide time series. Mol Syst Biol (in press).

I downloaded the data called "Raw & processed gene expression data", which
yields a file called `idea_tall_expression_data.tsv`.

## Kemmeren TFKO

see `data/kemmeren`

[https://deleteome.holstegelab.nl/](https://deleteome.holstegelab.nl/)

> Kemmeren P, Sameith K, van de Pasch LA, Benschop JJ, Lenstra TL, Margaritis T, O'Duibhir E, Apweiler E, van Wageningen S, Ko CW, et al. 2014. Large-scale genetic perturbations reveal regulatory networks and an abundance of gene-specific repressors. Cell 157: 740–752.doi:10.1016/j.cell.2014.02.054

I downloaded the files:

- `deleteome_all_mutants_ex_wt_var_controls.txt`: this is all of the data __except__ those that had variable WT controls

- `deleteome_all_svd_transformed.txt`: transformed gene expression changes of all mutants with the similarity to the slow 
growth profile removed. "For details see "Cell cycle population effects in perturbation studies", Molecular Systems Biology 2014"

Note: I joined the M values adjusted by removing the first principle component,
which the authors published a paper on showing that this can be confidently
associated with a 'slow growth' signature in a population of the cells.
Additionally, there are some genes for which there are multiple probes, and 
some probes which are noted in the original data to different loci which have
since been merged in the annotation set. In both cases, I deduplicated the data
before uploading by selecting, for a locus with multiple probes, the maximum M,
M-adjusted (Madj), A and minimum p-value.


## Hu et al TFKO re-analyzed by Reimand et al

See `data/hu`

[https://pubmed.ncbi.nlm.nih.gov/20385592/](https://pubmed.ncbi.nlm.nih.gov/20385592/)

> Reimand J, Vaquerizas JM, Todd AE, Vilo J, Luscombe NM. Comprehensive reanalysis of transcription factor knockout expression data in Saccharomyces cerevisiae reveals many new targets. Nucleic Acids Res. 2010 Aug;38(14):4768-77. doi: 10.1093/nar/gkq232. Epub 2010 Apr 12. PMID: 20385592; PMCID: PMC2919724.

This data was downloaded from the supplement 1B and 1C of the paper cited above

# Promoter sets

There are two legacy promoter sets. The first is one Yiming Kang produced in
the Brent lab. I do not have the code for this. The other is from the Mitra
lab and is the intergenic (not open reading frame) regions from sacCer3. This
code is available upon request
