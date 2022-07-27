### Clinical data analysis files

##### All files required for clinical data analysis are in this folder.

The files are:

Script used to load in and standardize VCF files. Paths will need to be updated by user
- vcf_load.ipynb

Variant files from all callers.  
These were generated using the vcf_load script to load in vcf files from each sample from each tool.  
Individual variant files are available upon request. Variant files from each tool should be loaded in, combined and saved as 'compare.callers.vcfs.csv' before loading into analysis files.  
- freebayes.variants.csv
- hc.variants.csv
- ivar.variants.csv
- lofreq.variants.csv
- mutect2.variants.csv
- timo.variants.csv
- varscan.variants.csv

Metadata for all samples and samples with high coverage, used in all downstream analyses.  
- metadata.csv
- highcov.samples.csv

Metadata and accession IDs for all samples. This is also Table S2 in the manuscript.  
- TableS2.xlsx
