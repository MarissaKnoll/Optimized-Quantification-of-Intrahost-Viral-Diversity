## Synthetic Data Analysis

#### All files required for synthetic data analysis are available here

The files are:  

Script used for loading in VCF files
- vcf_load-synflu.ipynb

Combined variant calls for all callers for all samples.  
- flu.synthetic.vcfs.csv. 

Script used for marking variant calls as TP, FP or FN.  
- synthetic_golden.ipynb 

Input files for synthetic_golden.ipynb.  
1. 'golden' csv to mark TPs and FNs
2. metadata for matching sequencing file names with viral load and expected AF
- golden_vcf.csv
- flu_metadata.csv

Output file after marking - input of all synthetic data analyses. 
- flu.synthetic.afdata.csv


