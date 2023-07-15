# RNAtracker

RNAtracker is a computational workflow for the identification of allele-specific RNA stability (asRS) and allele-specific RNA transcription (asRT) genes. 

To reproduce the results in our paper, please follow the steps below. 

### Preliminary steps
1. Download compressed folder of RDDs (RNA-DNA differences) identified in each sample from https://drive.google.com/file/d/1RwRx_6qN2Vs5hMs4Q-QPDu76gcXxKICL/view?usp=sharing. Extract the folder contents with: `tar -zxvf RNAtracker_input_files.tar.gz`
This will give you a folder called `input_files/` which holds the RDDs that were identified in each cell line/sample, as well as the following files:

    * `ALL_COUNT.dat`: has the total number of reads corresponding to each sample, which is information RNAtracker needs for filtering SNVs and calling ASE events.
    * `GENCODE_CNV_filtered.RData`: predicted CNV regions for each cell line (this file was written using R/4.1.0)
    * `ENCODE_cell_line_names`: list of the 16 ENCODE cell lines that we have data for

3. Make sure you have all the neccessary Python packages. The easiest way is to create a conda environment using the provided RNAtracker.yml file: `conda env create -f RNAtracker.yml`

This will create a conda environment named RNAtracker that you should activate before running the python scripts.

### Filter RDD files for reliable SNVs
The RDD files contain RNA-DNA mismatches identified at all positions that have single-nucleotide variants according to dbSNP153. However, we want to filter these RDDs to only include positions that meet the criteria elucidated in Methods. This step is run seperately for each cell line and for each time point. 

The command to run this step is: `python RNAtracker_filter_snvs.py --rdd ${rdd1} ${rdd2} -c ALL_COUNT.dat --prefix ${input_dir}`

Here, `$rdd1` and `$rdd2` are the RDD files for the two replicates in the same cell line at the same timepoint. Since this command needs to be run for each cell line and timepoint seperately, we provide some wrapper scripts for job submission, although these may need to be modified depending on your HPC. 
1. Generate job array with `source submit_filter_snvs_ja1.sh`
2. Submit jobs with `qsub submit_filter_snvs_ja2.sh`

### Identify ASE events
The previous step will output several files, among which is a pickle object containing the filter SNVs that will be used for ASE testing. 

The command to run this step is: `python RNAtracker_ase.py --in-prefix ${dir_ase} --prefix ${out_dir} -p 0.05 --nsigma 1.95996 --filter-test-set 2 10 --exclude-intronic --adjustP fdr_bh`

As before, we provide the wrapper script `submit_find_ASE.sh` for your convenience. Make sure `filtergenes.py` and `libloess_ar.py` are in the same directory as `RNAtracker_ase.py`. 

### Categorize genes as asRS or asRT 
Finally, we can identify our asRS, asRT, and mixed genes/variants. This R script was run in R/4.1.0: `Rscript define_gene_types.R`

This step will output a dataframe, each row of which corresponds to a specific gene at each time point in a particular cell line. 
  * Gene: HGNC gene symbol
  * ASEornonASE: ASE or nonASE
  * #SNV: number of testable SNVs in gene
  * #SNV-SIG: number of ASE SNVs
  * #SNV-AR: number of SNVs with delta_AR <  1.9599 * snv_AR_sigma 
  * Coverage: average normalized coverage across replicates
  * Expected Sigma: Estimated standard deviation of SNV allelic ratios for gene 
  * Allelic_Ratio (by ASE SNVs): average allelic ratio across all ASE SNVs
  * Allelic Ratio (by testable SNVs): average allelic ratio across all testable SNVs
  * SNVs: info on testable SNVs
  * time: timepoint (0h, 2h, or 6h)
  * dataset: cell line
  * gene_label: asRS, asRT, mixed, or not categorized

If you are interested in getting information corresponding to SNVs, you can separate the `SNVs` column into a separate row per SNV (`;` delimited). Each of these rows can be further separated into the following columns (`,` delimited):
  * chr: SNV chromosome
  * pos: SNV position (hg38 coordinates)
  * strand: SNV strand
  * region: genomic region
  * AR: major allele counts/total counts
  * ref_allelic_ratio: reference allele counts/total counts
  * snv_AR_delta: AR - Allelic Ratio (by testable SNVs)
  * snv_AR_sigma: Estimated variance of normalized read counts 
  * major: average major allele normalized counts (across the two replicates)
  * minor: average minor allele normalized counts (across the two replicates)
  * p_orig: SNV ASE p-value
  * p_adj: SNV ASE p-value after BH-adjustment
  * nraw: raw counts (rep1 ref allele count/rep1 alt allele count/rep2 ref allele count/rep2 alt allele count)
  * nnorm: normalized counts (rep1 ref allele count/rep1 alt allele count/rep2 ref allele count/rep2 alt allele count)

You can uncomment the command `all_ase_full = all_ase_full %>% separate_rows(SNVs, sep = ';') %>% filter(SNVs != "")  %>% separate(SNVs, into = c('coord', 'region', 'AR', 'ref_allelic_ratio', 'snv_AR_delta', 'snv_AR_sigma', 'major', 'minor', 'p_orig', 'p_adj', 'nraw', 'nnorm'), sep = ',', convert = TRUE) %>% separate(coord, into = c('chr', 'pos', 'strand'), sep = ':')` in Rscript define_gene_types.R to get these SNV columns. 