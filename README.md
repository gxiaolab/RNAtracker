# RNAtracker

RNAtracker is a computational workflow for the identification of allele-specific RNA stability (asRS) and allele-specific RNA transcription (asRT) genes. 

To reproduce the gene categorization results in our paper, please follow the steps below. The first two scripts read in .txt files of initial hyperparameter values from the data/ folder that were estimated using MLE. Please refer to Supplemental Methods for additional details. 

### System Requirements
All steps only require R (scripts have been tested on R version 4.1.0). No non-standard hardware is required. The following R packages are required: 
- data.table (1.14.2)
- dplyr (1.1.2)
- VGAM (1.1-10)
- scales (1.3.0)
- R.utils (2.10.1)
- tidyr (1.1.3)

The package versions are listed above, but not strictly required to be the same. 

After cloning the repository (`git clone https://github.com/gxiaolab/RNAtracker`), simply run the following scripts from the `scripts/` folder. Downloading the repo should take no more than a minute.


### RNAtracker gene categorization 

1. Determine gene ASE status at 0h

```
Rscript 1_preprocess_0h.R $CELL_LINE $INPUT_DIR $RESULTS_DIR 
Rscript 1_preprocess_0h.R GM12878 ../data/ ../results/
```

2. Categorize genes as asRS, asRT, mixed, or non-ASE. We require the estimated posterior probability of the assigned gene state to be ≥ 0.95 or to be above the first tercile of that state's distribution of posterior probabilities across all genes predicted to be that state. 

```
Rscript 2_categorize_gene.R $CELL_LINE $INPUT_DIR $RESULTS_DIR 
Rscript 2_categorize_gene.R GM12878 ../data/ ../results/
```


3. For timepoints in gene categorizations that reflect ASE (e.g. 0h, 2h, and 6h for asRT; 2h and 6h for asRS), we require at least half of the SNVs for each gene to exhibit allelic imbalance in a consistent direction (e.g. the ref allele has either increased/decreased expression across both replicate 1 and replicate 2). This script will output a .txt file of the final gene categorizations and their SNVs. 

```
Rscript 3_post_categorization_filter.R $CELL_LINE $INPUT_DIR $RESULTS_DIR 
Rscript 3_post_categorization_filter.R GM12878 ../data ../results
```

The above steps should take no more than half an hour to run on a desktop computer. 

### Misc. scripts 
See `misc_scripts/`. 

### ActD RNA-Seq

1. identify validated asRS genes based on allelic imbalance
   ```
        Rscript validate_by_BEAPR_ASE_SNV_null.R GM12878 GM12878gDNA_S3_L006
        Rscript validate_by_BEAPR_ASE_SNV_null.R HCT116 HCT116gDNA_S2_L006
        Rscript validate_by_BEAPR_ASE_SNV_null.R MCF-7 MCF7gDNA_S1_L006
   ```

3. Get number/proportion of validated genes/snvs:
        ```
        Rscript plot_prop_validated_events.R
        ```


### Prime editing
1. run perbase
    ```
        source run_perbase_splice_spec.ja1.sh
        qsub run_perbase_splice_spec.ja2.sh
     ```
1. Combine perbase counts
    ```
        Rscript collect_perbase_splice_spec_counts.R
    ```
3. Plot perbase counts
    ```
        Rscript plot_perbase_spliced_spec_sig_variants.R
     ```

Finally, additional scripts (mostly related to plotting figures) can be found in `misc_scripts/plot_scripts`. 

