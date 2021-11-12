# qiSPI

relative quantification of *in vitro* digested polypeptides and proteins

## sample list
The user must provide Mascot search results and quantification results obtained by Mascot Distiller.
All information must be provided in the `sample_list.csv` (see an example below). You can edit the sample list using, e.g., MS Excel. In any case, make sure to save it as file with comma-separated values (`.csv`) and NOT as `.xlsx` notebook!

| protein_name | substrateID | substrateSeq | digestTime | biological_replicate | order | raw_file | search_result | quantitation_result | peptide_quantification|
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| OPN-C | OPN-C | Hu_OPN_C.fasta | 0 | 1 | Ref | W_Soh_211116_130921_Goe_NC1_0h_R1.raw | F031080.csv | OPN-C_table_matches.txt | OPN-C_table_pept_int.csv |
| OPN-C	| OPN-C	| Hu_OPN_C.fasta | 2	| 2	| C12| 	W_Soh_211116_130921_Goe_OC3_2h_R1.raw	| F031092.csv	| OPN-C_table_matches.txt |	OPN-C_table_pept_int.csv|
| aSyn | P37840aSyn	| P37840aSyn.fasta	| 3	| 1	| C2	| W_Soh_160821_170921_Goe_aSyn_A1_3h_R1.raw | F031252.csv	| aSyn_table_matches.txt| 	aSyn_table_pept_int.csv |

A few general remarks:
- For demonstration purpose, I only provided a few example entries of the sample list. However, it is crucial that you **always provide the full kinetics**, *i.e.*, all files that were processed with Mascot Distiller at once.
- You can put any number of proteins / experiments in the sample list. iqSPI also allows processing multiple proteins/peptides at once. The proteins that are processed together are specified in `SOURCE/config.yaml` (see below).

#### protein_name
The name of the protein/polypeptide to which the respective list entry corresponds. Choose a short and comprehensive name and avoid spaces or special characters.

#### substrateID
The ID under which your quantification results should appear in the final kinetics / quantitative DB. Can but does not have to be identical to *protein_name*.

#### substrateSeq
If you are processing proteins, save their sequence in a single-entry `.fasta` file and deploy it in `INPUT/sequences/`. Specify the name of the `.fasta` file in the `substrateSeq` column of the sample list.
Alternatively, you can paste the entire protein/peptide sequence into the corresponding list entry (only recommended for polypeptides).

#### digestTime
Time point after which the digestion was stopped. **Please provide the time in hours!** For instance, if the digestion time is 15 min, enter 0.25. **Do not put any units in this column!**  
Zero-hours time points are treated as control measurements. Please enter *0* as time point for the control measurements.

#### biological_replicate
Number of the biological (NOT technical!) replicate. In the final kinetics, the mean over all technical replicates is calculated, whereas biological replicates are displayed separately. Please have a look at the full sample list (`INPUT/sample_list.csv`) for clarification.

#### order
Order in which the .raw files were processed in Mascot Distiller. **This column is important for matching the time points to the correct intensities!**. The first sample processed should be labelled *Ref*. All subsequent samples should be labelled as *C1, C2, C3, ....*. Please have a look at the full sample list (`INPUT/sample_list.csv`) for clarification.

#### raw_file
Provide the **full** name (including .raw suffix) of the .raw file for the respective sample. Note that you do NOT have to provide the .raw files for the pipeline. The name of the .raw file is solely required as information during sample processing.

#### search_result
Provide the **full** name (including .csv suffix) of the sample's search result file. Deploy all search result files in the `INPUT/search_results/` folder.

#### quantitation_result
Put the file name of the table matches in this entry. Table matches should be provided as `.txt`, `.csv` or `.html` files. Note that the pipeline is tested for the `.txt` file format.  
Deploy the respective file in `INPUT/quantitation_results/`. **They should not contain any headers!**  
If you are processing multiple proteins/polypeptides at once, make sure to rename the files so that each protein/polypeptide is processed using a distinct quantification file. Please have a look at the full sample list (`INPUT/sample_list.csv`) for clarification.

#### peptide_quantification
Put the file name of the peptide quantification file in this entry. Table matches should be provided as `.txt`, `.csv` files. Note that the pipeline is tested for the `.csv` file format.  
Deploy the respective file in `INPUT/quantitation_results/`.  **They should not contain any headers!**  
If you are processing multiple proteins/polypeptides at once, make sure to rename the files so that each protein/polypeptide is processed using a distinct quantification file. Please have a look at the full sample list (`INPUT/sample_list.csv`) for clarification.

## config
Please have a look at `INPUT/config.yaml`.
```
protein_name:
  - aSyn
  - OPN-C
KK: 0.05
IONscore: 20
thrDif: 0.3
thrDifpcp: 0.3
PCPthresh: 50
PSPthresh: 90
rm_tmp: "no"
```

Please enter the names of all proteins that should be processed together under `protein_name` following the syntax that is given.
**The protein names have to be identical to those in the protein_name column in the sample list!**  

Next, the config file contains hyperparameters for sample processing. Leave them unchanged unless you have a reason to change them.
The last parameter, `rm_tmp` specifies whether temporarily generated files (unfiltered kinetics for each sample) should be removed upon pipeline completion. This is recommended since intermediate files are redundant and occupy a lot of disk space. If you want to remove them, enter `"yes"`. Make sure to use `"` since plain `yes` or  `no` are treated as Boolean in `.yaml` file format.

## execution
### Remarks on I/O
Please make sure to follow exactly the instructions given above.
The user must provide as an input:
- search result files deployed in `INPUT/search_results`
- quantification results (without header!) in `INPUT/quantitation_results`
- protein/peptide sequences in `INPUT/sequences`

As an output, you will get most importantly:
- `finalKinetics.csv` in the `OUTPUT/protein_name/` subfolder containing information about peptide sequence, mean intensities over technical replicate at each time point, product type and positions, among others.
- `finalKinetics.pdf` in the `OUTPUT/protein_name/PLOTS/` subfolder containing a visualisation of all kinetics sorted decreasingly by intensity.
- `KineticsDB.csv` file in the `OUTPUT/` folder containing the quantification results of **all proteins that were processed simultaneously** sorted by peptide/biological replicate

### Snakemake and Conda
qiSPI relies on [Conda](https://docs.conda.io/en/latest/) and Snakemake.
In order to install Conda, click on this [link](https://docs.conda.io/en/latest/miniconda.html) and follow the installation guidelines for your respective operating system.  
After installing Conda, you need to install Snakemake. The Snakemake installation procedure is described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Briefly, open the terminal on your computer and paste the following lines sequentially:  
`conda install -n base -c conda-forge mamba`  
`conda activate base`  
`mamba create -c conda-forge -c bioconda -n snakemake snakemake`  
Additionally, you might need to run `conda update conda`.

Download this repository as a .zip file (click on *Code* at the upper right corner of this repository --> Download ZIP), move the .zip file in the desired directory on your computer and unpack it.
Open the terminal in this directory and enter: `conda activate snakemake`.

### qiSPI execution
Make sure that you are in the correct directory of your terminal. After entering `pwd` into the terminal, it should display `./qiSPI` or `./qiSPI-main`.  

The pipeline can be executed by pasting `snakemake --use-conda --cores all -R parse_input` into the terminal. The progress of the pipeline execution should appear in your terminal window.
In case you have installed an older version of Conda/Snakemake and encounter an error when executing the pipeline, try executing
`snakemake --use-conda --cores all -R parse_input --conda-frontend conda`.

After your jobs finished, enter `conda deactivate` in order to terminate your Conda environment.

In the `OUTPUT/` directory, you should find subfolders containing the quantification results for each protein/polypeptide. Additionally, there should be a `KineticsDB.csv` file containing the quantification results of **all proteins that were processed simultaneously** sorted by peptide/biological replicate. In case you rerun the pipeline with different settings, make sure to re-name `KineticsDB.csv` since it will be overwritten upon re-execution.
