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

## execution
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

The pipeline can be executed by pasting `snakemake --use-conda --cores all -R parse_input` into the terminal. The progress of the pipeline execution should appear in your terminal window.
In case you have installed an older version of Conda/Snakemake and encounter an error when executing the pipeline, try executing
`snakemake --use-conda --cores all -R parse_input --conda-frontend conda`.

After your jobs finished, enter `conda deactivate` in order to terminate your Conda environment.
