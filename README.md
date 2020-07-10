# About
This repository is used to conduct the analysis for "The dynamics of B-cell repertoires in response to SARS-CoV-2 infection in COVID-19 patients with different disease severity."
It investigates receptor compostions depending on disease severity, expansion of BCR clonal lineages over time, sharing of BCRs among individuals, and the emergence of cross-reactivity from a SARS-CoV-2 response to a SARS response.

## Dependencies

Most dependences are listed in `env.yml`.
You can set up a [conda](https://docs.conda.io/en/latest/) environment with

```bash
$ conda env create -f env.yml
```
and activate the new environment (including all dependencies) with
```bash
$ conda activate covid-bcr
```

Other necessary software includes:
 - [R](https://www.r-project.org/)
 - [IGoR](https://github.com/qmarcou/IGoR)
 - [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/)

## Pipeline and preparing data for analyses

### Raw sequence processing

`scripts/presto_pipeline.sh` processes the data with [pRESTO](https://presto.readthedocs.io/en/stable/overview.html).
Pairs of raw sequences are assembled and filtered with a QScore of 30, V primers are masked and the C primer is cut, and finally sequences are deduplicated. For example,

```bash
bash presto_pipeline.sh -s SAMPLENAME-REPLICATE -afmc
```

where `-s SAMPLENAME-REPLICATE` specifies the sample, e.g. `-s S1-0`, and `-afmc` will (a)ssemble the paired-end reads, (f)ilter sequences for quality, (m)ask primers, and (c)ollapse sequences.

### Assemble FASTA file for annotating

Use `scripts/assemble_fastas_for_abstar.py` to annotate replicate, timepoint, severity, and patient information in headers and combine sample FASTAs into one larger FASTA for each patient.
This script takes in the directories where the data is stored and a directory to which the combined FASTAs will be saved.

```bash
python assemble_fastas_for_abstar.py --dirs dir1 dir2 dir3 --save_dir annotated_fastas_dir
```

The output FASTAs will be called PATIENTID.fasta, where PATIENTID is an integer <img src="https://render.githubusercontent.com/render/math?math=\large z \in [1,19]">.

### Annotating sequences

[abstar](https://github.com/briney/abstar) should have been installed when creating the conda environment. To annotate sequences, simply use

```bash
abstar -i /PATH/TO/INPUT_FILE.fasta -o /PATH/TO/OUTPUT_DIRECTORY -t /PATH/TO/TEMP_DIRECTORY
```

for each file or 

```bash
abstar -i /PATH/TO/INPUT_DIRECTORY -o /PATH/TO/OUTPUT_DIRECTORY -t /PATH/TO/TEMP_DIRECTORY
```

to run abstar iteratively over all files in a directory.

### Error correction, data filtering, creating lineages, and creating input for SONIA

Use `scripts/abstar_pipeline.py` to error correct and filter sequences

```bash
python abstar_pipeline.py --annotations PATH/TO/ABSTAR_OUTPUT.json PATH/TO/SAVE/FILTERED_ANNOTATIONS.json
```

create lineages

```bash
python abstar_pipeline.py --lineages PATH/TO/FILTERED_ANNOTATIONS.json PATH/TO/SAVE/LINEAGES.json
```

and create input for SONIA

```bash
python abstar_pipeline.py --sonia PATH/TO/LINEAGES.json PATH/TO/SAVE/SONIA_INPUT.csv
```

See the code or Methods in [Montague et al., Dynamics of B-cell repertoires and emergence of cross-reactive responses in COVID-19 patients with different disease severity]() for more details.

### Create input for IGoR

Use `scripts/assemble_fasta_for_igor.py` to create the FASTA file used as input for IGoR.

```bash
python assemble_fasta_for_igor.py --infiles PATH/TO/LINEAGES/* --outfile PATH/TO/SAVE/IGOR_INPUT.fasta
```

### Prepare input for expansion analysis

Use `scripts/expansion_data_analysis.py` to create the .csv file used as input for the expansion analysis.

```bash
python expansion_data_analysis.py --lineages PATH/TO/LINEAGE_FILE.json --outfile PATH/TO/SAVE/LINEAGE_COUNTS.csv
```

### Obtain CDR3 anchors and genomic references

All genomic references and CDR3 anchors for IGoR and SONIA are in `igor_input/`. If you would like to prepare your own CDR3 anchors and genomic references, use `scripts/cdr3_anchors_and_references.py`.

```bash
python cdr3_anchors_and_references.py --indir PATH/TO/ANCHORS_AND_REFERENCES_DIR
```

## Analysis

### IGoR model

Genomic references and CDR3 anchors for [IGoR](https://github.com/qmarcou/IGoR) are available in `igor_input/`. Use `scripts/run_igor.sh` to make the IGoR model.

```bash
bash run_igor.sh PATH/TO/WORKING_DIRECTORY PATH/TO/IGOR_INPUT.fasta BATCHNAME
```

### SONIA models

IGoR model output and CDR3 anchors for [SONIA](https://github.com/statbiophys/SONIA) are available in `sonia_input/`.
To assemble all the individual patient data into cohorts quickly, `cat` them together and then delete the unnecessary header lines that are not the first line.

To create a SONIA model, execute

```bash
sonia-infer --set_custom_model_VDJ PATH/TO/SONIA_INPUT.csv --sonia_model leftright --epochs 150 \
            --independent genes --seq_index 0 --v_mask_index 1 --j_mask_index 2 \
            --infile PATH/TO/SONIA_INPUT.csv -o PATH/TO/SAVE/SONIA_MODEL --lines_to_skip 1 \
            --n_gen_seqs 500000
```

If you want to add some determinism and replicability when creating a SONIA model, you can additionally use the `--seed` option in `sonia-infer` which will lead to SONIA models with nearly identical selection factors.

To generate sequences, execute

```bash
sonia-generate --set_custom_model_VDJ PATH/TO/SONIA_MODEL --sonia_model leftright --ppost \
               --N NUM_SEQS_TO_GENERATE --delimiter_out , --outfile PATH/TO/SAVE/GENERATED_SEQS.csv
```

To evaluate sequences from data, execute

```bash
sonia-evaluate --set_custom_model_VDJ PATH/TO/SONIA_MODEL --sonia_model leftright --ppost \
               --infile PATH/TO/SONIA_INPUT.csv --lines_to_skip 1 --delimiter_out , \
               --outfile PATH/TO/SAVE/DATA_EVALUATIONS.csv
```

To evaluate generated sequences, execute

```bash
sonia-evaluate --set_custom_model_VDJ PATH/TO/SONIA_MODEL --sonia_model leftright --ppost \
               --infile PATH/TO/INPUT.csv --delimiter_out , --outfile PATH/TO/SAVE/GEN_EVALUATIONS.csv
```


### Differential sequence features

Use `scripts/abstar_stats.py` to obtain statistics of gene usage, HCDR3 length, and insertion and deletion profiles of nonsingletons and progenitors in productive lineages with at least three unique sequences across all timepoints.

```bash
python abstar_stats.py --infile PATH/TO/LINEAGES.json --outfile PATH/TO/SAVE/STATISTICS.json
```

For plotting these statistics and performing ANOVA (with boxplot visualizations), see `notebooks/sequence_features_plotting.ipynb`.

### Expansion analysis and overlap with known Abs

See `notebooks/expansion_and_overlap_analysis.ipynb`.

### Sharing analysis

<img src="https://render.githubusercontent.com/render/math?math=\large P_{post}"> obtained from evaluating data sequences using inferred SONIA models.
See `sequence_features_plotting.ipynb`.
See Methods in [Montague et al., Dynamics of B-cell repertoires and emergence of cross-reactive responses in COVID-19 patients with different disease severity]() for calculating the null hypothesis for <img src="https://render.githubusercontent.com/render/math?math=\large P_{share}">.

## References

1. Montague et al., Dynamics of B-cell repertoires and emergence of cross-reactive responses in COVID-19 patients with different disease severity, (2020)
2. Vander Heiden, J.A., Yaari, G., Uduman, M., Stern, J.N., O'Connor, K.C., Hafler, D.A., Vigneault, F., and Kleinstein, S.H. (2014). pRESTO: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. Bioinformatics 30, 1930-1932.
3. Briney, B., and Burton, D.R. (2018). Massively scalable genetic analysis of antibody repertoires. bioRxiv 10.1101/447813, 447813.
4. Marcou, Q., Mora, T., and Walczak, A.M. (2018). High-throughput immune repertoire analysis with IGoR. Nat Commun 9, 561.
5. Sethna, Z., Isacchini, G., Dupic, T., Mora, T., Walczak, A.M., and Elhanati, Y. (2020). Population variability in the generation and thymic selection of T-cell repertoires. bioRxiv 10.1101/2020.01.08.899682, 2020.2001.2008.899682.

## Contact

Any issues or questions should be addressed to [me](mailto:zacander.mon@gmail.com).

## License

Free use of this analysis is granted under the terms of the GNU General Public License version 3 (GPLv3).
