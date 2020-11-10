[Back](../README.md) to main doc.

---

Run single cell analysis
-----------

### 1) download files (from SODAR)

First, make sure you have the following files accessible (or download them from SODAR):

- fastq files
- ISA-tab files
- feature reference (if some samples used feature barcoding)

Note: The ISA-tab assay file should contain the columns **Parameter Value[Library name mRNA]** and **Parameter Value[Library name sample tag]**.
`Library name sample tag` is only filled, when feature barcoding was used for a sample.
Samples with and without feature barcodes can be mixed.

### 2) edit the config file

After creating a [working directory](../README.md#running-the-pipeline)

```
path/to/sea-snap working_dir
```

edit the config file:

```
vim sc_config.yaml
```

Set the **in_path_pattern** to the fastq files and the **transcriptome** as well as **gtf** files for [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) that are in the reference package.
If the experiment used feature barcoding, set the path to the [**Feature Reference CSV**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

### 3) create a sample info file

This extracts information from the fastq file paths about the samples using the [in_path_pattern](prepare_input.md).

The ISA-tab assay file can also be used to extract meta-information about the samples.
It should contain columns **Parameter Value[Library name mRNA]** and **Parameter Value[Library name sample tag]**.

```
./sea-snap sample_info --from sodar --input <a_isa_assay_file> --config_files sc_config.yaml
```

### 4) run the pipeline

```
./sea-snap sc --slurm c
```

### 5) create a jupyter notebook from snippets

(this is under development; feel free to edit and add snippets)

```
./sea-snap sc l create_ipynb
```