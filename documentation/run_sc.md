[Back](../README.md) to main doc.

---

Run single cell analysis
-----------

### 1) edit the config file

```
vim sc_config.yaml
```

Set the **in_path_pattern** and the **transcriptome** and **gtf** files for [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) that are in the reference package.
If the experiment used feature barcoding, set the path to the [**Feature Reference CSV**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

### 2) create a sample info file

First get the ISA-tab assay file. 
That can be used to extract meta-information about the samples.
It should contain columns **Parameter Value[Library name mRNA]** and **Parameter Value[Library name sample tag]**.

```
./sea-snap sample_info --from sodar --input <a_isa_assay_file>
```

### 3) run the pipeline

```
./sea-snap sc --slurm c
```

### 4) create a jupyter notebook from snippets

(this is under development; feel free to edit and add snippets)

```
./sea-snap sc l create_ipynb
```