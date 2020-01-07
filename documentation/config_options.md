[Back](../README.md) to main doc.

---

Config options
--------------

---

### general

|     Config keyword          |        Description                                                                |
| --------------------------- | --------------------------------------------------------------------------------- |
| **`organism_defaults`**     | Name of a file in `SeA-SnaP/defaults/` where default values for an organism are defined |
|                             |                                                                                   |
| **`organism:`**             | **(overwrite `organism_defaults`)**                                               |
| \|---`name`                 | name of the organism, e.g. "human"                                                |
| \|---`genus`                | name of the genus, e.g. "Homo Sapiens"                                            |
| \|---`taxon`                | taxon number, e.g. 9606                                                           |
| \|---`files:`               |                                                                                   |
| ...\|---`genome`            | path to a genome file (.fa); required                                             |
| ...\|---`transcriptome`     | path to a transcriptome file; if empty "" it will be automatically generated      |
| ...\|---`gtf`               | path to a gtf file with genome annotation; required                               |
| ...\|---`bed`               | path to a bed file; only required for `infer_experiment`                          |
| ...\|---`seqc_gtf`          | path to a gtf file; only required for `qc`                                        |
| \|---`star_index`           | path to a folder with indices for STAR; auto-generated if empty                   |
| \|---`salmon_index`         | path to a folder with indices for Salmon; auto-generated if empty                 |
| \|---`R:`                   |                                                                                   |
| ...\|---`annotations`       | annotation string for R's `AnnotationDbi`, e.g. "org.Hs.eg.db"                    |
|                             |                                                                                   |
|**`pipeline_param:`**        | **(general pipeline settings)**                                                   |
| \|---`out_path_pattern`     | path pattern for output files. Wildcards can be used inside braces `{...}`.<br>Available wildcards: `{step}`, `{extension}` and<br>*---mapping*: `{sample}`, `{mate}`, `{batch}`, `{flowcell}`, `{lane}`, `{library}`<br>*---DE*: `{contrast}`, `{mapping}`<br>*default mapping: `mapping/{step}/{sample}.{mate}/out/{step}.{sample}.{mate}.{extension}`*<br>*default DE: `DE/{contrast}/{step}/out/{step}.{contrast}.{extension}`* |
| \|---`log_path_pattern`     | path pattern for log files. Wildcards can be used inside braces `{...}`.<br>Available wildcards: `{step}`, `{extension}` and<br>*---mapping*: `{sample}`, `{mate}`, `{batch}`, `{flowcell}`, `{lane}`, `{library}`<br>*---DE*: `{contrast}`, `{mapping}`<br>*default mapping: `mapping/{step}/{sample}.{mate}/report/{step}.{sample}.{mate}.{extension}`*<br>*default DE: `DE/{contrast}/{step}/report/{step}.{contrast}.{extension}`* |
| \|---`in_path_pattern`      | path pattern for input files. Wildcards can be used inside braces `{...}`.<br>Available wildcards mapping: `{sample}`, `{mate}`, `{batch}`, `{flowcell}`, `{lane}`, `{library}`<br>Available wildcards DE: same as `out_path_pattern` for mapping<br>*default mapping: `../input/{sample}/{sample}.{mate}`*<br>*default DE: `mapping/{step}/{sample}.{mate}/out/{step}.{sample}.{mate}.{extension}`* |
| \|---`report_snippets`      | directory containing Rmd snippets; *default: `SeA-SnaP/report/`*                  |
| \|---`input_choice:`        | (set choices for the `choose_input()` path handler method)                        |
| ...\|---`mapping`           | For DE: list of rules to use as an input for `DESeq2`; first entry used if no wildcard `{mapping}` was set in the `out_path_pattern`<br>Options: `"import_gene_counts"` for input from STARs gene counts, `import_featurecounts` for input from FeatureCounts and `"import_sf"` for input from Salmon sf files |

### mapping pipeline

|     Config keyword          |        Description                                                                |
| --------------------------- | --------------------------------------------------------------------------------- |
|**`pipeline_param:`**        | **(general pipeline settings)**                                                   |
| \|---`mapping_results`      | list of options which algorithm to use. Options: `"salmon-transcript_counts"` for Salmon. `"star-gene_counts"` for STAR. |
| \|---`QC_results`           | list of options which QC steps to perform. Options: `"fastqc"`, `"dupradar"`, `"infer_experiment"` |
|                             |                                                                                   |
| **`rule_options:`**         | **(set parameters for rules)**                                                    |
| \|---`star:`                |                                                                                   |
| ...\|---`cmd_opt`           | a string with additional command line options for STAR                            |
| ...\|---`trim`              | trim fastq files? Options: "yes" or "no"                                          |
| \|---`star_index:`          |                                                                                   |
| ...\|---`cmd_opt`           | a string with additional command line options for STAR index generation; e.g. set "--sjdbOverhang \<read len -1\>" |
| \|---`salmon:`              |                                                                                   |
| ...\|---`cmd_opt`           | a string with additional command line options for Salmon                          |
| ...\|---`trim`              | trim fastq files? Options: "yes" or "no"                                          |
| \|---`salmon_index:`        |                                                                                   |
| ...\|---`cmd_opt`           | a string with additional command line options for Salmon index generation         |

### DE pipeline

|     Config keyword          |        Description                                                                |
| --------------------------- | --------------------------------------------------------------------------------- |
|**`experiment:`**            | **(settings about the experiment)**                                               |
| \|---`covariate_file:`      |                                                                                   |
| ...\|---`star`              | covariate file to use for input from STAR; *default: covariate_file.txt*          |
| ...\|---`salmon`            | covariate file to use for input from Salmon; *default: covariate_file.txt*        |
| \|---`design_formula`       | design formula to be used by `DESeq2`; *default: "~ group"*                       |
| \|---`columns:`             | (define level order for specific columns; *default: empty*)                       |
| ...\|--- \<column name\>    | list with level names of the column in the wished order (first level will be the reference by default) |
|                             |                                                                                   |
|**`filters:`**               | **(set filters for input data)**                                                  |
| \|---`low_counts`           | exclude genes with counts lower than <x>; *default: 0*                            |
| \|---`experiment_blacklist` | exclude certain entries of the covariate file from analysis<br>given as a dictionary of the form: {<column name>: [<level name>, ...]}<br>e.g. exclude samples: {"group": ["sample1", "sample1"]}<br>*default: {}* |
| \|---`experiment_whitelist` | only allow certain entries of the covariate file for analysis<br>given as a dictionary of the form: {<column name>: [<level name>, ...]} |
|                             |                                                                                   |
|**`contrasts:`**             | **(define which contrasts to produce and how)**                                   |
| \|---`contrast_list:`       | List of contrast definitions (see following); *default: empty*                    |
| ...\|---(list entry)        |                                                                                   |
| ......\|---`title`          | name of contrast, e.g. "nonclassical vs classical"                                |
| ......\|---`ratio`          |                                                                                   |
| .........\|---`column`      | column in covariate file, e.g. "condition"                                        |
| .........\|---`numerator`   | numerator of the contrast (a level of `column`), e.g. "nonclassical"              |
| .........\|---`denominator` | denominator of the contrast (a level of `column`), e.g. "classical"               |
| ......\|---`coef`           | alt. to `ratio`; the coefficient of DESeq2 results, e.g. "condition_classical_vs_nonclassical" |
| ......\|---`vector`         | alt. to `ratio`; a list with entries corresponding to columns in the design matrix, defining the linear combination, e.g. [1,1,0,-1,0] |
| ......\|---`goseq`          | whether to run GO and KEGG enrichment analysis with `goseq`; "true" or "false"    |
| ......\|---`cluster_profiler` | options for running MSigDB GSEA or ORA with `cluster_profiler`                  |
| .........\|---`run`         | whether to run `cluster_profiler`; "true" or "false"                              |
| .........\|---`category`    | list of MSigDb categories to test; e.g. ["H","C1","C2"]; if `category` is not set, tests all categories |
| .........\|---`type`        | whether to run gene set enrichment analysis or overrepresentation analysis; options: "gsea" or "ora"; default: "gsea" |
| ......\|---`...`            | any key from `defaults` (overwrite them for this contrast)                        |
| \|---`defaults:`            |                                                                                   |
| ...\|---`max_p_adj`         | FDR cutoff 'alpha' for DESeq2's results function; *default: 0.1*                  |
| ...\|---`ranking_by`        | rank results by (column in results table): `log2FoldChange` for log2 fold change (the effect size estimate), `pvalue` for the p-value, `padj` for the multiple testing corrected p-value |
| ...\|---`ranking_order`       | R expression for ordering (with `ranking_by`) with 'x' as input, e.g. "-abs(x)"  |
| ...\|---`results_parameters:` |                                                                                  |
| ......\|---`lfcThreshold`   | test for log fold change higher than <x>; *default: 0*                             |
| ......\|---`altHypothesis`  | alternative hypothesis of the test<br>Options: `greater`, `less`, `greaterAbs`, `lessAbs` (see [results](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results))<br>*default: greaterAbs* |
| ......\|---`independentFiltering` | perform independent filtering; "yes" or "no"                                 |
| ...\|---`lfcShrink_parameters:` |                                                                                |
| ......\|---`type`           | algorithm to use for log fold change shrinkage; Options: `none`, `apeglm`, `ashr`, `normal`<br>*default: "none"* |
| ...\|---`GO:`               |                                                                                    |
| ......\|---`fdr_threshold`  | FDR threshold to determine which results to use for functional annotation; *default: 0.1* |
|                             |                                                                                    |
|**`report:`**                | **(define which snippets to include in the report)**                               |
| \|---`report_snippets`      | List of report snippets (Rmd files) in the `report/` directory. Snippets will be appended in the order defined in this list (see section [Adding Rmd Snippets](#adding-rmd-snippets)) |
| \|---`defaults:`            |                                                                                    |
| ...\|---`contrast`          | default list of report snippets (Rmd files) added to each contrast in the report   |

---

[Back](../README.md) to main doc.