This is the documentation for the rules. Rather then give specific file
names (or patterns) as input and output, the goal is to describe *how* the
rule is working (what is the logic behind the code) and *what* is
needed as input and output, what type of data is produced etc.

[[_TOC_]]

# rule: `all`

 * **Short description:** The `all` rule depends on the generated
   `report.Rmd` file and on the files returned by `get_inputs_all()`,
   which have been added using `inputs.append()`

# rule: `TxDb_from_GTF`

 * **Short description:** The GTF file provided by the species definition
   is converted to a TxDB, which can be used from R.

 * **Requirements:** `GenomicFeatures`

 * **Context:** run once (contrast `all`)

# rule: `import_gene_counts`, `import_featurecounts`, `import_sf`

 * **Short description:** Import STAR / feature count / Salmon counts from
   the files listed in the covariate file. Produces a DESeq2 object saved
   as RDS data.

 * **Input:** Files defined in the covariate file, containing the gene
   or transcript counts produced by the respective method. Furthermore, it
   a model is required.
 
 * **Output:** a DESeq2 object 

 * **Requirements:** `DESeq2`

 * **Context:** run once (contrast `all`)

# rule: `export_raw_counts`

 * **Short description:** Save the collected raw counts as XLSX for
   external usage.

 * **Input:** DESeq2 object created by one of the import rules.

 * **Output:** An Excel (XLSX) file with raw counts.

 * **Requirements:** `DESeq2`, `writexl`

 * **Context:** run once (contrast `all`)

# rule: `annotation`

 * **Short description:** use the PrimaryID of the project to add annotation from
organism genome annotation databases such as org.Hs.eg.db.

 * **Input:** PrimaryID, read from rownames of the `DESeq2/out/deseq2.rds` object

 * **Output:** data frame stored in annotations, with columns PrimaryID, ENSEMBL,
SYMBOL, ENTREZID, REFSEQ and GENENAME: `annotation/out/annotation.{rds,csv}`.

 * **Requirements:** orthomapper

 * **Context:** run once (contrast `all`)

# rule: `goseq`
 
 * **Short description:** run basic gene set enrichment analysis (ORA /
   hypergeometric test) on KEGG and GO gene sets.

 * **Output:** Two RDS files, one for KEGG, one for GO, with the results
   of the analysis.

 * **Context:** run once for each contrast

# rule: `cluster_profiler`

 * **Short description:** run clusterProfiler using predefined sets of
   genes (see cluster profiler configuration).

 * **Output:** 
   * an RDS file. The object stored is a list with one element for each
     type of test to run; these elements are lists which one element for
     each gene set to run.
   * a directory (`cluster_profiler.csv/`) which contains the results of
     the tests.

 * **Depends on rules:** `contrast`

 * **Requirements:** annotation DBI package, msigdbr, clusterProfiler
 
 * **Context:** run once for each contrast


# rule: `tmod_dbs`

 * **Short description:** Based on tmod database configuration read from
   the YAML configuration file (`DE_config.yaml`), compile and store the
   gene set databases used by tmod for gene set enrichments. Furthermore,
   use the annotation to produce a mapping between feature PrimaryIDs of the
   project and the gene primaryID of the database.

 * **Context:** run once 

 * **Input:** output of the annotation rule, name of the yaml configuration
   file produced by snakemake in the `onstart:` block.
 * **Output:** 
     * RDS object: named list containing the compiled databases. Each
       element is itself a named list with details of database
       configuration as well as a tmod object which can be directly used by
       tmod. This file can be quite large.
     * RDS object: list containing two elements:
         * "maps", a list of mappings
         * "dbs", a list associating a mapping with each of the databases
           (note that since many databases can use the same mapping, there
           will be fewer mappings than databases)
 * **Requirements:** orthomapper, tmod

 * **Depends on rules:** `annotation`

# rule: `contrasts_full`

 * **Short description:** Same as the rule `contrast`, except that the
   output produced is not truncated or filtered, but instead contains the
   full list of *all* results for all genes, as this is required for a
   reasonable gene set enrichment analysis.

 * **Input:** contrast definitions, `DESeq2/out/deseq2.rds` file.

 * **Output:** DataFrame object with results for each gene.

 * **Depends on rules:** `DESeq2`

 * **Context:** run once for each contrast


# rule: `tmod`

 * **Short description:** This runs gene set enrichment for each contrast
   and each database.

 * **Context:** run once for each contrast

 * **Input:** results of the differential gene expression analysis for each
   contrast; tmod databases object; annotation.
 * **Output:** for a given contrast, an object which is a list with one
   element for each database that was tested. Furthermore, an object that
   contains the ordered list of primary identifiers of genes sorted by the
   one of the sorting options. This object is a list with an element for
   each database, each of which contains an element of 


 * **Depends on rules:** `tmod_dbs`, `contrasts_full`


# rule: `tmod_summary`

 * **Short description:** Gather the objects for all contrasts for which
   tmod was run and restructure the mapping, grouping the results by
   database rather than contrast. The rationale is that for visualization
   we want to put the different contrasts next to each other within a
   section for a given database.

 * **Context:** run once 

 * **Input:** results of the tmod analysis run in the `tmod` step.
 * **Output:** object holding (i) a list of all the results *by database*,
   and (ii) a mapping between contrast IDs used by the pipeline and human
   defined contrast names (as taken from the yaml configuration file).

 * **Depends on rules:** `tmod`


