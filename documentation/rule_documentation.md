This is the documentation for the rules. Rather then give specific file
names (or patterns) as input and output, the goal is to describe *what* is
needed as input and output, and what type of data is produced..

[[_TOC_]]

# rule: `annotation`

 * **Short description:** use the PrimaryID of the project to add annotation from
organism genome annotation databases such as org.Hs.eg.db.

 * **Context:** run once 

 * **Input:** PrimaryID, read from rownames of the `DESeq2/out/deseq2.rds` object

 * **Output:** data frame stored in annotations, with columns PrimaryID, ENSEMBL,
SYMBOL, ENTREZID, REFSEQ and GENENAME: `annotation/out/annotation.{rds,csv}`.

 * **Requirements:** orthomapper

 * **Depends on rules:** `DESeq2`

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

 * **Context:** run once for each contrast

 * **Input:** contrast definitions, `DESeq2/out/deseq2.rds` file.
 * **Output:** DataFrame object with results for each gene.

 * **Depends on rules:** `DESeq2`

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


