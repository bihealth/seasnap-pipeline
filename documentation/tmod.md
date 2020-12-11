[[_TOC_]]


# Quick start to using tmod 

 1. Copy the `tmod` section from file
    `sea-snap/defaults/DE_config_defaults.yaml` into your `DE_config.yaml`
    file in the directory from which you run sea-snap. Remove the tmod
    databases you don't want to test gene enrichments for.
 2. Under each contrast for which you want to run tmod (section `contrasts:
    contrasts_list` in the `DE_config.yaml` file), add the following line:

     ```
     tmod: true
     ```

 2. In the `report: report_snippets` section of `DE_config.yaml`, enter the
    snippet 

      ```
      - Functional:
         - tmod.Rmd
      ```

 3. Optionally, copy the section `report: snippet_parameters: tmod` from
    the file `sea-snap/defaults/DE_config_defaults.yaml` into the file
    `DE_config.yaml`.
 4. Run the DE pipeline.

That's it, this will run tmod and include the relevant section in the
report.

# Running tmod in the pipeline

## Basics

To see results of gene set enrichment using tmod, you need to do two
things:

 1. Specify that tmod should be run in the `DE_config.yaml` file
 2. Specify that the tmod report snippet should be included in the report
    Rmd file.

## General configuration

tmod configuration is a section of the `DE_config.yaml` file. The following
parameters can be configured:

|Parameter      |Default value  |Explanation                                                      |
|------         |-------        |---------------------                                            |
|`sort_by`      |pval           |Sorting key: what to sort the genes by. Possible values: pval, lfc, optionally with the suffix \_p (only positively regulated) or \_n (only negatively regulated). Sorting keys may be combined using comma (,). If more than one key is provided, the enrichments will be calculating for each sorting key.|
|`tmod_db_path` | `./`          |Path prefix to where the database files are stored. By default paths will be relative to the main directory.    |
|`databases`    | see below     |Specification of the databases to be used (see below)            |

Notes:

 * If you set `tmod_db_path`, then this path will be prefixed to all files explicitely mentioned in the tmod block.

## Snippet configuration

There are several parameters that control how the results of tmod
enrichment analyses are displayed. First, there are parameters that conrol
the output in the individual contrast sections.
These are set up in the `DE_config.yaml`
file, section `snippet_parameters: tmod_contrast:`

|Parameter                        |Default value|Explanation |
|-----------                      |--------     |---------   |
|`res_auc_thr`|0.65  |Minimum AUC value to be reported in the results table|
|`res_pval_thr`|0.01 |Maximum p-value to be reported in the results table|

If you wish to see more results, increase `res_pval_thr` and decrease
`res_auc_thr`.

The enrichment results for all contrasts are summarised in the functional
analysis section.  These are set up in the `DE_config.yaml`
file, section `snippet_parameters: tmod:`

|Parameter                        |Default value|Explanation |
|-----------                      |--------     |---------   |
|`fig_qval_max`| 0.01 |Maximum q-value to be shown on the panel plots|
|`fig_auc_min` | 0.55 | Minimum AUC value to be shown on the panel plots|
|`fig_n_max`   | 35   |Maximum number of gene sets to be shown on the panel plots|
|`fig_n_min`   | 10   |Minimum number of gene sets to be shown on the panel plots|
|`n_evid`      | 5    | Number of top gene sets to be shown on the evidence plots|


# Databases in tmod

## Introduction

To run tmod, you need a database, which defines mapping between genes and
gene sets. 

A tmod database is a set of defined gene sets. For example, gene ontology
(GO) is a database defining sets of genes called "GO categories".

By default, tmod includes a database based on co-expression profiles from
human blood samples, focused on immune respones.
Another database worth noting is the MsigDB.
These two databases are available without much ado in sea-snap, you do not
need to include additional files.

Furthermore, you can provide a number of databases in one of the accepted
formats, including RDS, CSV files and more. The tmod facility in sea-snap
will import these files and run the analysis accordingly.

Each database, including the default ones (tmod and msigdb) will require
some kind of a minimal setup, described further in the document.


## Standard pipeline databases

By default, a number of gene set databases will be pre-defined for you. These
will correspond to the built-in databases tmod and various subsets of the
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) database, including
gene ontologies etc. 

The only issue you need to know if you want to run tmod with the default
db's is that the default databases are defined using human gene symbols
([HGNC](https://www.genenames.org/) symbols) or human Entrez IDs. So, what
if you have another organism, e.g. mouse?

For msigdb, we use the msigdbr package that contains mappings between human
IDs of the MSigDB and a number of other organisms (notably including
mouse). To use these mappings, simply make sure that there is no `taxonID` in the
`tmod: databases: ...` section. Then, the taxonID from 
the `DE_config.yaml` file, sections `organism_defaults:` / `organism:` will be
used. sea-snap will then attempt to use this taxon for getting the database
from the `msigdbr` package, and the genes will have native Entrez IDs.

Alternatively, modify the 
`tmod: databases: ...:` configuration, inserting a `taxonID:` indicating
the taxon ID of the organism that you are analysing. E.g., if your samples
are from mouse, use `taxonID: 10090`. The mappings to that organism
contained in the msigdbr package will be used.

The tmod database does not by itself contain such mappings. However, if you are
using a standard organism (eg. mouse) and if the entrez IDs and an annotation
DBI are defined for that organism in the config file (which they normally
should), you should be fine. The package orthomapper will map your IDs to
the IDs of the database. If not, you need to provide a mapping between your
IDs and the HGNC IDs yourself and need to read the [gory details](#mapping)
below. You also need to read on if you plan to use custom DBs.

Following databases are predefined in
`sea-snap/defaults/DE_config_defaults.yaml`, section `tmod`:

|Database ID            |Description                                                        |
|--------               |-----------                                                        |
|`tmod`                 |Transcriptional modules based on gene co-expression in human blood |
|`msigdb_go_bp`         |Gene ontology / Biological Process from MSigDB                     |
|`msigdb_reactome`      |Reactome pathways from MSigDB                                      |
|`msigdb_hallmark`      |Hallmark gene sets from MSigDB                                     |
|`msigdb_mir`           |MIR target gene sets from MSigDB                                   |


## Quick start for using a custom DB

...soon...


## ID mapping

The central issue with tmod and databases is mapping of the keys (genes).
Tmod databases refer to the genes by one of the many identifiers. For
example, the default databases msigdb and tmod use HGNC (human) identifiers
(symbols), e.g. ANKRD22 or TNF. Other databases may use something else
entirely.

The definition of the database in the yaml configuration file of sea-snap
must indicate how to translate the identifiers from these which are available
in the results file (like entrez, ensembl or plain "symbol") to the
identifiers used by the respective database.

## Subsetting

Some databases, like MSigDB, are very large and contain many irrelevant
entries. It makes sense to either split them (and analyse them
independently) or at least use only a portion of them. For this, use
subsetting.


## Configuration of tmod databases

Each database should have an entry in the tmod block of the yaml config file.
Each entry is a yaml block that at the minimum has the elements "name" and
"file".  Note that the database blocks are a yaml list (i.e., each entry starts
with a dash; see examples below). Also, if you define databases field in 

File is either a path to a file containing the mapping between genes and gene
sets (in various formats, see below), or one of two special keywords: `tmod` or
`msigdb`, which correspond to the "built-in" databases. 

Furthermore, a database definition block can contain a number of parameters
which define a user-readable title, description for the report, subsets to be
used (by default, all of the database will be used) and options for mapping between gene ids in sea-snap 
and the gene ids in the database.

|Parameter                |Explanation                                                                                    |
|----------               |--------------                                                                                 |
|file                     |**mandatory:** Path to db file or one of the two special keywords: `msigdb` or `tmod`          |
|format                   |RDS, CSV, TSV or XML; mandatory if `file` is a real file name; see below                       |
|name                     |short variable name for the databse (no spaces, no special characters)                         |
|title                    |Human readable title                                                                           |
|description              |Description of the database for humans                                                         |
|subset                   |Definition of a subset; [see below](#mapping)                                                  |
|mapping                  |Definition of a mapping; [see below](#mapping)                                                 |
|primaryID                |Type of the primary gene identifier used by the database                                       |
|taxonID                  |the taxonID to which the gene identifiers in the database refer to; by default, it is taken over from the `organism:` section of `DE_config.yaml`|
|annotationDBI            |AnnotationDBI to translate entrez IDs into the required primaryID; [see below](#mapping)       |
|sort\_by                 |Same as for general tmod options, use this to replace the sorting order for a particular db    |

Two databases are predefined and included in the pipeline implicitely:
MSigDB (via rmsigdb package) and tmod (via tmod package). You do not have
to provide a file name to use them, instead use `subset` 
([see below](#subsetting)) to define which
parts of these databases you want to use.

### Database formats

 * **RDS**: A `tmod-class` object saved as an `RDS` file. See documentation for `tmod` on how to create these.

 * **TSV,CSV**: A table with one row corresponding to each gene set. First
    column contains the module ID, second column module description ("Title").
    Last column contains list of gene ids separated by semicolons. The tables
    should have a header row.
 
 * **TSVBYCOL**, **CSVBYCOL**: A table with one column corresponding to one
   gene set. The gene IDs are given in separate rows. Column headers serve
   as gene set IDs and Titles.

 * **XML**: Database in the MSigDB XML format.


### Subsetting

tmod databases can have (internally) a number of fields (columns) associated
with gene sets. For example, all databases will have the fields "ID" (ID of a gene set) and "Title" (short description) 
defined. The database MSigDB ("msigdb" keyword) contains the additional fields
"Category" and "Subcategory" which can be used for subsetting.

The value of the `subset` parameter is a string with comma-separated values for
fields. The conditions are combined using the logical "and" operator, for example

```
subset="Category=C5,Subcategory=BP"
```

select gene sets from the database which are in category C5 *and* subcategory BP.


Here is an example:

```
    - name: msigdb_go_bp
      file: msigdb
      title: "GO Biological Process (MSigDB)"
      description:
        GO Biological Process definitions from
        the Molecular Signatures DB
        (https://www.gsea-msigdb.org/gsea/msigdb/).
      taxonID: 9606
      primaryID: ENTREZID
      subset: "Category=C5,Subcategory=BP"
```

The above definition is based on the pre-defined msigdb database, but
includes only a subset of its contents: namely only the gene sets from
category C5 and subcategory BP.

## Mapping

Mapping between the target (i.e. database) gene IDs and the source (i.e.
RNA-seq project) gene IDs relies on ortholog mapping and translation of entrez
IDs to the specified primaryID of the database. That is, if the target
(database) uses, for example, human gene IDs, but the sequenced samples were
from mouse, we map the source mouse gene IDs to orthologous target human IDs.

For this, the orthology information is derived from
[NCBI](https://ncbiinsights.ncbi.nlm.nih.gov/2018/02/27/gene_orthologs-file-gene-ftp/).
Furthermore, to convert Entrez gene IDs into the primaryID of the database (e.g. "SYMBOL" required
by the `tmod` database), an AnnotationDBI is required (e.g. `org.Hs.eg.db`).
*However*, since the annotation uses [orthomapper](https://cubi-gitlab.bihealth.org/january.weiner/orthomapper),
the annotation DBI should be correctly guessed for most of the common organisms for which a
BioConductor annotation package (such as `org.Hs.eg.db`) is available.
Refer to orthomapper documentation for details.

What exactly happens when databases are interrogated for mapping? There are
several options for correctly defining the mapping. Below, *source organism* denotes
the organism which was the source of the sequencing reads, defined in the `organism` field of the 
pipeline configuration file, `DE_config.yaml`:

 1. The tmod database defines a `taxonID` field which is identical to the taxon
		of the source organism, and a `primaryID` (such as `ENTREZ`) which is
		present in the annotation of the count data matrix generated by the `annotation` rule.
 
 2. The `taxonID` of the tmod database is *different* from the taxon of the source organism. In that case,
		`orthomap()` function from the
		[orthomapper](https://cubi-gitlab.bihealth.org/january.weiner/orthomapper)
		package will be used to find orthologs between the genes of the source
		organism and the genes in the target databse. Then, `entrez_annotate()` from the same package
		will be used to translate the ENTREZID's returned by `orthomap()` to the IDs defined by `primaryID`
		in the database definition.

If the `taxonID` or `primaryID` are not defined for a tmod database, or if the
`primaryID` cannot be found in the respective org.db (such as org.Hs.eg.db for
taxon 9606), then a default mapping will be used and the script will issue a
warning ("Hoping for the best!"). The default mapping is simply the gene symbols derived from the
`annotation` rule.

## Example configuration 

Here are some examples of DB configuration:

``` yaml
tmod:
  databases:
    - name: msigdb_mir
      file: msigdb
      title: "MIR targets (MSigDB)"
      description:
        MIR targets from
        the Molecular Signatures DB
        (https://www.gsea-msigdb.org/gsea/msigdb/).
      taxonID: 9606
      primaryID: entrez
      subset: "Category=C3,Subcategory=MIR"
    - name: tmod
      file: tmod
      title: "Co-expression gene sets (tmod)"
      description:
        Gene sets derived from clustering expression profiles from human blood
        collected for various immune conditions. These gene sets are included
        in the tmod package by default. Check tmod documentation for further 
        information.
      taxonID: 9606
      primaryID: SYMBOL
      annotationDBI: org.Hs.eg.db
```
    
Above, two databases were defined: `msigdb_mir`, which is a subset of the
msgidb, and `tmod`, the database which comes with the tmod package. Both
databases are based on human genes (`taxonID: 9606`), but while msigdb uses
entrez IDs (`primaryID: entrez`), tmod uses symbols (`primaryID: SYMBOL`).

## Default configuration
