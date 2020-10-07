## A bunch of functions which tmod requires
library(glue)

#' Volcano plots with ggplot
ggvolcano <- function(lfc, pvals, symbol=NULL, xlim=NULL, ylim=NULL, pval.thr=.05, lfc.thr=1) {

  require(tidyverse)
  require(ggplot2)

  p.t <- pval.thr
  l.t <- lfc.thr

  foo <- data.frame(log2FoldChange=lfc, padj=pvals)
  foo <- foo %>% replace_na(list(padj=1, lfc=0))

  if(is.null(xlim)) xlim <- range(foo$log2FoldChange)
  if(is.null(ylim)) ylim <- range(foo$padj)

  if(!is.null(symbol)) {
    foo$Symbol <- symbol
  } else {
    foo$Symbol <- 1:nrow(foo)
  }

  foo$color <- ifelse(foo$log2FoldChange < 0, "#3333cc33", "#cc333333")
  foo$color[ foo$padj > p.t | abs(foo$log2FoldChange) < l.t ] <- "#33333333"
  foo$labels <- ""
  foo <- foo %>% arrange(sign(log2FoldChange), padj)
  foo$labels[1:20] <- foo$Symbol[1:20]
  foo <- foo %>% arrange(-sign(log2FoldChange), padj)
  foo$labels[1:20] <- foo$Symbol[1:20]

  log10rev_trans <-
    scales::trans_new("log10rev", function(x) -log10(x), function(x) 10^-x, scales::log_breaks(10), domain = c(1e-150, Inf))

  ggplot(foo, aes(x=log2FoldChange, y=padj, col=padj < 0.05 & abs(log2FoldChange) > 1 )) +
    geom_point(size=2, col=foo$color) +
    theme_minimal() +
    theme(legend.position="none") + 
    ylim(ylim[1], ylim[2]) +
    xlim(xlim[1], xlim[2]) +
    scale_y_continuous(trans=log10rev_trans) +
    ggrepel::geom_label_repel(aes(label=labels), col=gsub("..$", "", foo$color)) +
    xlab("Log 2 Fold Change") + ylab("FDR")


}


#' Filter and format tmod results for showing with DT::datatable()
tmod_filter_format_res <- function(res, pval.thr=.05, AUC.thr=.65) {

  if(is.null(res)) { return(res) }

  require(tidyverse)
  require(DT)

  res <- res %>% 
    mutate(Title=tmod_labels(ID, Title, max.length=Inf, nlines=1)) %>%
    dplyr::filter(adj.P.Val < pval.thr & AUC > AUC.thr) %>%
    dplyr::arrange(adj.P.Val, AUC) %>%
    dplyr::select(ID, Title, N=N1, AUC, p.value=P.Value, FDR=adj.P.Val)

  res
}

#' Break the labels into multiple lines, replace redundant words at the
#' beginning, replace underscore by space, return a named vector of new
#' labels
tmod_labels <- function(ids, labels, max.length=20, nlines=2) {
  #message(sprintf("length of ids: %d", length(ids)))
  if(length(ids) < 1) return(as.character(NULL))

  #ts <- tmodSummary(res, ...)
  #ids    <- ts$ID
  #labels <- ts$Title
  labels[ is.na(labels) ] <- ""

  labels <- gsub("_", " ", labels)
  labels.split <- strsplit(labels, '[[:space:]]')

  while(length(labels) > 1 && length(labels.split[[1]]) > 1 && length(unique(sapply(labels.split, `[`, 1))) == 1) { labels.split <- lapply(labels.split, `[`, -1) }
  lab.l <- lapply(labels.split, function(x) cumsum(nchar(x) + 1))

  labels.new <- lapply(1:length(lab.l), function(i) {
    x   <- lab.l[[i]]
    lab <- labels.split[[i]]

    n <- max.length - 1
    while(n > 0) {
      i1 <- which(x > max.length)[1]
      if(is.na(i1) || i1 == length(x)) {
        n <- 0
      } else {
        lab[i1] <- paste0(lab[i1], '\n')
        x <- x - x[i1]
        n <- n - 1
      }
    }
    return(lab)
  })

  labels.new <- sapply(labels.new, paste0, collapse=" ")
  labels.new <- gsub("\\n  *", "\\\n", labels.new)
  names(labels.new) <- ids
  return(labels.new)
}


#' This function selects modules to be shown. The modules are first cut-off
#' by minumum p-value and maximum AUC (over all comparisons in res). 
#' Then, if the max.n is specified and the number of
#' modules is larger, the modules are ordered by AUC and the
#' remaining modules are removed. If the number of modules is smaller than
#' min.n, then the modules which do not fullfill the AUC threshold, but which
#' do fullfill the p-value threshold are inserted. 
#' @param res tmod resuls (a list of data frames)
#' @param qval.thr q-value threshold (only matters if sel_by="auc" or max.n=Inf)
#' @param auc.thr AUC-value threshold (only matters if sel_by="p-value" or max.n=Inf)
#' @param max.n maximum number of modules to report
#' @return character vector with gene sets IDs. Attributes include notes
#'         about how the modules were chosen (whether vector of modules was shortened,
#'         extended or unchanged).

tmod_mod_sel <- function(res, qval.thr=.01, auc.thr=.65, max.n=25, min.n=max.n, res_sum=NULL,
  effect.col="AUC", pval.col="adj.P.Val") {

  if(is.null(max.n)) max.n <- 25
  if(is.null(min.n)) max.n <- 25
  if(is.null(auc.thr)) auc.thr  <- .65
  if(is.null(qval.thr)) qval.thr <- .01

  ret <- list(qval.thr=qval.thr, auc.thr=auc.thr, max.n=max.n, min.n=min.n)
  if(is.null(res_sum)) {
    res_sum <- tmodSummary(res, pval.col=pval.col, effect.col=effect.col)
    res_sum <- asS3(res_sum)
  }

  ret$full_summary <- res_sum

	if(is.null(res_sum) || nrow(res_sum) == 0) {
    ret$note <- "no results"
    return(ret)
  }

  ## lowest p-value for each row (gene set)
  res_sum$min_qvals <- res_sum %>% select_at(vars(starts_with("q.")))   %>% apply(1, function(x) min(x, na.rm=T))
  ## highest AUC value for each row (gene set)
  res_sum$max_aucs  <- res_sum %>% select_at(vars(starts_with(effect.col))) %>% apply(1, function(x) max(x, na.rm=T))

  ret$full_summary <- res_sum

  NN <- sum(res_sum$min_qvals < qval.thr & res_sum$max_aucs >= auc.thr)
  Nsign <- sum(res_sum$min_qvals < qval.thr)
  ret <- c(ret, list(#full_summary=res_sum, 
    results_n=NN, results_sign=Nsign))

  if(Nsign < 2) {
    ret$note <- "less than 2 significant results"
    return(ret)
  }

  if(NN < min.n) {
    # add ids which pass the q-value but not the AUC to fill up the results
    res_sum <- res_sum %>% dplyr::filter(min_qvals < qval.thr) %>% arrange(desc(max_aucs))
    if(nrow(res_sum) < min.n) min.n <- nrow(res_sum)
    ret$sel_ids <- res_sum$ID[1:min.n]
    ret$extended <- min.n - NN
    return(ret)
  }

  # we can filter, since NN >= min.n
  res_sum <- res_sum %>% dplyr::filter(max_aucs >= auc.thr & min_qvals < qval.thr) %>% arrange(desc(max_aucs))

  if(NN > max.n) {
    # remove the results with lowest AUC
    ret$sel_ids <- res_sum$ID[1:max.n]
    ret$shortened <- NN - max.n
    return(ret)
  }

  ret$sel_ids <- res_sum$ID
  return(ret)
}





#' Convert a data frame to a tmod object
#'
#' Each line corresponds to a gene.
#' All columns except for one in the file are assumed to contain module descriptions
#' @param df data frame with an arbitrary number of columns. 
#' @param gene_id_col which column holds the gene identifiers
#' @param module_id_col which column holds the module identifiers
#' @param module_title_col which column holds the module titles
#' @return a tmod object
df2tmod <- function(df, gene_id_col=ncol(df), module_id_col=1, module_title_col=2) {

  require(tmod, quietly=TRUE)

  gene_ids <- df[[gene_id_col]]
  module_ids <- df[[module_id_col]]
  m2g <- tapply(gene_ids, module_ids, list)
  df[[gene_id_col]] <- NULL
  df <- df[!duplicated(df[[module_id_col]]), ]
  colnames(df)[module_id_col] <- "ID"
  colnames(df)[module_title_col] <- "Title"

  makeTmod(modules=df, modules2genes=m2g)
}


#' Convert data from msigdbr package to tmod
#'
#' Use the Msigdb from the msigdbr package to generate a tmod object
#' @return a tmod object
msig2tmod <- function(taxon=NULL) {

  if(is.null(taxon)) stop("taxon is null")

  if(!require(msigdbr, quietly=TRUE)) stop("Package msigdbr not installed! Cannot proceed")

  organism <- "Homo sapiens"

  if(!is.null(taxon)) {
    if(!require(orthomapper, quietly=TRUE)) {
      stop("Package orthomapper not installed! Cannot proceed")
    }
    sptab <- orthomapper::speciesDBITable()
    if(!taxon %in% sptab$taxonID) {
      stop(sprintf("Taxon %s not in msigdbr!  use another taxonID (human or mouse, preferably)", taxon))
    }
    organism <- sptab[["species"]][ match(taxon, sptab[["taxonID"]]) ]
    if(!organism %in% msigdbr::msigdbr_show_species()) {
      stop(sprintf("Species %s not in msigdbr! Use another species (human or mouse, preferably)", organism))
    }
  }

  #df <- as.data.frame(msigdbr:::msigdbr_genesets)
  if(is.null(organism)) {
    warning("organism is null")
  }
  message("msigdb, reading organism: ", organism)
  df <- as.data.frame(msigdbr::msigdbr(species=organism))
  df <- df[ , c("gs_name", "gs_id", "gs_cat", "gs_subcat", "entrez_gene") ]
  colnames(df) <- c("Title", "ID", "Category", "Subcategory", "GeneID")

  df2tmod(df, gene_id_col=ncol(df), module_id_col=2, module_title_col=1)
}

#' Convert a tmod object to a data frame with one row per gene set
#' 
#' @param sep separator to join the gene ids
#' @param x a tmod object
#' @return a data frame with one row per module
tmod2df <- function(x, sep=";") {
  require(tmod, quietly=TRUE)
  g <- sapply(x$MODULES2GENES, function(xx) paste(xx, collapse=sep))
  cbind(x$MODULES, Genes=g, stringsAsFactors=FALSE)
}

#' Convert a CSV file to tmod format
#'
#' Read a CSV file, parse it and return a tmod object. The file should
#' contain one line per each gene set. Last column of the file is assumed to
#' be a character string with gene ID's of the genes in the gene set,
#' separated by `gene_sep` (by default, ";").
#' @param file path to a CSV file
#' @param gene_sep pattern separating genes (last column of the file)
#' @param header passed on to read.csv()
#' @param stringsAsFactors passed to read.csv
#' @param sep default column separator, passed to read.csv
#' @param ... all further arguments are passed to read.csv
#' @return a tmod object
csv2tmod <- function(file, gene_sep=";", header=TRUE, stringsAsFactors=FALSE, sep=",", ...) {

  df <- read.csv(file, header=header, stringsAsFactors=stringsAsFactors, sep=sep, ...)
  n <- ncol(df)

  genes <- df[ , n, drop=TRUE]

  modules <- df[ , -n]
  colnames(modules)[1] <- "ID"
  colnames(modules)[2] <- "Title"

  m2g <- strsplit(genes, split=gene_sep)
  names(m2g) <- modules[["ID"]]

  makeTmod(modules=modules, modules2genes=m2g)
}

tsv2tmod <- function(file, gene_sep=";", header=TRUE, stringsAsFactors=FALSE, sep="\t", ...) {
  csv2tmod(file=file, gene_sep=gene_sep, header=header, stringsAsFactors=stringsAsFactors, sep=sep, ...)
}

#' Convert a by-column CSV to tmod format
#' 
#' Convert a CSV in which each column corresponds to a gene set. The header
#' of the file contains module IDs.
#' @param file path to a CSV file
#' @param header passed on to read.csv(). If FALSE, automatic gene set
#'        names will be generated.
#' @param row.names passed to read.csv. NULL enforces numbered row names.
#' @param stringsAsFactors passed to read.csv
#' @param sep default column separator, passed to read.csv
#' @param ... all further arguments are passed to read.csv
#' @return a tmod object
#' @export
csv_bycol2tmod <- function(file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, sep=",", ...) {
  df <- read.csv(file, header=header, stringsAsFactors=stringsAsFactors, row.names=row.names, sep=sep, ...)

  Nm <- ncol(df)
  if(!header) { colnames(df) <- paste0("GID", 1:Nm) }
  mnames <- colnames(df)

  m2g <- lapply(1:Nm, 2, function(i) df[,i])
  names(m2g) <- mnames
  modules <- data.frame(ID=mnames, Title=mnames)

  makeTmod(modules=modules, modules2genes=m2g)
}


#' Subset a tmod object 
#'
#' Subset a tmod object using a parseable definition of subsets.
#'
#' @param subset definition of the subset
#' @return a tmod object
subset_tmod <- function(x, subset=NULL) {
  if(is.null(subset)) return(x)
  m <- x$MODULES

  s <- strsplit(subset, split=" *, *")[[1]]

  s <- lapply(s, function(x) unlist(strsplit(x, split=" *= *")))
  s <- lapply(s, function(xx) {
      xx <- gsub("^ *", "", xx)
      xx <- gsub(" *$", "", xx)
      xx
    })
  sel <- lapply(s, function(x) m[[ x[[1]] ]] == x[[2]])
  sel <- Reduce(`&`, sel)

  x[sel]
}


#' fill out missing information
.fill_missing <- function(db, i) {

  if(is.null(db$file)) 
    stop(sprintf("processing yaml file config/tmod/databases, entry %d: `file` field is mandatory!", i))
  if(is.null(db$name)) db$name <- paste0("DB", i)
  if(is.null(db$title)) db$title <- db$name

  db
}

not_implemented <- function(message) {
  stop(paste0(message, ": not implemented yet"))
}

#' Read a file and convert it to a tmod object
#'
#' Read a file and convert it to a tmod object
#'
#' # Possible database formats:
#' * **RDS**: A `tmod-class` object saved as an `RDS` file. See documentation for `tmod` on how to create these.
#' 
#' * **TSV,CSV**: A table with one row corresponding to each gene set. First
#'    column contains the module ID, second column module description ("Title").
#'    Last column contains list of gene ids separated by semicolons. The tables
#'    should have a header row.
#' 
#' * **TSVBYCOL**, **CSVBYCOL**: A table with one column corresponding to one
#'   gene set. The gene IDs are given in separate rows. Column headers serve
#'   as gene set IDs and Titles.
#' 
#' * **XML**: Database in the MSigDB XML format.
#' @param format file format, can be one of the following: CSV, RDS, XML,
#'               TSV, CSVBYCOL, XLSX
#' @return a tmod object
tmod_read_file <- function(x, format) {
  message(sprintf("tmod db: reading file %s, format %s", x, format))
  if(!file.exists(x)) stop(sprintf("tmod db: file %s does not exist", x))
  ret <- switch(format,
    CSV=csv2tmod(x),
    CSVBYCOL=csv_bycol2tmod(x),
    RDS=readRDS(x),
    XLSX=not_implemented("Reading XLSX tmod database format"),
    XML=not_implemented("Reading XML tmod database format"),
    TSV=tsv2tmod(x))

  return(ret)
}

#' Based on a config file, build the tmod databases
#' 
#' @param config configuration object – a list with values (the tmod:
#'        section from the yaml configuration file)
#' @return object containing all tmod databases
process_dbs <- function(config) {

  config_full <- config
  config <- config$tmod
  require(tmod, quietly=TRUE)

  dbs <- config$databases
  if(is.null(dbs)) {
    message("job tmod defined, but no databases configured!")
    warning("job tmod defined, but no databases configured!")
    return(NULL)
  }

  ## fill in name and title, check for the file field
  dbs <- lapply(1:length(dbs), function(i) .fill_missing(dbs[[i]], i)) 

  dbs.names <- sapply(dbs, function(x) x$name)
  msig    <- NULL
  msig_mm <- NULL
  tmod    <- NULL

  if(is.null(config$file_path)) config$file_path <- "./"
  
  # we define this function inline to avoid repeatedly reading the same
  # databses, especially msigdb which is large
  process_entry <- function(x) {

    dbobj <- NULL
    x$PROCESSED <- TRUE
    if(is.null(x$taxonID)) { 
      x$taxonID <- config_full$organism$taxon 
    }
    
    # two special keywords: msigdb and tmod define databases configured
    # from within the script
    if(x$file == "msigdb") {
      if(is.null(msig)) {
        message("reading msigdb")
        msig <<- msig2tmod(taxon=x$taxonID)
      }
      dbobj <- msig
      if(is.null(x$PrimaryID)) { x$PrimaryID <- "ENTREZID" }
    } else if(x$file == "tmod") {
      if(is.null(tmod)) {
        data("tmod", envir=environment())
        tmod <<- tmod
      }
      dbobj <- tmod
      if(is.null(x$PrimaryID)) { x$PrimaryID <- "SYMBOL" }
    # if a file is provided, a format is required
    } else if(is.null(x$format)) {
      stop("Processing tmod db configuration: if file path provided, format must not be empty")
    } else {
      x$file <- file.path(config$file_path, x$file)
      dbobj <- tmod_read_file(x$file, x$format)
    }

    # we need to make sure that the above code worked
    if(is.null(dbobj)) stop("tmod/process_dbs: dbobj is NULL, this should not happen")

    if(!is.null(x$subset)) dbobj <- subset_tmod(dbobj, x$subset)

    x$dbobj <- dbobj

    if(nrow(dbobj$MODULES) < 1L) {
      warning(sprintf("process_dbs: database %s contains no gene sets and will be removed", x$name))
      x$dbobj <- NULL
    } else {
      message(sprintf("process_dbs: processed database %s successfully", x$name))
    }
    x
  }

  res <- lapply(dbs, process_entry)
  names(res) <- dbs.names

  res <- res[ !sapply(res, function(x) is.null(x$dbobj)) ]

  message(sprintf("process_dbs: %d tmod databases read successfully", length(res)))

  return(res)
}



#' Map from entrez to selected column of an annDBI
#'
#' Map from entrez to selected column of an annDBI.
#' If `toupper(d$primaryID) %in% c("ENTREZ", "ENTREZID")`, simply return the annotation as is.
#' If not, find the column corresponding to primaryID in the
#' `d$annotationDBI` and translate the ENTREZ IDs to this column.
#' @param d list defining the database
#' @param annotation a named vector of entrez IDs
#' @return named vector with column replacing the entrez
map_primary_id <- function(d, annotation) {

  names(annotation) <- annotation

  if(is.null(d$primaryID)) {
    warning(sprintf("tmod database %s: no primary ID defined, using entrez", d$name))
    return(annotation)
  }

  if(toupper(d$primaryID) %in% c("ENTREZ", "ENTREZID")) 
      return(annotation)

  db.name <- d$annotationDBI

  if(is.null(db.name)) {
    warning(sprintf("tmod database %s: no annotation DBI defined, using entrez", d$name))
    return(annotation)
  }

  ## we are stopping here because this problem is easy to repair, so we
  ## want to throw it into the users face
  if(!require(db.name, character.only=TRUE)) 
    stop(sprintf("Cannot load package %s (not installed?)", db.name))

  message(sprintf("tmod database %s: translating from ENTREZID to %s", d$name, d$primaryID))
  dbi <- eval(parse(text=db.name))
  annotation[is.na(annotation)] <- "NA" # otherwise mapIds returns a list!
  newannot <- mapIds(dbi, annotation, column="ENTREZID", keytype=d$primaryID, multiVals="first")
  names(newannot) <- names(annotation)
  return(newannot)
}

#' Find a mapping for a defined database
#'
#' Find a mapping for a defined database
#' @param d a list; definition of a database mapping with elements "name",
#'          "primaryID", "taxonID" and "annotationDBI"
#' @param ann a data frame with columns "PrimaryID", "entrez", "symbol" etc
#' @param def.annot default annotation to be returned upon failure
#' @return a named vector
map_db <- function(config, d, ann, def.annot) {
  require(orthomapper, quietly=TRUE)
  message(sprintf("Preparing mapping information for %s", d$name))

  if(is.null(d$taxonID)) {
    warning(sprintf("tmod database %s: taxon ID field empty, using symbols and hoping for the best!", d$name))
    return(def.annot)
  }

  # simple case: both database and data set belong to the same organism
  if(d$taxonID == config$organism$taxon) {
    message(sprintf("tmod databse %s: same taxon (%s) as organism", d$name, d$taxonID))

    if(is.null(d$primaryID)) {
      warning(sprintf("tmod database %s: primaryID not defined, using symbols and hoping for the best!", d$name))
      return(def.annot)
    }

    if(!d$primaryID %in% colnames(ann)) {
      warning(sprintf("tmod database %s: primaryID (%s) not in the annotation, using symbols and hoping for the best!", 
        d$name, d$primaryID))
      return(def.annot)
    }

    annot <- ann[[ d$primaryID ]]
    names(annot) <- ann[[ "PrimaryID" ]]
    #annot <- map_primary_id(d, annot)
    return(annot)
  }

  # orthology mapping required because
  # d$taxonID != config$organism$taxon

  if(is.null(ann$ENTREZID)) {
    warning(sprintf("tmod database %s: entrez numbers missing from annotation, hoping for the best!", d$name))
    return(def.annot)
  }

  # retrieving orthologs of genes in the RNA-Seq project in the database
  # taxon; we use the ENTREZID for that
  message("Calling orthomapper::orthomap")
  annot <- orthomapper::orthomap(ann$ENTREZID, config$organism$taxon, d$taxonID)[,2]

  PID <- toupper(d$primaryID)

  if(!PID %in% c("ENTREZ", "ENTREZID")) {
    message("Calling orthomapper::entrez_annotate")
    annot <- orthomapper::entrez_annotate(annot, taxon=d$taxonID, 
      column=PID, orgdb=d$annotationDBI)[[PID]]
  }

  names(annot) <- ann[[ "PrimaryID" ]]
  return(annot)
}

#' Create a mapping for all databases
#'
#' Create a mapping for all defined tmod databases
#' @param config global configuration object 
#' @param dbs configuration object for the tmod databases
#' @param annotation_file RDS file containing the source annotation for the project
#' @param orthologs_file SQLite database 
#' @return a data frame with mapping for each database
get_mapping <- function(config, dbs, annotation_file) {

  # connect the required information
  ann <- readRDS(annotation_file)

  # default annotation
  def.annot <- ann$SYMBOL
  names(def.annot) <- ann$PrimaryID

# prepare a parameter table: taxonID, primaryID and annotationDBI
  .getfield <- function(l, field) {
    if(is.null(l[[field]])) return(NA)
    else return(l[[field]])
  }

  ## pt lists all the databases, their parameters and the ID of the mapping
  ## we do not need a mapping for *each* DB, only one for each unique
  ## combination of parameters: taxon, primaryID and annotationDBI
  pt <- data.frame(
    dbID=sapply(dbs, .getfield, "name"),
    taxonID=sapply(dbs, .getfield, "taxonID"),
    primaryID=sapply(dbs, .getfield, "primaryID"),
    annotationDBI=sapply(dbs, .getfield, "annotationDBI"), stringsAsFactors=FALSE)
  pt$ID <- apply(pt[,-1], 1, paste, collapse=".")

  ## running mapping for unique combinations of taxonID, primaryID and annotationDBI
  ids <- unique(pt$ID)
  maps <- lapply(ids, function(id) {
    x <- pt[ match(id, pt$ID), ]
    x$name <- paste(pt$dbID[ pt$ID == id ], collapse=", ")
    map_db(config, x, ann, def.annot)
  })

  names(maps) <- ids

  ## some statistics – helps to diagnose problems with mapping
  lapply(names(maps), function(n) 
    message(sprintf("tmod mapping %s: %d keys, %d NA's",
      n, length(maps[[n]]), sum(is.na(maps[[n]])))))

  db2map <- pt$ID
  names(db2map) <- pt$dbID

  ret <- list(maps=maps,
    pt=pt,
    dbs=db2map)
  return(ret)
}


#' Calculate 95% CI and the MSD statistic from a DE object
get_msd <- function(de, ci=.95) {

  lfc <- de[["log2FoldChange"]]
  se  <- de[["lfcSE"]]

  qq <- (1 + ci)/2

}

#' Create an ordered gene list for a given contrast
#' 
#' The sort keys are a character string of sort keys joined by commas, 
#' e.g. 'pval,pval_n,pval_p'. Sort keys are defined in the database object,
#' config and the default parameter; the first non-null value from these three
#' sources is used.
#' @param de differential expression data frame (form DESeq2)
#' @param config main sea-snap config object (a list)
#' @return a list with one element for each sort key
#' @export
get_ordered_genelist <- function(de, config, default="pval") {

  ## take the first non-NULL value
  sort_by  <- c(config$tmod$sort_by, default)[1]

  sort_by <- strsplit(sort_by, " *, *")[[1]]
  gene.ids <- rownames(de)

  ret <- lapply(sort_by, function(sby) {

    message(glue("Sorting by {sby}"))
    pval <- de[["pvalue"]]
    lfc  <- de[["log2FoldChange"]]

    sel <- NULL

    # suffix _p or _n means we have to sort only the positive lfc or
    # negative lfc, respectively
    
    # for enrichment in positively regulated genes, we replace all
    # negatively regulated genes by NA
    if(grepl("_p", sby)) {
      sel <- lfc < 0
      sby <- gsub("_p", "", sby)
    }

    # for enrichment in negatively regulated genes, we replace all
    # positively regulated genes by NA
    if(grepl("_n", sby)) {
      sel <- lfc > 0
      sby <- gsub("_n", "", sby)
    }

    if(!is.null(sel)) {
      pval[sel] <- NA
      lfc[sel]  <- NA
    }

    ret <- switch(sby,
      pval=order(pval),
      lfc=order(abs(lfc), decreasing=TRUE),
      stop(glue("incorrect sort_by value: {sby}")))

    # if sorting a part of the data, the rest should be resampled
    if(!is.null(sel)) {
      resamp <- ret %in% which(sel)
      ret[resamp] <- sample(ret[resamp])
    }
    gene.ids[ret]
  })

  names(ret) <- sort_by
  ret
}


#' Wrapper around tmod 
#'
#' Wrapper around tmod 
#' @param db object holding a tmod database in db$dbobj
#' @param genelists a list of character vectors; each character vector is a sorted gene list mapped to the IDs of the respective database
#' @param config object holding the global configuration
#' @param de results of differential gene expression analysis
#' @param db.map mapping for the databases
#' @return a list containing tmod results for each gene list: data frame – raw tmod
#'         results, no filtering, no sorting
#' @export
run_tmod <- function(db, config, genelists, db.map) {
  require(tmod)
  message(sprintf("Running tmod with db %s", db$name))

  name <- db$name

  gl <- genelists[[name]]

  tmod.func <- "tmodCERNOtest"
  
  res <- lapply(gl, function(gene.ids) {
    #gene.ids.db <- mapping[gene.ids]
    do.call(tmod.func, list(gene.ids, qval=Inf, order.by="n", mset=db$dbobj))
  })
  return(res)
}

#' Generate mapped ordered lists of genes for enrichment analysis
#'
#' For a given database, map all genelists which correspond to that
#' database with the respective database IDs.
#'
#' The ordered_genelist parameter
#' @param db database object from configuration
#' @param ordered genelist a list with ordered gene lists
#' @param db.map mapping object for databases
#' @return a list of named ordered character vectors. Names are the
#' original IDs, values are the IDs in the database.
map_genelists <- function(db, ordered_genelist, db.map) {

  name <- db$name
  mapping.id <- db.map$dbs[[ name ]]
  mapping <- db.map$maps[[ mapping.id ]]

  ret <- lapply(ordered_genelist, function(gene_ids) {
    ret <- mapping[gene_ids]
    names(ret) <- gene_ids
    return(ret)
  })

  return(ret)
}



#' Reformat tmod results
#'
#' Reformat tmod results
#' @param res results as returned from tmod
#' @param pval.thr adjusted p value threshold
#' @param auc.thr AUC threshold
#' @return a dataframe filtered and ordered by p-value
reformat_res <- function(x, pval.thr=.05, auc.thr=.55) {
  if(is.null(x)) return(x)
  x <- x[ order(x[["P.Value"]]), ]
  x <- x[ x$adj.P.Val < pval.thr & x$AUC > auc.thr, , drop=FALSE]
  x 
}
