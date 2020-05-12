require(magrittr)

formatResultsForFctAnn <- function(
    res, 
    fromID=c("GENCODE", "ENSEMBL", "ENTREZID", "SYMBOL"), 
    toID=c("ENTREZID", "ENSEMBL", "SYMBOL"),
    annotation=org.Hs.eg.db::org.Hs.eg.db,
    multiVals=c("first", "filter", "maxBase", "maxAbsLFC", "minP", "random"),
    na.rm=TRUE
) {
    fromID <- match.arg(fromID, c("GENCODE", "ENSEMBL", "ENTREZID", "SYMBOL"))
    toId   <- match.arg(toID, c("ENSEMBL", "ENTREZID", "SYMBOL"))
    if (missing(res) || is.null(res) ||
        !(is(res, "DESeqResults") || is(res, "DataFrame") || is.data.frame(res)) ||
        nrow(res)<1 || !all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res)))
        stop("Missing or illegal results table")
    
    df <- as.data.frame(res)
    if (!("gene_id" %in% colnames(df))) {
        if (is.null(rownames(df)))
            stop("No gene IDs in results table")
        df$gene_id <- rownames(df)
    }
    if(any(is.na(df$gene_id) | duplicated(df$gene_id)))
        stop("Duplicated gene IDs in results table")
    df$gene_id <- as.character(df$gene_id)
    
    if (na.rm) {
        df <- df[!is.na(df$log2FoldChange) & !is.na(df$pvalue) & !is.na(df$padj),]
        if (nrow(df) == 0)
            stop("No result left after filtering")
    }
    
    gencode <- (fromID == "GENCODE")
    if (gencode) {
        df$gene_id <- sub("\\.[0-9]+$", "", df$gene_id)
        fromID <- "ENSEMBL"
    }
    
    if (fromID != toID) {
        m <- AnnotationDbi::select(annotation, keys=df$gene_id, keytype=fromID, columns=c(fromID, toID))
        stopifnot(all(m[[fromID]] %in% df$gene_id))
        m <- m[!is.na(m[[toID]]),]
    
        m <- split(match(m[[fromID]], df$gene_id), m[[toID]])
        if (multiVals=="filter") m <- m[lengths(m)==1]
        m1 <- unlist(m[lengths(m)==1])
        m2 <- m[lengths(m)>1]
        if (length(m2) > 0) {
            if (multiVals=="first")     m2 <- sapply(m2, function(i) i[1])
            if (multiVals=="maxBase")   m2 <- sapply(m2, function(i) i[which.max(df$baseMean[i])])
            if (multiVals=="minP")      m2 <- sapply(m2, function(i) i[which.min(df$pvalue[i])])
            if (multiVals=="maxAbsLFC") m2 <- sapply(m2, function(i) i[which.max(df$log2FoldChange[i])])
            if (multiVals=="random")    m2 <- sapply(m2, function(i) i[sample.int(length(1), 1)])
        }
    
        m <- sort(c(m1, m2))
        df <- df[m,]
        df[[toID]] <- names(m)
    }

    if (is(res, "DataFrame")) {
        df <- S4Vectors::DataFrame(df)
        S4Vectors::metadata(df) <- S4Vectors::metadata(res)
        tmp <- S4Vectors::DataFrame(
            type=c("gene_id", toID),
            description=c(
                sprintf("Original IDs: %s", ifelse(gencode, "GENCODE", fromID)),
                sprintf("IDs mapping using: %s", multiVals)
            ),
            row.names=c("gene_id", toID)
        )
        tmp <- tmp[!(rownames(tmp) %in% colnames(res)),]
        if (nrow(tmp) > 0) {
            tmp <- rbind(S4Vectors::mcols(res), tmp)
            tmp <- tmp[colnames(df),]
            S4Vectors::mcols(df) <- tmp
        }
    }
    if (is(res, "DESeqResults")) df <- DESeq2::DESeqResults(df, priorInfo=DESeq2::priorInfo(res))
    df
}

goAnnotToTable <- function(
    gl=NULL, keytype="ENTREZID", 
    annotation=org.Hs.eg.db::org.Hs.eg.db,
    category=c("BP", "CC", "MF")
) {
    if (is.null(gl)) gl <- AnnotationDbi::keys(annotation, keytype=keytype)
    
    tmp <- AnnotationDbi::select(annotation, keys=gl, keytype=keytype, columns=c(keytype, "GOALL"))
    tmp <- tmp[tmp$ONTOLOGYALL %in% category,]
    tmp <- unique(tmp[,c(keytype, "GOALL")])
    colnames(tmp) <- c("GENE", "TERM")
    
    query <- "SELECT go_id, term, ontology, definition FROM go_term"
    query <- paste(query, "WHERE", paste(sprintf("ontology='%s'", category), collapse=" OR "))
    go <- DBI::dbGetQuery(GO.db::GO_dbconn(), query)
    colnames(go) <- c("TERM", "NAME", "CATEGORY", "definition")
    
    go <- go[go$TERM %in% tmp$TERM,]
    
    list(TERM2GENE=tmp[,c("TERM", "GENE")], TERM2NAME=go)
}

msigdbAnnotToTable <- function(
    gl=NULL, keytype="ENTREZID",
    annotation=org.Hs.eg.db::org.Hs.eg.db,
    category=c("H", "C2")
) {
    tmp <- msigdbr::msigdbr(species=AnnotationDbi::species(annotation)) %>%
        dplyr::filter(gs_cat %in% category) %>%
        dplyr::select(TERM=gs_id, GENE=entrez_gene, NAME=gs_name, CATEGORY=gs_cat, gs_subcat, sources) %>%
        dplyr::distinct() %>%
        as.data.frame()
    
    list(
        TERM2GENE=tmp[,c("TERM", "GENE")], 
        TERM2NAME=unique(tmp[,c("TERM", "NAME", "CATEGORY", "gs_subcat", "sources")])
    )
}

userAnnotToTable <- function(
    filename, 
    fromID=c("GENCODE", "ENSEMBL", "ENTREZID", "SYMBOL"), 
    toID=c("ENTREZID", "ENSEMBL", "SYMBOL"),
    annotation=org.Hs.eg.db::org.Hs.eg.db
) {
    fromID <- match.arg(fromID, c("GENCODE", "ENSEMBL", "ENTREZID", "SYMBOL"))
    toId   <- match.arg(toID, c("ENSEMBL", "ENTREZID", "SYMBOL"))
    if (missing(filename) || is.null(filename) || !is.character(filename) ||
        length(filename)!=1 || any(is.na(filename) | filename=="")) 
        stop("Missing or illegal filename")
    df <- read.table(filename, sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE)
    if (!all(c(fromID, "TERM", "NAME", "CATEGORY") %in% colnames(df)) || nrow(df)<1)
        stop("Missing or illegal TERM2GENE table ", filename)
    
    # Remove ENSEMBL version numbers
    if (fromID == "GENCODE") {
	df$fromID <- sub("\\.[0-9]+$", "", df$fromID)
    	fromID <- "ENSEMBL"
    }

    # Create preliminary TERM2GENE table
    TERM2GENE <- unique(na.omit(df[,c("TERM", fromID)]))
    TERM2GENE <- TERM2GENE[TERM2GENE[["TERM"]]!="" & TERM2GENE[[fromID]]!="",]

    # Convert gene ids
    if (fromID != toID) {
        mappings <- AnnotationDbi::select(
            annotation,
            keys=unique(as.character(na.omit(df[[fromID]]))),
            keytype=fromID, columns=c(fromID, toID)
        )
        mappings <- mappings[!is.na(mappings[[fromID]]) & mappings[[fromID]]!="" & !is.na(mappings[[toID]]) & mappings[[toID]]!="",c(fromID, toID)]
	colnames(mappings)[2] <- "GENE"
        TERM2GENE <- TERM2GENE %>% dplyr::left_join(mappings, by=fromID)
	TERM2GENE <- TERM2GENE[,c("TERM", "GENE", fromID)]
	TERM2GENE <- unique(na.omit(TERM2GENE))
	TERM2GENE <- TERM2GENE[TERM2GENE[["TERM"]]!="" & TERM2GENE[["GENE"]]!="" & TERM2GENE[[fromID]]!="",]
    } else {
        TERM2GENE <- TERM2GENE[,c("TERM", fromID, fromID)]
    	colnames(TERM2GENE)[2] <- "GENE"
    }
    
    # Create TERM2NAME
    TERM2NAME <- unique(na.omit(df[,c("TERM", "NAME", "CATEGORY")]))
    if (any(duplicated(TERM2NAME)))
	stop("Duplicated term definition")

    list(TERM2GENE=TERM2GENE, TERM2NAME=TERM2NAME)
}
