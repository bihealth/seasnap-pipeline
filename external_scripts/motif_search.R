#' @param peaks.meta a data frame with  columns Chr, Start and End
#' @param chr character vector with chromosomes to write out
#' @param file.pattern file pattern suitable for glue; use '{chr}' to denote chromosome
#' @param out.dir output directory
write_peaks <- function(peaks.meta, chr_v, file.pattern, out.dir) {
  for(chr in chr_v) {
    fn <- glue(file.pattern)
    fc <- file(file.path(out.dir, fn))
    lines <- with(peaks.meta %>% filter(Chr == chr), glue("{Start} {End}\n"))
    writeLines(lines, fc)
    close(fc)
  }
}



#' Read all results from a dreme results directory
#' 
#' dn_motif_n is a list, with one element for each contrast / each
#' subdirectory of the motif.dir. Each of these elements itself is a list
#' containing elements `up`, `down` and `neg` which correspond to numbers of
#' peaks which were written to the respective fasta files in the
#' motif.dir/contrast directories.
#' @param motif.dir path to the directory
#' @param dn_motif_n counts of motifs (see Details)
#' @param contrasts Contrasts for which to read the results. Must be a subset of names(dn_motif_n)
read_dreme_results_all <- function(motif.dir, dn_motif_n, contrasts) {

  ret <- imap_dfr(dn_motif_n[contrasts], ~ { 
      dir <- file.path(motif.dir, .y)
      res <- read_dreme_result(dir, .x$up, .x$down, .x$neg)
      cbind(Contrast=.y, res)
  })

  return(ret)
}


#' Read all tomtom results for all contrasts
#'
#' 
read_tomtom_results_all <- function(motif.dir, jaspar.dir, contrasts) {

  ret <- map_dfr(contrasts, ~ {
      tomtom <- read_tomtom_result(file.path(motif.dir, .x), jaspar.dir)
      tomtom %>% mutate(Contrast=.x)
  })
  return(ret)
}



#' Read dreme results for a given contrast
#'
#' @param path path to the folder with dreme results
#' @param tot_up,tot_down total numbers of peaks in the Up and Down data sets
read_dreme_result <- function(path, tot_up, tot_down, tot_neg) {

  message(path)
	f_up   <- file.path(path, "TPup_dreme", "dreme.txt")
	f_down <- file.path(path, "TPdown_dreme", "dreme.txt")

	dreme_res_up   <- read_dreme(f_up)   %>% mutate(Direction="Up",   Pos = as.numeric(Pos) / tot_up)
	dreme_res_down <- read_dreme(f_down) %>% mutate(Direction="Down", Pos = as.numeric(Pos) / tot_down)

	ret <- rbind(dreme_res_up, dreme_res_down) %>% 
    mutate(Neg = as.numeric(Neg) / tot_neg, Enrichment=Pos/Neg, No=gsub("DREME-", "", No))
  ret <- ret %>%  mutate(Logo= file.path(path, sprintf("TP%s_dreme", tolower(Direction)), sprintf("m%02dnc_%s.png", as.numeric(No), Motif))) 
  return(ret)
}

#' Read tomtom results for a given contrast
read_tomtom_result <- function(path, jaspar.dir) {

	f_up   <- file.path(path, "TPup_tomtom", "tomtom.txt")
	f_down <- file.path(path, "TPdown_tomtom", "tomtom.txt")

	tomtom <- rbind(
		data.frame(read.table(f_up, comment.char="", sep="\t", header=T, stringsAsFactors=F), Direction="Up"),
		data.frame(read.table(f_down, comment.char="", sep="\t", header=T, stringsAsFactors=F), Direction="Down"))

	tomtom$Annot <- sapply(tomtom$Target.ID, function(id) {
		read_motif(file.path(jaspar.dir, paste0(id, ".meme")))[3]
  })
  tomtom <- tomtom %>% filter(q.value < .01) %>% rename(Annotation=Annot, Query=X.Query.ID)

  tomtom
}


#' Parse a tomtom results file and augment it by motif annotation
#' @param file tomtom results file (TSV)
#' @param jaspar.dir directory with jaspar motifs
#' @return data frame
read_tomtom <- function(file, jaspar.dir=NULL) {
  require(dplyr)

 res <- read.table(file, comment.char="#", sep="\t", header=TRUE, stringsAsFactors=FALSE) 
 if(!is.null(jaspar.dir)) {
   res$Annotation <- map_chr(res$Target_ID, ~ read_motif(file.path(jaspar.dir, paste0(.x, ".meme")))[3])
 }
 return(res)
}



#' Compile dreme and tomtom output
#'
motifsearch_summary <- function(dreme_res, tomtom_res, motif.dir, enr.thr=2) {

  dreme_sel <- dreme_res %>% dplyr::filter(Enrichment > enr.thr) %>% 
    mutate(Logo=sprintf("![%s](%s)", Motif, Logo)) %>% 
    dplyr::select(Contrast, Direction, Motif, Pval, Logo=Logo)

  dreme_sel$Annotation <- sapply(1:nrow(dreme_sel), function(i) {
    ann <- tomtom_res %>% dplyr::filter(Contrast == dreme_sel$Contrast[i],
    Direction == dreme_sel$Direction[i], Query == dreme_sel$Motif[i]) %>%
    pull(Annotation) %>% unique %>% paste(collapse = ", ")

  })

  dreme_sel <- dreme_sel %>% dplyr::filter(Annotation != "") %>% 
    mutate(Direction = ifelse(duplicated(paste0(Contrast, Direction)), ",,", Direction),
           Contrast  = ifelse(duplicated(Contrast), ",,", Contrast))


}


#' Prepare dreme results table
#' @param dreme_res data frame produced by read_dreme()
#' @param topn How many of the top hits to show for each direction / contrast combination
dreme_res_table <- function(dreme_res, topn=5) {

  .dreme_res <- dreme_res %>% 
    mutate(dc=paste(Direction, Contrast), i=1:n())  %>% group_by(dc) %>% top_n(-topn, i) %>% ungroup %>%
    mutate(Logo=sprintf("![%s](%s)", Motif, Logo)) %>%
    dplyr::select(Contrast, Direction, Motif, Logo, Pos:E) 
  .tmp <- paste(.dreme_res$Contrast, .dreme_res$Direction)
  .dreme_res$Contrast[ duplicated(.dreme_res$Contrast) ] <- ",,"
  .dreme_res$Direction[ duplicated(.tmp) ] <- ",,"

  return(.dreme_res)
}


#' Prepare the tomtom results table
tomtom_res_table <- function(tomtom_res) {
  ret <- tomtom_res %>% dplyr::select(Contrast, Direction, Query, Target.ID, Annotation, p.value, q.value)
  ret$Direction[ duplicated(paste0(tomtom_res$Contrast, tomtom_res$Direction)) ] <- ",,"
  ret$Contrast[ duplicated(tomtom_res$Contrast) ] <- ",,"
  ret
}


reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

#' Calculate enrichment of x in y
#' @param x a logical vector
calcE <- function(x, y) {
  if(length(x) != length(y)) stop("length x and y differ")
  b <- sum(x & y, na.rm=T)
  n <- sum(y, na.rm=T)
  B <- sum(x, na.rm=T)
  N <- length(x)

  c(b=b, n=n, B=B, N=N, E=(b/n)/(B/N))
}


myrender <- function(fn, open=T, ...) {

  fb <- gsub("\\.rmd$", "", fn, ignore.case=T)
  fn2 <- paste0(fb, "_html.rmd")
  file.copy(fn, fn2, overwrite=T)
  res <- render(fn2, output_format="html_document", envir=globalenv(), ...)
  system(sprintf("google-chrome %s", res))
  res
}

#' message + sprintf
smessage <- function(...) message(sprintf(...))


#' Static cache mechanism
#'
#' Checks for the existence of the specified file. If it exists, it is
#' loaded with readRDS. If not, expr is evaluated and the result is stored in
#' the specified path using saveRDS. In any case, the object stored is
#' returned.
#' @param path path to an RDS file
#' @param expr expression which is the "recipe" for the object
#' @param name (optional) name of the variable to report in messages
#' @param quiet if TRUE, no messages are shown
x_static_cache <- function(name, path, expr, quiet=FALSE) {

  if(file.exists(path)) {
    ret <- readRDS(path)
    if(!quiet) { 
      message(sprintf("%s: Loading existing version from file %s, created on %s", name, path, attr(ret, "creation time")))
    }
  } else {
    if(!quiet) {
      message(sprintf("%s: Generating object and saving to file %s", name, path))
    }

    ret <- expr
    attr(ret, "creation time") <- Sys.time() 
    saveRDS(ret, file=path)
  }

  return(ret)
}


#' get coverage of transcription start sites
coverageTSS <- function(reads, genes, width=3000, chromosome="chr20", plot=TRUE, strand=NULL,
  min.length=-Inf, max.length=Inf, start=TRUE, ...) {
	## ranges of reads
	reads.gr <- granges(reads)
  reads.w <- width(reads.gr)
  reads.gr <- reads.gr[ reads.w > min.length & reads.w < max.length]
	reads_cov <- coverage(reads.gr)[[chromosome]]
	n <- length(reads.gr)

	## transcription start positions
	TSreg <- flank(genes, width=width/2, start=start, both=TRUE)
	TSreg <- TSreg[ seqnames(TSreg) == chromosome ]

	pos <- cbind(start(TSreg), end(TSreg))
	g.strand <- as.character(decode(strand(TSreg)))
  ngenes <- length(TSreg)

	if(!is.null(strand)) { sel <- g.strand == strand } else { sel <- TRUE }
  print(sum(sel))

	foo.B <- sapply((1:ngenes)[sel], function(.) { catf("\r%d", .) ; reads_cov[pos[.,1]:pos[.,2]] })
  cat("\n")
	foo2.B <- sapply(foo.B, decode)

	#iqrs <- t(apply(foo2.B, 2, IQR))
	#qs <- t(apply(foo2.B, 2, function(.) quantile(., c(0.25, 0.75))))
	#toshow <- iqrs > 1
  toshow <- TRUE
	scores.B <- rowSums(foo2.B[,toshow])/ncol(foo2.B[,toshow])

	if(plot) plot(scores.B, type="l", bty="n", ...)
 	ret <- data.frame(N=1:width - width/2, score=scores.B)

	invisible(ret)
}

#' read dreme results from a single file
#' @param file file with the results
read_dreme <- function(file) {
  require(dplyr)
  con <- file(file)
  rl <- readLines(con)
  close(con)

  rl.motifs <- rl[ grep("^MOTIF", rl) ] %>% map_chr(~ strsplit(.x, split=" ")[[1]][2])
  rl.no <- rl[ grep("^MOTIF", rl) ] %>% map_chr(~ strsplit(.x, split=" ")[[1]][3])
  rl.matches <- rl[ grep("^# BEST", rl) ] %>% { gsub("^# BEST  *", "", .) } %>% { read.table(text=.) } %>% as_tibble
  colnames(rl.matches) <- c("Best", "RCWord", "Pos", "Neg", "Pval", "E")
  rl.matches <- rl.matches %>% mutate(Motif=rl.motifs, No=rl.no) %>% dplyr::select(Motif, Best:E, No)
  rl.matches
}

#' Read a motif from a jaspar file
#' @param file path to motif file
read_motif <- function(file) {
  con <- file(file)
  rl <- readLines(con)
  close(con)

  ret <- rl[ grep("^MOTIF ", rl) ] %>% map(~ strsplit(.x, split=" ")[[1]]) %>% { Reduce(rbind, .) }
  return(ret)

}


#' 
interleave <- function(vec1, vec2) {

  if(length(vec1) != length(vec2)) stop("vectors don't have equal length")
  
  as.vector(matrix(c(vec1, vec2), nrow=2, byrow=T))


}


#' Calculate and plot autocorrelation for fragment size distribution
autoCplot <- function(sizes, frag.min=50, frag.max=1000, lag.max=1000, log=TRUE, plot=TRUE, ...) {

  sizes.t <- table(sizes)
  sizes.n <- as.numeric(names(sizes.t))
  sizes.v <- as.numeric(sizes.t)

  sel <- sizes.n %in% frag.min:frag.max

  if(log) sizes.v <- log(sizes.v)
  trend <- lm(sizes.v[sel] ~ sizes.n[sel])
  sizes.r <- trend$residuals

  #sel <- frag.min:frag.max

  if(plot) {
    par(mfrow=c(1,2))
    if(log) ylab <- "log(Fragment count)"
    else ylab <- "Fragment count"
    plot(sizes.n[sel], sizes.v[sel], type="l", bty="n", xlab="Fragment size", ylab=ylab)
    abline(trend, col="red")
    par(new=TRUE)
    plot(sizes.n[sel], sizes.r, col="blue", type="l", axes=F, xlab="", ylab="", bty="n", ylim=c(min(sizes.r), 2 * max(sizes.r)))
    legend("topright", c("Fragment distribution", "Trend", "Residuals"), lty=1, col=c("black", "red", "blue"), bty="n")
  }

  ret <- acf(sizes.r, lag.max=lag.max, plot=plot, bty="n", ...)

  invisible(list(acf=ret, sizes.n=sizes.n, sizes.r=sizes.r, sizes.v=sizes.v))
}

## convenience function for testing different parameter sets
try_dba <- function(dba, summits=50, overlap=1, mset=NULL) {
  name <- sprintf("s.%d_o.%d", summits, overlap)
  print(name)
  .dba.2 <- dba.count(dba, summits=summits, minOverlap=overlap)
  .dba.2 <- dba.contrast(.dba.2, categories=DBA_TISSUE)
  .dba.2 <- dba.analyze(.dba.2)
  .dba.rep <- dba.report(.dba.2, th=1)
  .dba.rep.an <- as.GRanges(annotatePeak(.dba.rep, TxDb=txdb))
  .dba.rep.an$symbol <- genemap$symbol[ match(.dba.rep.an$geneId, genemap$EG) ]
  cerno <- tmodCERNOtest(toupper(.dba.rep.an$symbol), mset=mset)
  res <- list()
  res[[name]] <- cerno
  res
}


volcano <- function(p, lfc, lfc.thr=1, p.thr=0.01, ylim=NULL) {

  if(is.null(ylim)) ylim <- rev(range(p))
  plot(lfc, p, ylim=ylim, log="y", bty="n", pch=19, col="#33333333",
    xlab="log(Fold)", ylab="FDR")

  sel <- p < p.thr & abs(lfc) > lfc.thr
  points(lfc[sel], p[sel], pch=19, col="#cc000033")
}

c_ <- function (label, asis=TRUE) 
{
		if(asis) label <- as.character(substitute(label))
    prefix <- unlist(strsplit(label, "_"))[1]
    message(prefix)
    if (!exists("counters")) 
        counters <<- list()
    if (is.null(counters[[prefix]])) 
        counters[[prefix]] <<- c()
    if (is.na(match(label, counters[[prefix]]))) {
        counters[[prefix]] <<- c(counters[[prefix]], label)
    }
    return(match(label, counters[[prefix]]))
}

