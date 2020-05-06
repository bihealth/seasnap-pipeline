# Recipies

Specific recipies for common tasks.

## Applying the DE pipeline to a matrix of counts

Say, you have a matrix with counts that you have gotten from somewhere,
e.g. GEO, and would like to run the pipeline on that matrix. 

**Solution:** create an input that looks like gene counts from STAR.

Here is how you do it:

 1. In your report directory, create a new directory `input/star`. Why
    star? We will pretend that the counts come from star step of the
    mapping pipeline.

 2. For each sample in your matrix, create a file in `input/star`. Note the
    following:

    - do not use dots in the sample name
    - the file name should be `star.{sampleid}.tsv`
    - the file should contain two columns: gene identifier (preferably ENSEMBL) and integer counts, tab separated
    - prepend the data with four rows (first four rows of a STAR file
      contain summary information and are ignored by the pipeline)

 3. Create the covariates file with sea-snap:

       `sea-snap covariate_file star tsv`

Here is how the first two steps can be done in R, assuming that `mtx` is
the count matrix and that its row names are ENSEMBL identifiers:

```
## get rid of the dots in sample IDs
colnames(mtx) <- gsub("\\.", "_", colnames(mtx))

dir.create(file.path("input", "star"), recursive=TRUE)
mtx2 <- rbind(matrix(0, nrow=4, ncol=ncol(mtx)), mtx)

for(i in 1:ncol(mtx)) {
  f.out <- file.path("input", "star", paste0("star.", colnames(mtx)[i], ".tsv")) 
  write.table(cbind(rownames(mtx2), mtx2[,i]), 
    col.names=FALSE,
    row.names=FALSE,
    quote=FALSE,
    sep="\t",
    file=f.out)
}
```

