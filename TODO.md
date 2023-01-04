The roadmap ahead:

1) merging Eric's branches and the ATAC_seq branch
2) making sure sea-snap runs smoothly with newer versions of snakemake (it
doesn't currently) --> DONE for mapping pipeline if snakemake version = 7.19.1
3) making the step from drmaa to cluster profiles --> mostly DONE for mapping pipeline 
4) optimizing resource requirements in `mapping_pipeline.snake` (check `mem` vs `mem_per_cpu`) and adding them to `DE_pipeline.snake`
5) allowing sample-specific indices for bwa / salmon / kallisto?
6) adapting the `conda_env.yaml` files
