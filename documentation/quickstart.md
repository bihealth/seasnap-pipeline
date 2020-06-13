# Before you start




# Mapping pipeline

 1. Create a directory for running the pipeline (you can use `sea-snap.py working_directory DIR` for that).

 2. Make sure that the directory contains the following objects:

   * symbolic link to the sea-snap.py executable (call it `sea-snap` and
     not `sea-snap.py` because the former looks nicer)
   * symbolic link to the input directory with the FASTQ files
   * do yourself a favor and make a symbolic link to the sea-snap directory
     (call it `sea-snap-dir`)
   * the file `mapping_config.yaml` (copy it from the sea-snap directory if missing)
   * the file `cluster_config.json` (copy it from the sea-snap directory if missing)

 3. Edit the file `mapping_config.yaml`. In most cases, the comments are
    self-explanatory. Here a few hints:

    * in general, take a look at the default configuration file in
      `sea-snap-dir/defaults/mapping_config_defaults.yaml`
    * organism: take a look at `sea-snap-dir/defaults` to see what organism
      files are available.
    * Wildcards: the Python/snakemake wildcards are identifiers in curly
      braces, such as `{this}`.
    * under `pipeline_param`, you should have three wildcard infested paths
      defined. These are somewhat tricky, so here are two examples, each
      with example of the resulting *real* paths:
      
         ```
         ## Example 1.
         ## Note that the '.fastq.gz' extension in input file is defined elsewhere
         ## The files should have extensions '.fastq.gz', '.fastq', '.fq.gz' or '.fq'
         ## If they don't, see below

         ## example input files:
         ## input/S1254/S1254_R1.fastq.gz
         ## input/S1254/S1254_R2.fastq.gz
         in_path_pattern: input/{sample}/{sample}_{mate}

         ## example output file:
         ## notice how the `{mate}` wildcard got replaced by `all_mates`
         ## mapping/star/out/S1254_all_mates/star.S1254.all_mates.bam

         ## matching out/log path patterns:
         out_path_pattern: mapping/{step}/out/{sample}_{mate}/{step}.{sample}.{mate}.{extension}
         log_path_pattern: mapping/{step}/report/{sample}_{mate}/{step}.{sample}.{mate}.{extension}

         ## Note that the `{sample}` wilcard can be eventually 
         ## replaced by something like `all_samples.all_mates`
         ## mapping/multiqc/out/all_samples.all_mates/multiqc.all_samples.all_mates.qc_report.html

         ## Example 2.
         ## example input file:
         ## input/1509/1509_R1_001.fastq.gz
         in_path_pattern: input/{sample}/{sample}_R1_001

         ## matching out/log path patterns (no mates here):
				 out_path_pattern: mapping/{step}/out/{sample}/{step}.{sample}.{extension}
				 log_path_pattern: mapping/{step}/report/{sample}/{step}.{sample}.{extension}
         ```

    * One important point about the above paths: each output and log file name 
      (the last
      element of the path, i.e. for example '{step}.{sample}.{extension}'
      in the example above) *must* include *all* wildcards used *anywhere*
      in the path.
    * also under `pipeline_param`, you should choose which steps (e.g.
      dupradar, star counts or multiqc) are to be run
    * edit other options and parameters to your taste, but basically you
      ready to go.

 4. Generate the file `samples.info` using the command `sea-snap.py sample_info`

    * If you have a custom extension to your read files, then use the
      `--add_ext` option.
    * A common problem is that the `in_path_pattern` defined in the
      `mapping_config.yaml` file does not match the actual file paths. 
    * After successfully running this step, make sure that the
      automatically generated `samples.info` file is correct.

 5. Run the pipeline with `./sea-snap mapping l` on a local machine or
    `./sea-snap mapping c` on the computing cluster. To run on SLURM nodes,
    run `./sea-snap mapping --slurm c`.
 
 6. So, where is the report? Where is the QC? Just take a look at
    `mapping/multiqc/all_samples.all_mates/out/multiqc.all_samples.all_mates.qc_report.html`,
    it probably contains all the information you are looking for.
    Individual results (fastqc, dupradar) are found under `mapping/fastqc`,
    `mapping/dupradar` etc. Just mind that you want the `out` subdirectory
    and not the `report`; the `report` is the pipeline report (ie. the
    error messages, logs etc.).

# DE pipepine

You probably want to run this pipeline from the same directory in which you
have run the mapping pipeline, as this makes it so much easier.

Be warned, DE pipeline is much more tricky to run than the mapping
pipeline. But you already know it.

## Quick start to the DE pipeline

 1. Edit the `DE_config.yaml` file. Make sure that the `in_path_pattern`
    matches files which were produced by the mapping pipeline. Which files
    are these depends on whether you want to take the count data from STAR, 


## Initial steps



## Configuration file `DE_config.yaml`


### The report snippet section

Here is a weird one from the yaml config file:

```
report:
  report_snippets:
    - contrast:
        - __list__: __contrasts__
```

Actually it is not weird, but it is just sea-snap pipeline extension to the
YAML format.

This construction essentially means: for each contrast, run the defaults
for the contrasts which are defined below, like this (also in the report
  section):

```
report:
  report_snippets:
  # ...
  defaults:
    contrast_list:
      - Init_code.Rmd
      - MA_plot.Rmd
      - Result_table.Rmd
      - Goseq_GO_table.Rmd
```

Having thus set up the report, the code snippets will be run once for each contrast.


## Final report

Once the mapping pipeline finishes, it has not generated an HTML report, but an
Rmarkdown file which you can now edit at will. You will find it (depending
on how you defined the input path) under something like
`DE/all/report/out/report.all.Rmd`. 

The full path to the project is hard-encoded in the Rmd file, so if you copy it
from the cluster, even if you copy the whole project directory, it will
most likely not run correctly. On the other hand, you can copy it literally
anywhere on the same filesystem and you will be able to render it
correctly, which is something.

Render the project from command line with

```
Rscript -e 'rmarkdown::render("report.all.Rmd")'
```
