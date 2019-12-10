[Back](../README.html) to main doc.

---

Run mapping
-----------

---

<p align="center">
  <img width="400" src="pictures/run_pipeline/run_mapping.svg" />
</p>


---

### 1) edit the config file

After creating a working directory with the file `mapping_config.yaml` you have to edit that file to set your configuration.

SeA-SnaP parses an internal config file with default values from the git repository first and the files you edit in your working directory will overwrite the default values.
In this way the config can be kept simple and you do not need to specify all values.

Many configuration options are available described in a [`separate section`](config_options.html).
(Hint: it might also be helpful to look into the default configuration file in the [git repository](../defaults/mapping_config_defaults.yaml)).

Here, we describe the minimal settings required:

**required settings**

In the config file which was copied to your working directory values that need to be filled are marked with `### FILL IN ###`.
For other sections of the config file `__options__` gives hints to which further adjustments can be made.
But you can also completely delete these other sections, if not required.

At least you need to specify:

- An organism after `organism_defaults:`. This will load the default settings for this organism, which you can overwrite by filling `organism:`. Note: some paths to static GTF, FA, etc. files on the CUBI cluster are set by default. Quite likely you will need to overwrite them. If there are no defaults for your organism, you can also leave `organism_defaults:` blank and only fill `organism:`.
- The folder structure of your input (fastq) files.

Folder structures for input and output files in SeA-SnaP are specified by *path patterns with wildcards*, e.g.:

```
input_files/{sample}/some/path/{sample}.{mate}
```

See [`path patterns`](path_patterns.html) for a detailed description.

See the section [`Prepare input`](prepare_input.html) on how to provide fastq input to the mapping pipeline and how to import from SODAR.

---

### 2) create a sample info file

You will need a file providing sample information to the pipeline, including whether a stranded or unstranded protocol was used for library preparation and whether reads are paired or not.
SeA-SnaP provides a helper function to automatically create a list of samples from your *input folder* containing the fastq files:

```
./sea-snap sample_info
```

This will parse your input folder structure as defined in the config file `mapping_config.yaml`, extract relevant information and write it into a YAML file called `sample_info.yaml`.
Such a file is required by the *mapping pipeline*.

The `sea-snap sample_info` helper also provides options to convert between `tsv` (or `csv`) and `yaml`.
Therefore, you can also manually create the file as a table and then convert it into the required YAML file.
It needs columns `paired_end_extensions` (e.g. `['R1','R2']` or `['']`), `read_extension` (e.g. `.fastq.gz`) and `stranded` (e.g. `unstranded` or `forward`) for each sample ID.

There is also a separate option to import from SODAR, which will additionally take an ISA-tab file as input and merge results from parsing the input folder structure for fastq files and from parsing the ISA-tab file.

*Note: when using auto-generation, check the file after it was created and edit if necessary!*

---

### 3) execute the pipeline

To run the mapping pipeline locally you can use:

```
./sea-snap mapping l
```

This is a wrapper for Snakemake.
You can add any [`snakemake options`](https://snakemake.readthedocs.io/en/stable/executable.html#all-options) after `l` and they will be passed to snakemake.

E.g. run

```
./sea-snap mapping l -nr
```

for a dry run to see which steps will be executed and why.
Other helpful options might be `-R`, `-U` and `-O` to run only specific parts of the pipeline.

---

To run the pipeline on the cluster, type

```
./sea-snap mapping c
```

instead of `l`.

Also note the file `cluster_config.json`, which was added to the working directory, when setting it up with the helper function.
It will be used if SeA-SnaP runs with the cluster option.
In the file, resources for the different pipeline steps are defined.

In addition, the file contains an object `__set_run_command__` where default [`snakemake options`](https://snakemake.readthedocs.io/en/stable/executable.html#CLUSTER) are defined for running the pipeline on the cluster.
You can adapt these to your needs.
E.g. per default `--drmaa` is used to submit jobs to other nodes.

When a job was submitted to the cluster you can run

```
tail -f pipeline_log.out
```

to follow the progress (press `ctrl-c` to exit).

A list of all options for the `./sea-snap` command is given in section [`SeA-SnaP options`](../README.html#sea-snap-options).

---

[Back](../README.html#running-the-pipeline) to main doc.