[Back](../README.html) to main doc.

---

Prepare input
-------------

---

### fastq files -- folder structure

SeA-SnaP can be adapted quite flexibly to read fastq files from different folder structures.
To this end the folder structure needs to be expressed as a valid input path pattern (see [path patterns](path_patterns.md)).

E.g.

```
../input/{sample}/{flowcell}.{lane}.{mate}
```

would be a valid option.

**Note: Wildcard values should not contain the characters `'/'` or `'.'`, since these are used to separate wildcards!**

Inside the `input_path_pattern` you can use:

- normal text
- wildcards in braces, e.g. `{sample}`, also with an optional regular expression after a comma, e.g. `{sample,[^_]+}`
- `'*'` and `'**'` can be used for any variable parts that are not used in the pipeline (unlike wildcard values). If there are multiple matches, the files will be pooled for a respective sample, e.g. running STAR on all those files together

`'*'` matches everything, except `'/'` and `'.'`.

The file extension, e.g. `.fastq.gz` is left out from the path pattern, because different files for different samples may be compressed differently and Sea-Snap determines that automatically.

---

### SODAR import

To import data from SODAR, use iRods as explained elsewhere to download the **fastq files**.
The `input_path_pattern` can then be adapted to match the folder structure.

E.g. for files like

```
inputs/YRJZTT-T1-RNA1-mRNA_seq1/YRJZTT-T1-RNA1-mRNA_seq1_R1_002.fastq.gz
```

put

```
inputs/{sample}/*_{mate,R1|R2}_*
```

to parse "YRJZTT-T1-RNA1-mRNA_seq1" as sample ID and "R1" as mate name.

Further, you can download an **ISA-tab file** from SODAR with meta information about the samples.

the `sample_info` wrapper has an option to parse both the folder structure and the ISA-tab file and merge the results into a sample_info file.
Run:

```
./sea-snap sample_info --from sodar --input path/to/isa-tab
```

*after* editing the path pattern in the config file.

---

[Back](../README.html) to main doc.