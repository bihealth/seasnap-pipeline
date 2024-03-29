[Back](../README.md) to main doc.

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

Inside the `in_path_pattern` you can use:

- normal text
- wildcards in braces, e.g. `{sample}`, also with an optional regular expression after a comma, e.g. `{sample,[^_]+}`
- `'*'` and `'**'` can be used for any variable parts that are not used in the pipeline (unlike wildcard values). If there are multiple matches, the files will be pooled for a respective sample, e.g. running STAR on all those files together

`'*'` matches everything, except `'/'` and `'.'`.

The file extension, e.g. `.fastq.gz` is left out from the path pattern, because different files for different samples may be compressed differently and Sea-Snap determines that automatically.

---

### SODAR import

<p align="center">
  <img width="600" src="pictures/import/sodar_import.svg" />
</p>

To import data from SODAR, use iRods to download the **fastq files**.
Further, you can download an **ISA-tab file** from SODAR with meta information about the samples.

The `input_path_pattern` in the config file should then be adapted to match the folder structure.

E.g. for files like

```
inputs/YRJZTT-T1-RNA1-mRNA_seq1/YRJZTT-T1-RNA1-mRNA_seq1_R1_002.fastq.gz
```

put

```
inputs/{sample}/*_{mate,R1|R2}_*
```

to parse "YRJZTT-T1-RNA1-mRNA_seq1" as sample ID and "R1" as mate name.

The `sample_info` wrapper has an option to parse both the folder structure and the ISA-tab file and merge the results into a sample_info file.
*After* editing the path pattern in the config file, run:

```
./sea-snap sample_info --from sodar --input path/to/isa-tab
```

A file called `sample_info.yaml` should be created, which contains information about the samples.
Check the file, edit other relevant sections of the config file and you are ready to run the pipeline.

**Note: The solution of parsing information from ISA-tab files is more or less temporary, since an API for SODAR is about to come at some point. 
Nevertheless, until then the sea-snap helper can be used to extract information from the file.**

Unfortunately, sometimes there are errors in the ISA-tab files, so please make sure to double-check the extracted information.

Further, the recognition of different ISA-tab column names and extraction of values is rudimentary at the moment, but can be extended in a separate configuration file that is located under [`tools/ISAtab_parse_conf.yaml`](../tools/ISAtab_parse_conf.yaml).
It contains entries of this form:

```
<field>:
  columns:
    - <regex>
    - <regex>
    - <regex>
  value:
    <regex>: <replacement string>
    <regex>: <replacement string>
```

For a field name (which will be added to `sample_info.yaml`) two keys are defined:

- `columns` contains regular expressions to recognize the name of the column that contains information for the respective field. Regular expressions in the list will be tried out consecutively until a match was found.
- `value` contains pairs of regular expressions and replacement strings, that will be used with pythons [`re.sub`](https://docs.python.org/3/library/re.html#re.sub) function. The substitution will be applied to all values in that column and thus allows trimming and formatting. Note: values will be converted to lower case before applying the substitution.

Thus, if you download an ISA-tab file, have a look at it and check whether the right regular expressions are included in the config.
Feel free to extend the config file and push it back to gitlab. In this way the configuration may grow and improve.

For an explanation how to export to SODAR see [`SODAR export`](export.md).

**Note: to use sample-specific STAR indices, add entries `star_index: /path/to/star_index/` (careful: must be compatible with STAR version in the pipeline!) or `genome: /path/to/genome.fa` and `gtf: /path/to/gtf` for auto-creation to the `sample_info.yaml`**

---

[Back](../README.md) to main doc.
