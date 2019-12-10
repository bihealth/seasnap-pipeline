[Back](../README.html) to main doc.

---

path patterns
-------------

---

A path pattern is a file path with wildcards in braces "{}".
This is the way `Snakemake` handles file paths.

In the pipeline path patterns have to be defined in the config file (`mapping_config.yaml` or `DE_config.yaml`).
Three patterns have to be defined:

- `in_path_pattern` -- for input files
- `out_path_pattern` -- for intermediate and output files
- `log_path_pattern` -- for log files (usually similar to out_path_pattern)

The following wildcards are available in path patterns:

- mapping pipeline
    - `in_path_pattern`: {sample}, {mate}, {batch}, {flowcell}, {lane}, {library}
        - only {sample} is required
    - `out/log_path_pattern`: {step}, {extension}, {sample}, {mate}, {batch}, {flowcell}, {lane}, {library}
        - {step}, {extension} and {sample} are required
- DE pipeline
    - `in_path_pattern`: {step}, {extension}, {sample}, {mate}, {batch}, {flowcell}, {lane}, {library}
        - only {step}, {extension} and {sample} are required
    - `out/log_path_pattern`: {step}, {extension}, {sample}, {mate}, {batch}, {flowcell}, {lane}, {library}, {contrast}
        - {step}, {extension} and {contrast} are required

**Note: {step} and {extension} are not read from the input files.
They are only filled directly by pipeline rules using methods of the path handler.**

`{step}` is the name of a rule and `{extension}` is the file extension of produced output files.
Both are required to save the output of different pipeline steps in unique locations.

**Note: {extension} can contain (in contrast to other wildcards) also `'.'`, e.g. "sorted.bam" would be a valid extension.**

All other wildcards should be used (if they are used) in all path patterns, because if the contained information is not used in different outputs there is currently no point in parsing it, and if it is not parsed, it is not available for the output.

For setting the `in_path_pattern` see also the section [`prepare input`](prepare_input.html)

---

[Back](../README.html) to main doc.