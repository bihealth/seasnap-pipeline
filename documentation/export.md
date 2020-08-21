[Back](../README.md) to main doc.

---

Exporting to a different folder structure
-----------------------------------------

---

SeA-SnaP produces many files as by-products of the analyses.
When passing on the results of an analysis only a part of these files may be required.
To this end there is an `export` rule, that can be run after the pipeline has completed.
E.g. run:

```
./sea-snap mapping l export
```

which will run the export- instead of the all rule.

There is a section `export:` in the config file that defines how files are copied into a new folder structure.

Alternatively, if `export: blueprint:` is not None, a file will be written in which source and destination paths are inserted into a command template (`export: blueprint: command:`) and such a block of commands is included for each source-dest pair.

**Per default this is configured to export files for SODAR**

The blueprint can be used to upload the files into `SODAR's landing zone`.
In the file `__SODAR__` has to be replaced with the project uuid to upload the results to.

To automatically replace `__SODAR__` and dispatch the upload to different threads, you can use the [CUBI SAK](https://cubi-gitlab.bihealth.org/CUBI/Pipelines/cubi-swiss-army-knife) tool with the command:

```
cubi-sak sea-snap itransfer-results
```

### export configuration

The config file contains a section `export:` that looks like this:

```
export:
  blueprint:
    file: SODAR_export_blueprint.txt
    command: |
      imkdir -p $(dirname {dest} )
      irsync -a -K {src} i:{dest}
  path_pattern:
    - __SODAR__/{sample}/{GENOME}/%Y_%m_%d/{files:ngs_mapping:out}/{step}/out/{step}.{sample}.{extension}
    - __SODAR__/{sample}/{GENOME}/%Y_%m_%d/{files:ngs_mapping:rep}/{step}/report.zip
    - __SODAR__/{sample}/{GENOME}/%Y_%m_%d/{files:gene_expression_quantification}/{step}/out/{step}.{sample}.{extension}
  ngs_mapping_out:
    - files: {step: star, extension: bam}
  ngs_mapping_rep:
    - dir: {step: star, log: true}
      compress: zip
  gene_expression_quantification:
    - files: {step: star, extension: gene_counts.tab}
    - files: {step: salmon, extension: sf}
```

SeA-SnaP will go through the list of path patterns under `path_patterns:` and try to create the files, while compiling the wildcard values from the config ({GENOME}), replacing [time formatting](https://docs.python.org/3/library/time.html#time.strftime) ("%Y_%m_%d") and filling additional wildcards from other entries under `export:`.

A special wildcard `{files:<A>[:<B>]}` is searched. It will be replaced by `A` and an entry `export: A[_B]:` is looked up. This entry contains instructions for the path handler to construct file paths based on wildcard values. E.g. if `files: {step:star, extension:bam}` is included, the path to the *bam* output files of *STAR* is constructed and this file is copied to the new location, while the {step} and {extension} wildcards of the destination file path under `path_pattern:` are also filled in respectively.

If the wildcard {sample} or {contrast} is in a file path under `path_pattern:` several file paths are constructed, expanding over samples or contrasts, respectively. They are copied to their destination with the {sample} or {contrast} wildcard filled correspondingly.

Whole folders can also be copied to new locations using `dir: {step:<step>}` instead of `files: {step:<step>, extension:<extension>}`.

Instead of only copying them, files and folders can also be compressed by adding a key-value pair `compression: <type>` to a dict.
Supported compression types are at the moment **zip** and **tar**.

---

[Back](../README.md) to main doc.