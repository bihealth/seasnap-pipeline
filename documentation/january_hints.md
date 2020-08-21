# January's hints on writing sea-snap rules

You write the new rules in the `*.snake` file (e.g. `DE_pipeline.snake`).

When writing new rule, remember to make sure that the name of the rule is
mentioned in all of the `file_path` calls (sections input, output etc.), as
well as the `pph.log` call at the end of the `run:` section.

## How Snakemake works

Snakemake is based on generating files (not on running tasks). The only
files that snakemake "wants" to produce are these which are the outputs of
rule "all". However, these require all the files mentioned in the "input"
section of the rule "all". To generate these, outputs of all the rules are
matched, thus generating a chain of dependencies.

## Contrasts or "all"


The way most of the file path generation functions (such as `file_path` or
`log`) are written, by default the generated file template will have a
`{contrast}` placeholder, which then will be matched if a file with the
specific `{contrast}` is requested.  However, these functions take the
parameter `contrast=all` which replaces the placeholder `{contrast}` by
`all`.

So, while

```
pph.file_path("whatever", "whatever.rds", contrast = "all")
```

generates something like

```
DE/all/whatever/out/whatever.rds
```

the following statement

```
pph.file_path("whatever", "whatever.rds")
```

generates

```
DE/{contrast}/whatever/out/whatever.rds
```

and will allow to match any of the files in which `{contrast}` stands for a
specific contrast.

Therefore, for contrast-independent rules, define parameter `contrast="all"` in all
the `file_path` calls *as well as* the `pph.log` call and the end of the
`run:` section, for example

```
    script_file = pph.log(log.out, snakemake_format(script),
        step="tmod_dbs", extension="R", contrast="all", **wildcards)
```

If no other rule depends on the new rule, the files produced by that rule
must be added to the required `input` for the rule `all`. For this, it is
necessary to modify the function `get_inputs_all()`. However, this
modification depends whether we want to have one output file for all
contrasts (thus, stored under "DE/all/"), or whether we want to have a
series of output files, one for each contrasts.

If a rule is to be run once for all contrasts, we add the following
statement to the definition of the `get_inputs_all` function (mind the
"all" parameter!).

```
  inputs.append(pph.file_path(step = "annotation", extension = "rds", contrast="all"))
```

If we want a series of output files, we add

```
  inputs.append(pph.expand_path(step = "tmod", extension = "rds", if_set = dict(tmod=True)))
```

The `if_set` parameter and the really weird syntax that follows basically
mean that if the configuration specifies that tmod should be run for all or
for that given contrast, a respective file will be added to inputs, thus
triggering the rule for that contrast and that step.




## Accessing the global config from R

The by far easiest way to do it is to directly read the yaml file which has
been created at the beginning of the snakemake file. This has the great
advantage of avoiding a multitude of placeholders (wildcards), mixing R and
Python code and the need to define them in Python; basically, we just use
a single placeholder and then a native R object "config" to access whatever
information we need.

This file (despite the misleading `step` parameter) is not the result of
the rule `pipeline_report`, it has been created explicitely using "onstart"
block in the snakemake file. Nonetheless, we specify it as an input
requisite:

```
rule whatever:
  input:
    config_file = pph.file_path(step="pipeline_report", extension="yaml", contrast="all")
  ...
```

Then, in the R script, we use `{input.config_file}` to access the global
configuration:

```
   run:
     script=textwrap.dedent(r"""
      conf.f   <- "{input.config_file}"
      config <- yaml::yaml.load_file(conf.f)
       ...
     """)
```

In general, I would avoid *any* wildcards except at the very beginning of
the script, just to separate the mixed R/snakemake/Python code from pure R code:

```
    script = textwrap.dedent(r"""
      {R_SESSION_INFO}

      res.file <- "{output}"
      script.d <- "{SCRIPTDIR}"
      conf.f   <- "{input.config_file}"

      ## no wildcards beyond this point

      config <- yaml::yaml.load_file(conf.f)

      source(file.path(script.d, "tmod_functions.R"))
 
      ## read the necessary files
      dbs <- process_dbs(config$tmod)

      saveRDS(dbs, file=res.file)
    """)



```

## Debugging rules

Since the pipeline is a mixture of Python, snakemake, yaml, R and shell,
debugging can be challenging at times.

Just a few points to remember:

 * for testing, use `./sea-snap DE l` 
 * you can put a `print()` Python statement in the code of a `run:` section
   and you will see the output when you run sea-snap
 * if your rule breaks within the R script, the actual R script (with all
   the wildcards filled in) has been saved in
   DE/{contrast}/{rule}/out/{...}.R, so you can just start an R session and
   source or copy&paste it for debugging.


## Annotation, gene IDs etc.


There is an issue with how the primary IDs are handled.

Basically, my work assumes that `rownames(dds)`, where `dds` is the object
generated by the `DESeq2` rule hold the ENSEMBL identifiers of the
features. This is dangerous, as it is hard to see how that is enforced, as
DESeq2 rule can take a multitude of inputs (including, for example, gene
counts). Respective care should be taken in the downstream code.

