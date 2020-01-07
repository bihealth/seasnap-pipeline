[Back](../run_DE.md) to run DE.

---

Config examples
---------------

---

### Set contrasts

In the config file under `contrasts: contrast_list:` you can define a list of contrasts that should be computed.

E.g.

```
contrast_list:
  - title: "A vs C"
    ratio:
      column: condition
      numerator: A
      denominator: C
  - title: "B vs C"
    ratio:
      column: condition
      numerator: B
      denominator: C
    ranking_order: "-x"
    results_parameters:
      altHypothesis: less
```

defines two contrasts "A vs C" and "B vs C".
"A", "B" and "C" are levels of a column "condition" in the covariate file.
In the second contrast, also defaults for `ranking_order` and `results_parameters` are overwritten.
The defaults are defined in the config file under `contrasts: defaults:`.

**conditional functional annotation analysis**

Analysis of e.g. GO and KEGG annotation enrichment can be optional for certain contrasts.
For example by adding a key-value pair `goseq: true` to a contrast specification, **goseq** is run for that contrast for GO and KEGG enrichment analysis.

At the level of the snakefile, the argument `if_set = dict(goseq=True)` of the `expand_path()` path handler function is used when collecting output files for the `all` rule.
While expanding over the contrasts, only output files will be requested, if a `dict(goseq=True)` is in the definition of the respective contrast.
This mechanism can also be used, when adding more analysis steps to the pipeline that are optional for contrasts.

If results of some analysis steps are not computed for all contrasts, also the report must adapt by not including sections for those steps for a given contrast.
See [`adding new snippets`](adding_rmd_snippets.md) on how to make snippet inclusion conditional.

**included functional annotation methods so far:**

- `goseq` for GO and KEGG over-representation analysis
- `cluster_profiler` for MSigDb over-representation or gene set enrichment analysis

### Configure report

In the config file under `report:` you can define a 'building plan' how to assemble the report from Rmd snippets.

E.g.

```
report:
  report_snippets:
    - Covariate_table.Rmd
    - PCA_plot.Rmd
    - SampleSimilarity_plot.Rmd
    - contrast:
        - __list__: __contrasts__
  defaults:
    contrast_list:
      - Init_code.Rmd
      - MA_plot.Rmd
      - Result_table.Rmd
```

Will first insert the snippets (from the `report/` directory) in this order:

- Covariate_table.Rmd
- PCA_plot.Rmd
- SampleSimilarity_plot.Rmd

Then, the pipeline will create a section "contrast" and add sub-sections for all (special keyword `__list__`) contrasts (special keyword `__contrast__`) defined in the config file.
For each contrast, report snippets will be added in this order defined under `defaults: contrast_list:`:

- Init_code.Rmd
- MA_plot.Rmd
- Result_table.Rmd

Instead of using `__contrast__`, contrasts can also be specified as a list, and defaults can be overwritten:

E.g.:

```
report:
  report_snippets:
    - Covariate_table.Rmd
    - PCA_plot.Rmd
    - SampleSimilarity_plot.Rmd
    - contrast:
      - __list__:
        - "A vs C"
        - "B vs C":
          - Result_table.Rmd
  defaults:
    contrast_list:
      - Init_code.Rmd
      - MA_plot.Rmd
      - Result_table.Rmd
```

Here, for the contrast "B vs C" the default snippet list is overwritten and only the snippet "Result_table.Rmd" is used.

Defaults for list entries (created by `__list__`) have the suffix `_list`. Defaults for *blocks* can also be defined and have no suffix.
This definition would be equivalent to the previous example:

```
report:
  report_snippets:
    - Covariate_table.Rmd
    - PCA_plot.Rmd
    - SampleSimilarity_plot.Rmd
    - contrast: __defaults__
  defaults:
    contrast_list:
      - Init_code.Rmd
      - MA_plot.Rmd
      - Result_table.Rmd
    contrast:
      - __list__:
        - "A vs C"
        - "B vs C":
          - Result_table.Rmd
```

Here, `contrast: __defaults__` would load `contrast:` from the defaults. There, in turn, via `__list__` the defaults from `contrast_list:` would be loaded for the entry "A vs C".

Note: *blocks* do not have to contain a `__list__`, they may simply be groups of snippets that are inserted under a common subsection.

Tip: when configuring the report, it can be helpful to copy, paste and edit the report configuration from the [default yaml file](../defaults/DE_config_defaults.yaml).


### Merging reports of several analyses

You might want to perform different analyses of the same data, using e.g. different design formulas or filtering.
To this end, it is possible to create a merged report of several analyses.

Assume for example that results were stored using these two `out_path_patterns`:

```
previous_analysis/{contrast}/{step}/out/{step}.{contrast}.{extension}
current_analysis/{contrast}/{step}/out/{step}.{contrast}.{extension}
```

If the second path pattern is the current path pattern, you may add the results of the first analysis to its report with the line:

```
report:
  merge: previous_analysis/{contrast}/{step}/out/{step}.{contrast}.{extension}
  report_snippets: ...
  defaults: ...
```

Then all snippets will be inserted twice, alternating between results of the two analyses.
Tags will be added to the sections: "-- analysis0" for the current analysis and "-- analysis1" for the added analysis.
Suffixes "_analysis\<x\>" will also be added to relevant variables in the report, to load and use the corresponding config files and lists of generated results files.
Therefore, if you wish to replace the tag in the section names, do a replace with e.g. "-- analysis0" by "-- \<my title\>" including the "--".

The value of merge can also be a *list* or a *dict*.
Especially if you don't need all snippets inserted twice, you can select which ones to add:

```
report:
  merge:
    pre_analysis: previous_analysis/{contrast}/{step}/out/{step}.{contrast}.{extension}
  report_snippets:
    - Covariate_table.Rmd
    - PCA_plot.Rmd
    - SampleSimilarity_plot.Rmd
    - contrast: __defaults__
    - pre_analysis:
        contrast: __defaults__
```

Here, only the section contrast is added twice to the report, once with the tag *analysis0* for the current analysis and once with the tag *pre_analysis* for the previous analysis.
In general, the value of merge can be a dictionary with multiple keys and either strings or lists as values.

**Note: using tags like `pre_analysis:` inside of blocks like `contrast:` is not recommended, since some global variables might be defined there, that are used across different snippets. The values of such variables might differ between analyses, though.**

*Note: combining path patterns with differing wildcards in one report was not tested yet and might result in error messages.*

---

[Back](../run_DE.md) to run DE.