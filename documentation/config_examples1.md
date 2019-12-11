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
    contrast:
      - Init_code.Rmd
      - MA_plot.Rmd
      - Result_table.Rmd
```

Will first insert the snippets (from the `report/` directory) in this order:

- Covariate_table.Rmd
- PCA_plot.Rmd
- SampleSimilarity_plot.Rmd

Then, the pipeline will create a section "contrast" and add sub-sections for all (special keyword `__list__`) contrasts (special keyword `__contrast__`) defined in the config file.
For each contrast, report snippets will be added in this order defined under `defaults: contrast:`:

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
    contrast:
      - Init_code.Rmd
      - MA_plot.Rmd
      - Result_table.Rmd
```

Here, for the contrast "B vs C" the default snippet list is overwritten and only the snippet "Result_table.Rmd" is used.

---

[Back](../run_DE.md) to run DE.