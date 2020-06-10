# Documentation of the various R snippets for the DE pipeline

To actually display any results in the R report, you need to include *R
snippets* in the `report: report_snippets:` section of the
`DE_config.yaml` file.

For some of the snippets, parameters can be specified in the `report:
snippet_parameters` section of the `DE_config.yaml` file. Default snippet
parameters are defined in the `defaults/DE_config_defaults.yaml` file.

Below are brief descriptions of R snippets for the DE pipeline and their
respective parameters.
To see all available snippets, take a look at the contents of the `report`
directory in the main sea-snap directory.

## Contrast snippets

These snippets can be included for each contrast.

### Results of tmod analysis

**Snippet:** `contrast/tmod_contrast.R`

**Description:**

**Parameters:**

|Parameter      |Description                                              |
|----------     |-------------                                            |
|`res_auc_thr`  |AUC threshold for reporting significant enrichments |
|`res_pval_thr` |q-value threshold for reporting significant enrichments |

## Functional annotation

The snippets in this section (directory `report/Functional`) summarize the
results of functional analyses performed on each contrast, or perform a
comparison between contrasts.

### Summary of tmod enrichment analysis

**Snippet:** `Functional/tmod.Rmd`

**Description:** for each of the pre-defined databases (rules `tmod_dbs`
and `tmod`), summarize and compare results for all defined contrasts. The
snippet produces a summary table showing how many significant enrichments
were found for a database in each contrast, a panel plot which compares
different contrasts for selected modules, and a selection of evidence
plots, which show details of enrichments (for example, how significantly
regulated genes are driving an enrichment).

**Parameters:**

The `fig_...` parameters describe which gene sets are selected to be shown
on the panel plot figure summarizing the result. Since the figure compares
several different contrasts, the thresholds refer to the minimum or maximum
value found over all contrasts. For example, to satisfy the q-value
requirement, a gene set must achieve a q-value below the specified
threshold for at minimum one contrast.

|Parameter      |Description                                              |
|----------     |-------------                                            |
|`fig_qval_max` |q-value threshold over all contrasts                     |
|`fig_auc_max` |AUC threshold over all contrasts                     |
|`fig_n_min` |Minimum number of gene sets to show on the figure. If fewer gene sets satisfy the thresholds, the AUC threshold will be relaxed to achieve the `fig_n_min` value.|
|`fig_n_min` |Maximum number of gene sets to show on the figure. If more gene sets satisfy the thresholds, the AUC threshold will be constrained to achieve the `fig_n_max` value.|
|`n_evid` |Number of gene sets for the evidence plots. The gene sets are first selected based on the q-value threshold specified by `fig_qval_max`, and then the top `n_evid` gene sets (by AUC) are selected for presentation.  |


**Notes:** The code of this snippet is quite convoluted. The reason for
that is that first, it needs to summarize information from all contrasts;
next, that for each database and contrasts multiple results can be present,
depending on the `tmod: sort_by` parameter defined for the `tmod` rule.

### Results of the cluster profiling analysis

**Snippet:** `contrast/cluster_profiler_summary.R`

**Description:** Snippet summarises the results of running cluster profiler
on the different contrasts. For each set of gene sets ("database") it
collects the results of the different contrasts and shows them next to each
other on a plot.

