Here is a HOWTO on running ATAC-Seq

Necessary conda packages to install:

 - macs2

Necessary R packages to install:

 - DiffBind


First, when when running the mapping pipeline (`sea-snap mapping ...`), you
need to configure a few things in the `mapping_config.yaml` file such that
both the BigWig (BW) files are produced and the macs2 program is run.

 1. In the `param:` section, add a section called ATAC_Seq:

     ```
     ATAC_Seq:
      - macs2
     ```

 2. In the `param:QC_results:` section, add `- bw_from_bed`.

Second, you need to configure the `DE_config.yaml` file for ATAC-Seq. The
outline of the procedure is as follows:

 * Use DiffBind to integrate the peaks and generate the peak counts. This
   takes a lot of time.
 * Convert the DiffBind object to raw counts, generate fake PeakIDs.
 * Import the raw counts to DESeq2. Thus, we can proceed with the
   contrasts.
   
