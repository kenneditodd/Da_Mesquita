# Cell Ranger Scripts
- **01_download_ref.sh** downloads the mouse reference dataset required for Cell Ranger
  - Links may be broken when updated. See https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest for the latest version.
- **02_count.sh** is a SLURM job that runs `cellranger count`
- **03_sample_loop.sh** submits 02_count.sh to the SLURM job scheduler with sample name passed to 02_count.sh as an argument. This allows for each sample to have it's own SLURM job submission.
- **04_get_web_summary.sh** will go into each output folder from `cellranger count` and copy the html reports to a single folder for easier access to view.
- **05_get_metrics.R** will merge all metrics_summary.csv reports into a single table.
