# Stratification Pipeline

This is the newer stratification pipeline with updated variant definitions.
Other than variant definitions the pipelines are pretty much the same.

## Running the pipeline

Symlink the `Assembly/` directory under `Assembly`.
This pipeline is primarily for simulated libraries.
Be sure to name it appropriately.
Add the symlinked directory name to the config file under `samples`.
Then just run snakemake as below.

```
# -s supplies the snakefile
# --configfile supplies the config file
# -j supplied threads
snakemake -s run_binning.smk --configfile run_binning.config -j 60 2>&1 | tee snakemake.err.txt
```

## run_binning.config

- samples
    - samples specifies which samples to evaluate
    - just add the name of the dir that's been added to `Assembly` before running the pipeline
