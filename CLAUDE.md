# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

CSI-GEP (Consensus and Scalable Interpretations of Gene Expression Programs) finds
gene expression programs (GEPs) in single-cell data. It runs GPU NMF (cNMF-style
consensus, using `torchnmf`) on two random halves of the cells and keeps the programs
that replicate across both halves. The rank `k` is picked automatically.

This checkout is a fork (`miltondp/CSI-GEP`, upstream `geeleherlab/CSI-GEP`) used as a
submodule of the `ai-comp-exp` superproject. Work on the `comp-exp` branch and never
commit to `main` (see the superproject's CLAUDE.md).

There is no package, build, test suite, or linter. The code is standalone scripts
driven by an LSF job script. Dependencies are pinned in `CSI-GEP/requirements_py.txt`
(Python, e.g. torch 1.13.1, scanpy 1.9.3, numpy 1.21) and `CSI-GEP/requirements_R.txt`.
Two Docker images are also published, `ghcr.io/geeleherlab/csi-gep_py` and
`csi-gep_r`, with the scripts under `/app/`.

## Running the pipeline

All paths are relative to the repo root; the scripts are run from there.

```bash
bsub < CSI-GEP/csigep_submit.bsub   # bare-metal (LSF)
bsub < Docker/csigep_submit.bsub    # same pipeline via singularity + the Docker images
```

The job script runs four stages; stages 2 and 3 fan out one LSF job per `k` in `CSI-GEP/K.txt`:

1. `python CSI-GEP/prepare.py -c counts.h5ad --output-dir OUT --numgenes 2000`
   splits the cells 50/50 (`random.seed(1234)`) into `split0`/`split1`. For each
   split it writes the TPM, TPM stats, and variance-scaled counts of the top
   overdispersed genes under `OUT/cNMF_split{0,1}/cnmf_tmp/`. Input must be a sparse
   h5ad.
2. `python CSI-GEP/cNMF.py --output-dir OUT --components K --n-iter 100` (GPU queue;
   falls back to CPU). For each split it runs `n_iter` NMF replicates and then the
   consensus step twice, with density thresholds `2.0` and `0.1`. It writes
   `cNMF_splitX.gene_spectra_score.k_K.dt_0_1.txt`, which is what the later stages
   read.
3. `Rscript --vanilla CSI-GEP/Jaccard.R K OUT/` compares every split0 GEP with every
   split1 GEP. For each Jaccard length (JL = top-N genes, 10..100 by 10), it takes each
   GEP's top-N genes, runs a Jaccard test (`jaccard.test.mca`), and Bonferroni-corrects
   the p-values. Significant pairs form a graph, and its edge-betweenness communities
   are counted. Output: `OUT/results/jaccardtest_results_K{K}.RDS`, plus a per-JL
   file for each `_JL{n}`.
4. `Rscript --vanilla CSI-GEP/GEP_analysis.R OUT/ DATASET GEP_DIR k1,k2,... RESCUE`
   builds community-count vs. rank curves, one per JL. For each curve it fits a loess
   and finds the elbow (`elbow` package), then picks the JL whose post-elbow slope is
   steepest, which gives the rank. If `RESCUE` is 1, it retests the unmatched GEPs at
   the next larger JL. It averages the gene spectra scores within each community to
   get the final GEPs.
   Output in `GEP_DIR`: `*_community_geps.RDS` (genes × GEPs), `*_labels.RDS`,
   `*_gepindex.RDS`, `*_optimized_parameters.RDS`, `*_Curves.pdf`,
   `output_summary.txt`.

Output directory paths are joined by string concatenation in R, so `OUT/` must end
with `/`.

### Gotchas

- `GEP_analysis.R` parses argument 5 as `as.numeric(args[5])`, so it must be `0` or
  `1`. `Docker/csigep_submit.bsub` passes `$rescue` (0/1). `CSI-GEP/csigep_submit.bsub`
  passes `FALSE`, which becomes `NA`. The README changelog describes an
  `as.logical(as.numeric(...))` fix, but that fix is not in this file.
- `prepare.py` and `cNMF.py` duplicate helper functions (`save_df_to_npz`,
  `fast_euclidean`, etc.). Change both copies.
- `cNMF.py` sets `NUMBA_CACHE_DIR=/tmp/`.
- `Jaccard.R` hardcodes `registerDoParallel(16)`. The job script requests 8 cores
  (`-n 8`).

## Other contents

- `simulations.py` generates a Splatter-like simulated dataset (250k cells × 25k genes,
  40 groups plus activity programs) into `simulated_data/`, which must not already exist.
  The dataset itself is hosted on OSF (https://osf.io/tknm2/) and is not in git.
- `CSI-GEP/output_for_benchmarking/` holds committed reference results for that dataset.
  The expected result is 41 GEPs at k = 47, JL = 20, with rescue off.
- `benchmarking/` compares CSI-GEP against SignatureAnalyzer-GPU (`sa_gpu/`), scVI
  (`scvi/`), and iNMF (`inmf/`, R/LIGER). Each has its own requirements file.
  `summary/` has the R scripts that assign cells to GEPs and compute the comparison
  metrics (contingency heatmaps, silhouette, MI, ARI). The step-by-step walkthrough is
  in `benchmarking/Readme.md`. Some steps need manual edits, such as
  commenting/uncommenting lines in `cmd_scvi.py` or setting k in `inmf_run.R`.
- HPC specifics (LSF queues `standard`/`gpu`/`large_mem`, the St. Jude singularity
  registry) are hardcoded in the `.bsub` files. Adapt them for other clusters.
