# Cherubism single-cell RNA-seq analysis

This repository contains the version-controlled R scripts used for quality
control, cell-state reassessment, exploratory trajectory analysis, myeloid
differential-expression and enrichment analyses, and candidate cell-cell
communication analysis.

## Study scope

- The computational workflow compares one cherubism (CB) sample with one
  non-cherubism control (CTRL) sample.
- The sampled material was cancellous jawbone. The tissue portion allocated to
  each experimental modality measured `1 x 1 cm`.
- Cells are not treated as independent biological replicates. Between-condition
  results are descriptive and hypothesis-generating.
- Participant-level clinical details are not stored in this public code
  repository because they are not required to execute the analysis.

## Local inputs

Patient-derived data are deliberately excluded from version control. Before
running the workflow, place the following local files in `data/`:

```text
data/Several original1.Rdata
data/Sce_all final with anno.Rdata
```

`Several original1.Rdata` must contain Seurat objects named `CB` and `CTRL`.
The second file contains the legacy annotation used only as an input to the
annotation-reassessment stage. Both files are ignored by Git.

The optional 10X export utility also expects local genome-reference files under
`resources/`; their paths can be overridden with the environment variables
documented in that script.

## Repository structure

```text
R/          numbered analysis scripts
config/     non-identifying analysis configuration and colour palettes
scripts/    stand-alone utility scripts
data/       local patient-derived inputs; ignored by Git
resources/  local reference resources; large files ignored by Git
```

Analysis products are generated locally under `objects/`, `results/`,
`figures/`, `logs/`, and `qa/`. These directories are excluded from the public
repository.

## Software environment

The workflow was developed with R 4.4.x. Package versions are recorded in
`renv.lock`. From the repository root, restore the environment with:

```r
source("00_install_dependencies.R")
```

## Running the workflow

Run commands from the repository root:

```r
source("run_all.R")
```

The main runner rebuilds the quality-controlled base objects. Additional
numbered scripts in `R/` reproduce the annotation, pseudotime, enrichment,
CellChat, and NicheNet analyses. Scripts use project-relative paths and write
only to ignored local output directories.

## Important limitations

The source Seurat objects contain filtered cell matrices rather than unfiltered
droplet matrices. DecontX therefore runs without empty droplets as an external
background. The workflow records this limitation and retains both the original
and corrected count assays for audit.

This repository does not include raw scRNA-seq data, WES files, clinical images,
or directly identifiable participant information.
