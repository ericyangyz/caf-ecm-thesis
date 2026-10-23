# CAF matrix organization and tumor cell migration

Code for my undergraduate honors thesis in the Yeh Lab at UNC.

Pancreatic tumor cell lines were imaged for 24 h on cell-free matrix made by CAFs from
patient donors (proCAF or restCAF by DeCAF). The scripts here measure how aligned the
fibronectin fibers are, whether cells move along the fiber right under them, and whether
matrix organization relates to the matrix proteome.

Raw images, tracking output and proteomics data are not included.

## Folders

| folder | what it does | thesis figure |
|---|---|---|
| `decaf/` | ECM genes among the DeCAF gene pairs | Fig 1 |
| `ecm_orientation/` | fiber orientation in the 20x fibronectin images (OrientationPy vs Mathematica), matrix thickness and nuclei QC | Figs 2, 3 |
| `fiber_maps/` | per-pixel fiber orientation for every frame of the movies | input for Fig 4 |
| `guidance/` | angle between each cell step and the fiber under it, guidance index G | Fig 4 |
| `proteomics/` | matrisome proteins vs matrix organization (limpa / limma) | Fig 6 |
| `figures/` | panels drawn at final size | Figs 1D, 4E-H, 5, 6B, S1-S3, Table S2 |

## Running

Everything is run from the repo root, with `data/` and `results/` folders next to the code.
The `.sbatch` files are for SLURM on UNC Longleaf.

Python: `conda env create -f environment.yml`.
R 4.5 with tidyverse, here, readxl, patchwork, ggrepel, ggpubr, pheatmap, RColorBrewer,
limma, limpa, AnnotationDbi and org.Hs.eg.db.

Migration part, in order:

1. `fiber_maps/split_channels.py`, then `fiber_maps/run_wholefield_pixel.sbatch <sample> <field>` for each field
2. `guidance/01_reduce_fiber_field.py` for each field (`sbatch --array=1-N%6 guidance/run_reduce.sbatch`)
3. `guidance/02_spots_to_trajectories.py`, `guidance/03_guidance.py`
4. `guidance/extract_steps.py`, `guidance/alignment_scales.py` (local fiber order R_local)
5. `guidance/04_delta_distribution.py`, `guidance/05_fig_delta.R`, `figures/fig4_panels.R`

Proteomics: `proteomics/01_design.R` to `04_figures.R`. `ecm_orientation/uma_result_analysis.Rmd`
has to be run first, since its per-donor fiber spread is the predictor.

## Notes

- 0.64 um per pixel, 30 min per frame, 48 frames.
- Image y runs down, fiber angles are measured with y up, so step angles use `atan2(-dy, dx)`.
- Group comparisons use field or donor means, not single cells (Fig S3 shows the per-track version).
- Cell tracks are from TrackMate (Kristoff Chen). Stiffness and marker data are from Ian McCabe.
- `fiber_maps/alignment_analysis_pixel.py` is a modified copy of `alignment_analysis.py` from
  [UMA-tools](https://github.com/alexdolskii/UMA-tools).
