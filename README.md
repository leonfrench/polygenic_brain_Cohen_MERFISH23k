# Polygenic gene-set tester for the Cohen et al. RT&T-AMP-MERFISH mouse brain

A Streamlit app for the whole-transcriptome RT&T-AMP-MERFISH data of Cohen et al. (2026, *Cell*): seven coronal mouse brain
hemisections, 276,426 cells and 22,930 measured genes. Enter a gene set; every cell gets an exact set-versus-complement AUROC from its
gene expression ranks, shown on the seven hemisections in two rows, anterior to posterior. The app also tests whether the gene set is
specifically expressed in each **cell type** and each **brain region**, using pseudobulk profiles as in the
[SEA-AD polygenic tester](https://github.com/leonfrench/SEA_AD_multi_polygenic_app). It is a cell-based relative of the
[spatial fetal app](https://github.com/leonfrench/spatial_fetal_app), and works like the apps for the
[Stereo-seq](https://github.com/leonfrench/polygenic_brain_Han_StereoSeq) and
[500-gene MERFISH](https://github.com/leonfrench/polygenic_brain_Yao_MERFISH500) mouse brain atlases. Samples 1–4 are the deposit's
whole-transcriptome (`WT`) set and samples 5–7 its isoform (`ISO`) set, of which only the genes are used; the ~10,000 isoforms the paper
also imaged are left out.

## Run locally

```bash
~/miniconda3/bin/conda create -y -p ./.conda-env -c conda-forge --override-channels python=3.12 pip
./.conda-env/bin/python -m pip install -r requirements-prep.txt pytest
./.conda-env/bin/streamlit run app.py
```

Open http://localhost:8501. `requirements.txt` is enough to run the app on the prepared data in `data/rttamp-genes/`;
`requirements-prep.txt` adds h5py for preprocessing. `RTTAMP_APP_DATA=/path/to/export` selects a different prepared dataset.

## Deploy on Streamlit Community Cloud

- Main file: `app.py`. Python 3.12 is the development version. `requirements.txt` is the dependency file Community Cloud installs;
  do not add an `environment.yml`, which it would install instead.
- The prepared data (488 MB: 403 MB of cells in `data/rttamp-genes/` and 85 MB of bins in `data/rttamp-genes-bins178um/`; every file is
  under GitHub's 100 MB limit, the largest is 27 MB) is committed in plain git, so no Git LFS is needed.
- Memory: Community Cloud allows up to 2.7 GB. On a local server, a browser session that scored every cell and the bins, switched
  through colour modes and levels and then scored a second gene set peaked at a physical footprint of 0.39 GB (about 1 GB resident).
  Label columns load as categories, and the plot geometry is built once and shared by all sessions.
- The first page load sends about 4.6 MB of plot data to the browser.

## Using the app

- **Sections**: all seven samples in two rows (four and three), anterior to posterior and labelled by sample number: 1, 2, 4, 7, 5, 6, 3.
  Samples 1–4 are from one mouse, 5 and 6 from one mouse and 7 from another; the paper does not say whether the mouse of 1–4 is one
  of these. The files give no anterior–posterior coordinate, so preprocessing estimates one (see [Preprocessing](#preprocessing)):
  samples 1 and 2 sit at 4.14 mm (striatum and lateral septum), 4 and
  7 at 7.51 mm, 5 and 6 at 7.69 mm and 3 at 7.84 mm (hippocampus and thalamus). Sections at the same plane keep sample order. For 4 and
  7, which the registration cannot separate, the cells agree: both still have the tail of the striatum, but only 4 has the medial
  habenula, which in the atlas ends 75 µm before the striatal tail, so 7 lies behind 4. Samples 5 and 6 keep the habenula's posterior tip
  but have lost the striatal tail, the reverse of 7, so the sections are not all cut exactly in the atlas plane; 7 stays before them, as
  registered. Each is a hemisection shown at its CCF position (`ccf_z` across, `ccf_y` down), so dorsal is up and the midline (5.7 mm) is
  on the left, as in the paper's Fig. 2D; 2–11% of each section's cells lie just past the midline.
- **Genes**: mouse symbols by default (matched case-insensitively). Human or rhesus macaque symbols are converted to mouse orthologs
  with the bundled HomoloGene snapshot (`data/homologene_to_human.csv.gz`); one-to-many orthologs are all kept and unmapped symbols
  are reported. An optional background is the complete universe, including the targets; ranks are recomputed inside it.
  **Gene matching** lists the targets used and any symbols left out, and its download gives every gene of the scoring universe with
  its role (target or background). The universe is every measured gene, coding or not (see [Data](#data)). The paper left out about
  90 of the most abundant brain genes, so common markers such as `Snap25`, `Syt1`, `Mbp`, `Plp1`, `Ttr`, `Gad1`, `Npy`, `Cck`,
  `Tac1` and `Apoe` are reported as not in the dataset. The default gene set is the medial habenula markers of the
  [MERFISH version](https://github.com/leonfrench/polygenic_brain_Yao_MERFISH500) plus `Gpr151` (`Tac2 Chat Slc5a7 Chrnb3 Pou4f1 Gpr151`;
  the Stereo-seq version's set also has `Tac1`, which this panel lacks). It ranks `Cholinergic_TH` and MH first (AUROC 0.99). The
  atlas MH falls in samples 4 and 7 (284 cells), the two sections registered to the same plane, but only sample 4 has habenula neurons
  there (122 `Cholinergic_TH` cells); sample 7's cells there have no `Tac2`, `Chat`, `Slc5a7`, `Pou4f1` or `Gpr151`. Samples 5 and 6 hold
  the habenula's posterior tip (203 `Cholinergic_TH` cells), although their registered plane lies just past the atlas MH. The paper's
  rostral lateral septum module (Fig. 3F: `Glp1r Rasd1 Trpc4 Zic1 Zic2 Zic3 Zic4 Zic5`) ranks `Striatum_LSr` (0.97) and the LSr (0.98)
  first, in both anterior sections.
- **Colour cells by**: gene-set AUROC, cell type AUROC, region AUROC, then cell type or region (the two label AUROC modes colour each
  cell by its label's specificity AUROC from the tables). The cell type level (class / group / type) and region level (division /
  structure / substructure) set both the colours and the tables. A fourth option, before cell type and region, shows the gene-set AUROC
  of [square bins of cells](#binned-view) instead of the cells.
- **Cell types**: the source's 42 cell types (`cell_label`, the default level): the 41 of the paper's de novo clustering of the WT set
  (Fig. 4A,B), non-neuronal classes and neuronal subtypes such as `Cortex L2/3 IT` or `Striatum_LSr`, and `Striatum_other`, which only the
  ISO set's clustering has (585 cells, 532 of them in sample 7) and which the paper puts with the striatal neurons. Six WT types
  (`Striatum_D1`, `_D2`, `_LSr` and `_OT`, `Hippocampus_IG` and `Cholinergic_STR-HY`) do not occur in the posterior ISO sections. Two
  levels group them as the paper does: 23 cell groups (those of its ligand–receptor analysis, Fig. 4G,H, with 'Neu' spelled out, as in
  `Cortex IT excitatory neurons`) and 13 cell classes (Fig. 4A: `Neurons`, the 11 non-neuronal classes and `Other`, 9,270 cells; the WT
  set's 3,792 have a mixed marker profile, with some `Slc17a7`, `Gad2` and `Flt1`). The tooltip spells out `NFOLs`, `OPCs`, `SMCs` and
  `VLMCs`.
- **Plot**: one panel per section, all sharing the same CCF frame and scale in micrometres, on a black background. Dots are drawn 4%
  wider than the median nearest-cell spacing (~14 µm; `DOT_SCALE` in `app.py`). AUROC ranges differ between sections, so gene-set AUROC
  colours are relative within each section: each is scaled between its own 1st and 99th percentiles, which are shown above its panel,
  and the legend is in arbitrary units. This also absorbs the ISO sections' lower depth (see [Data](#data)). Hover gives the exact AUROC.
  The label AUROC modes (one pooled value per label) share one scale clipped at the 1st and 99th percentiles of the cells shown. The
  clip is `COLOUR_CLIP` in `app.py`. Colours only; the statistics are unchanged. Drag to pan, scroll to zoom, double-click to reset,
  click a legend entry to highlight that label, and save the current view as PNG. Zooming out stops at the fitted view; from there the
  wheel scrolls the page (a pinch or ctrl+wheel never zooms the page). Hover shows a cell's gene-set AUROC, then its cell type, its
  structure (name, then acronym in brackets) with the sample, and its substructure; the class, group or division is added when that
  level is chosen in the sidebar.
- **Tables**: cell type and region expression specificity across all sections (AUROC, p-value and the Benjamini–Hochberg adjusted
  pFDR), with CSV downloads, plus a download of every cell's score once scoring completes. The tables use pseudobulk profiles, not
  cell scores, so they appear as soon as the gene set is valid. Tick a row's checkbox to colour the plot by that table's specificity
  AUROC (cell type or region) and highlight the label; untick it to return to the previous colour mode. A table with a parent level
  shows it as the second column. The region table's **# sections** counts the sections with cells in each region's pseudobulk, so a
  region found in only one section, possibly a section-specific artifact, stands out.

## Data

Input is two files of the [Zenodo deposit](https://zenodo.org/records/20682863), saved in `data/` (not version controlled), plus two
[Allen Brain Cell Atlas](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html) parcellation tables and the atlas's 25 µm
annotation volume, which are downloaded into `data/Allen-CCF-2020/` when missing:

| File | Contents |
| --- | --- |
| `cellxgene_WT_labels.h5ad` | 146,863 cells of the whole-transcriptome samples 1–4 × 22,930 genes, raw counts (CSR); `obs`: `sample`, `cell_label`, `ccf_z`, `ccf_y`, `annot`, `global_x`, `global_y`, `volume`, `num_z`; `var`: `codebook`; `layers/logNormCounts` |
| `cellxgene_ISO_labels.h5ad` | 129,563 cells of the isoform samples 5–7 × 32,684 features (the 22,930 genes and 9,754 isoforms), raw counts (CSR); the same `obs` columns; `var`: `codebook` |
| `Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv`, `parcellation_term.csv` | `parcellation_index` → division, structure and substructure acronyms; their names and colours |
| `Allen-CCF-2020/20250331/annotation_25.nii.gz` | The CCF parcellation at 25 µm (528 × 320 × 456 voxels, each a parcellation index), for the section planes |

'WT' is the deposit's name for its whole-transcriptome samples (~23,000 genes), not wild type; its `ISO` set covers three other posterior
sections, imaged for the same genes and ~10,000 isoforms.

What the inspection showed:

- The WT file's `X` holds integer counts (1 to 784 per entry; median 834 transcripts and 592 detected genes per cell; 96.7% of cell × gene
  entries are zero), and `layers/logNormCounts` is exactly `log1p(X / total × 10⁴)`. The `raw` slot is not counts: it has non-integer
  values and 3 more genes (`Fbn2`, `Gm35188`, `H2-Eb1`), and is not used. Rows are not grouped by sample, and column indices are not sorted
  within rows. The anterior samples have more transcripts per cell (median ~1,020 against ~700).
- The ISO file holds one count matrix, without layers or a `raw` slot, so it is about a quarter the size of the WT file, which stores its
  110 million entries three times. Its genes are the WT file's 22,930 under the same names (`ISO_Codebook_1` and `_2`, 19,697 and 3,233
  genes, as in the WT codebooks); its isoforms (`ISO_Codebook_3`, Ensembl transcript IDs) hold 18.9 million of its 84.4 million entries and
  are dropped. Rows are grouped by sample, and `annot` is stored as floats, all whole numbers. Every cell has at least 50 gene transcripts.
- The ISO sections are shallower: across the 22,930 genes they are 97.8% zero, against 97.0% for the WT posterior sections, with median
  419, 408 and 283 genes detected per cell in samples 5, 6 and 7 against 538 and 479 in samples 3 and 4.
- `annot` is not an Allen structure ID: read as one, 32 of the WT file's 272 values are not in the ontology, and CA1 and the caudoputamen
  would have no cells. Read as ABC atlas parcellation indices, 89% of the `Hippocampus_CA1` cells fall in CA1, 95% of the `Hippocampus_DG`
  cells in DG, and the `Striatum_D1` cells in the caudoputamen, nucleus accumbens and olfactory tubercle.
- `ccf_z` spans 5.1–11.0 mm, on one side of the CCF midline (5.7 mm): each section is a hemisection, as the paper describes.
  `codebook` gives each gene's MERFISH run: 19,697 lower-expressed (Tier 1) and 3,233 higher-expressed (Tier 2) genes.
- Every cell has a cell type and at least 50 transcripts (the paper's filter), so all cells are kept. 3,148 cells (1.1%) lie outside
  the parcellation (`unassigned`).

## Preprocessing

```bash
./.conda-env/bin/python scripts/prepare_rttamp_cells.py --force
```

`--raw` (default `data/`) holds the two h5ad files, `--atlas` (default `data/`) the atlas tables and volume in the bucket's layout, and
`--reference` (default `data/reference`) the MGI marker lists; the tables, volume and lists are downloaded on the first run. Both count
matrices are read into memory: the run takes about 2 minutes, with peak memory of about 6 GB. `--force` replaces an existing output
directory once the new one is complete.

Steps:

1. **Cells**: every cell of both files, grouped by section and in source order within each. Display coordinates are `ccf_z, ccf_y` in
   micrometres. Neither file gives an anterior–posterior coordinate: the paper registered each section in 2D (affine, then B-spline) to
   one coronal slice of the Allen Reference Atlas, which the files do not record. Each cell's region was read from that slice's
   annotation, so the 25 µm coronal plane of the annotation volume whose labels match most of a section's cells is the slice it was
   registered to. 93–96% of cells match there (99.7–100% of those away from region boundaries), against 38–63% 250 µm away; one 25 µm
   step changes 5–13% of the labels. Samples 1 and 2, 4 and 7, and 5 and 6 were registered to the same atlas slices: cells within 3 µm
   of each other in the two sections of each pair carry the same region 98–99% of the time, as often as neighbouring cells within one
   section, against 36–56% between sections at different planes. Sections run anterior to posterior by that plane; ties keep sample
   order, which the cells support for 4 and 7 (see [Using the app](#using-the-app)).
2. **Genes**: the universe is every measured gene: 18,122 protein-coding and 4,516 lncRNA genes, 270 not in MGI and 11 others. The ISO
   file must measure the same genes; its isoform features are dropped. Non-coding genes are kept because the panel is already limited to
   protein-coding and lincRNA genes (Ensembl 102, intersected with the Allen scRNA-seq genes), unlike a whole-transcriptome universe such
   as Stereo-seq's. The paper also left out about 90 of the most abundant brain genes (at or above `Actb` in bulk RNA-seq, or with a
   maximum above that of `Gnrh1` across the Allen subclasses or of `Mef2c` across its clusters). Symbols are updated to current MGI
   symbols: 1,955 through unambiguous withdrawn symbols (`0610009B22Rik` → `Trappc2b`) and two that differ only in case (`ETV3L` →
   `Etv3l`, `ccdc198` → `Ccdc198`); no synonym guessing. The 95 unmapped genes (mostly clone-based names such as `AC103362.1`) and two
   ambiguous withdrawn symbols (`Adss`, `Lhfp`) keep their measured symbols. Eleven pairs of measured genes now name one MGI gene
   (`Exoc2` and `Gm29675`, `Plac9a` and `Plac9b`, …); their counts are summed before ranking, leaving 22,919 genes. Every cell's gene
   total is checked against its file after the merge.
3. **Ranks**: average ranks within each cell, with all undetected genes tied at the bottom, stored as doubled `uint16` values
   (format `sparse-doubled-average-ranks-v1`). Thirty-two random cells per section are checked against `scipy.stats.rankdata` and
   against a pairwise AUROC, and every batch is checked after a write/read round trip.
4. **Batches**: each section's cells are shuffled (seed 20261103) and split into equal batches of at most 16,384 cells
   (`--batch-cells`) and 12 million stored cell × gene entries (`--batch-entries`), in `progressive/<SECTION>/NNNNN.npz`.
5. **Annotations**: `annotations.json` holds labels and colours for each level. The source has no cell type colours. The groups fall
   into eight colour families, one per dark-surface categorical slot of a validated reference palette: cortex IT; thalamus,
   hypothalamus, midbrain and cholinergic; hippocampus; other cortical; inhibitory; striatum; glia; and vascular and ventricular (as
   the paper groups them in Fig. 2D). `Other` is grey. A family's types step through its hue at evenly spaced OKLCH lightness
   (0.50–0.80), dark to light in name order, so the cortex IT types go from `Cortex L2/3 IT` (dark) to `Cortex L6 IT` (light). Groups
   and classes take the colour of their most numerous type (`Neurons` is `Hypothalamus`'s orange), and every colour is checked for at
   least 3:1 contrast with the black plot. `annot` is each cell's parcellation index, read at three levels: division, structure and
   substructure (19, 122 and 271 labels). Regions use Allen CCF colours, with near-identical colours (such as cortical layers) spread
   in lightness so they stay distinguishable; `unassigned` is grey rather than the atlas black, the plot background.

Output (`data/rttamp-genes/`, 403 MB, of which 41 MB is pseudobulk; largest file 27 MB):

| File | Contents |
| --- | --- |
| `manifest.json` | Universe size, per-section cell counts and ranges, file set, mouse and estimated plane (`ccf_plane`, `ccf_x_um` and the share of cells it matches, `plane_match`), batch paths, input checksums, mapping counts, validation |
| `cells.parquet` | `cell_index`, `slice` (sample number), `source` (`WT` or `ISO`), `cell_id` (source name, as in the deposit's `cell_boundaries.parquet`), `h5ad_row` (row in its file), `rx, ry` (µm, `ccf_z, ccf_y`), `x, y` (µm, `global_x, global_y`), `volume`, `num_z`, `parcellation_index`, `cell_class`, `cell_group`, `cell_type`, `region_division`, `region_structure`, `region_substructure`, `draw_order`, `batch`, `batch_row`, `n_counts`, `n_genes_detected` |
| `genes.parquet` | `gene_index`, symbol, MGI ID, name, MGI feature type, `codebook` (the WT file's), merged `source_genes`, detection and total counts |
| `gene_detection.parquet` | Detected-cell counts per gene and section (for the sparsity note) |
| `pseudobulk/<SECTION>.npz` | For each annotation level: summed counts per (label, gene) for the labels present (`<level>_counts`, `uint32`), their global label indices (`<level>_index`, into `manifest.json` `pseudobulk.labels`) and cell counts (`<level>_cells`) |
| `gene_mapping.parquet` | All 22,930 measured genes with their MGI mapping (`mapping_method`), name, feature type and codebook |
| `annotations.json` | Titles, labels and colours for all six annotation levels, the names of the cell type abbreviations and of the regions |

| Section | File set | Mouse | Plane | Cells | Size | Batches |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 1 | WT | WT mouse | 4.14 mm | 31,385 | 4.7 × 6.3 mm | 3 |
| 2 | WT | WT mouse | 4.14 mm | 28,121 | 4.3 × 6.3 mm | 3 |
| 4 | WT | WT mouse | 7.51 mm | 42,364 | 5.5 × 7.4 mm | 3 |
| 7 | ISO | ISO mouse 2 | 7.51 mm | 40,569 | 5.7 × 7.3 mm | 3 |
| 5 | ISO | ISO mouse 1 | 7.69 mm | 47,649 | 5.7 × 7.4 mm | 3 |
| 6 | ISO | ISO mouse 1 | 7.69 mm | 41,345 | 5.7 × 7.2 mm | 3 |
| 3 | WT | WT mouse | 7.84 mm | 44,993 | 5.9 × 7.4 mm | 3 |
| **Total** | | | | **276,426** | | **21** |

## Binned view

```bash
./.conda-env/bin/python scripts/prepare_rttamp_bins.py --bin-um 178
```

Each section's cells are grouped by centroid into the squares of one lattice shared by all sections (they share the CCF frame), so
every bin is the same size and smaller sections simply have fewer bins. A bin's counts are the summed transcripts of its cells, ranked
within the bin as for cells. `--bin-um` sets the side; `--bins N` instead fits N squares across the widest section. 178 µm is the bin
size of the [Stereo-seq version](https://github.com/leonfrench/polygenic_brain_Han_StereoSeq), the other whole-transcriptome atlas, so
both bin the same area; its rule (the longest section extent over 40) would give 185 µm here. It gives 6,857 bins with a median of 39
cells (5,923 with at least 20 cells; 705 in sample 2 to 1,106 in sample 3). The cell dataset stores ranks, not counts, so both h5ad
files are read again, with genes merged as `gene_mapping.parquet` records; every cell's counts and every section's gene totals are
checked against the cell dataset, and 64 bins per section against `scipy.stats.rankdata`. The build takes about 80 seconds with peak
memory of up to 6.4 GB. It also writes per-bin counts and each cell's bin (both gitignored); the app reads only `manifest.json`,
`bins.parquet` (grid position, centre, cells, counts, and the most common label of each level with its share) and
`ranks/<SECTION>.npz` (the cell batch format, one row per bin), so only those are committed in `data/rttamp-genes-bins178um/` (85 MB),
with the manifest's file list trimmed to them.

When the bins were built from the app's cell dataset, **Colour cells by** gains *178µm bin AUROC*: the plot shows the bins instead
of the cells, each with a gene-set AUROC from its within-bin ranks, on one absolute scale for all sections from the lowest to the
highest bin, without clipping. Bins with fewer than 20 cells are grey and not scored. The tooltip gives the bin's cells and grid
position and the most common cell type and structure of its cells, plus the sidebar's cell type and region levels when they differ
(each first within its kind). Bin scores read one section's ranks at a time. If scoring fails, the page shows a reference that
matches the error in the server log. The tables are unchanged. `RTTAMP_APP_BINS=/path/to/bins` selects other bins; only the
per-section layout of the current script is offered.

## Methods

**Cell scores (plot).** For k target genes in a universe of G genes, a cell's score is
`AUROC = (sum of target ranks − k(k+1)/2) / (k(G − k))`, the probability that a target gene outranks a
non-target gene in that cell, with half credit for ties.

**Expression specificity (tables).** For the chosen level, each label's pseudobulk sums the transcript counts of its
cells across all sections and is converted to log₂(CPM + 1). Each gene is z-scored across the labels, removing its
average expression (genes that do not vary get z = 0). Within each label, genes are ranked by z-score (average ranks for
ties), and the AUROC compares the target genes with the rest of the universe: above 0.5 means the gene set is more
specific to that label than other genes. With a custom background, genes are ranked within that universe. Labels with
fewer than 20 cells are not profiled. Two-sided Mann–Whitney p-values over genes (normal approximation with tie and
continuity correction, matching `scipy.stats.mannwhitneyu`) are adjusted with Benjamini–Hochberg within each table
(pFDR).
Genes are not independent (co-expression), so the p-values are exploratory. This follows the
[SEA-AD polygenic tester](https://github.com/leonfrench/SEA_AD_multi_polygenic_app), which ranks genes by pseudobulk
z-scores in the same way; its p-values omit the tie correction.

## Tests

```bash
./.conda-env/bin/python -m pytest tests -q
node tests/test_cells.cjs
```

These cover rank encoding and scoring against SciPy and pairwise definitions, gene mapping (current, withdrawn and case-only symbols)
and the merging of genes that share a symbol, count reading and the dropping of isoform columns, the annotation volume reader and the
section planes, pseudobulk sums, the cell type groups and colours, specificity z-scores and AUROCs against `mannwhitneyu`, ortholog
mapping, the grid layouts, the score transport protocol, the canvas helpers, bin scoring order and square bin positions, and end-to-end
app runs on `data/rttamp-genes` and its bins in `data/rttamp-genes-bins178um`.

## Data credit

Cohen et al. *Whole-transcriptome-scale isoform-resolved spatial imaging of single cells in tissues.*
Cell 189, 6189–6207 (2026). [doi:10.1016/j.cell.2026.06.027](https://doi.org/10.1016/j.cell.2026.06.027).
Data, cell types and CCF registration: [Zenodo](https://zenodo.org/records/20682863) (CC BY 4.0).
Region names, colours and annotation volume: [Allen Brain Cell Atlas](https://brain-map.org/bkp/explore/abc-atlas) parcellation
(Allen-CCF-2020). Gene symbols: [MGI](https://www.informatics.jax.org/). Orthologs: HomoloGene.
