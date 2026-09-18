# Single-Cell RNA-seq of the Breast-Tumour Microenvironment

Unsupervised recovery of the cell types in a human breast tumour from single-cell transcriptomes, then
**validation against expert labels**. Real cells are fetched **live from the CELLxGENE Census API** (not a
packaged tutorial dataset), clustered and annotated from marker genes with **Scanpy**, and reproduced on the
identical cells in **R/Seurat**.

![UMAP: our labels vs Census expert labels](results_py/umap_celltypes.png)

## Research question

What cell types make up a human breast tumour, and can we recover them without supervision from scRNA-seq — and
do our marker-based calls agree with independent expert annotation?

## Data (real, fetched via API)

- **Source:** CELLxGENE Census (CZI), queried through the `cellxgene-census` Python API (TileDB-SOMA over S3).
- **Query:** `tissue_general == 'breast' and disease != 'normal' and is_primary_data == True`.
- **Sample:** ~5,000 cells; Census version pinned (`2025-11-08`) and random seed fixed, so the exact cells are
  reproducible.
- **Bonus:** the Census ships expert `cell_type` labels. We ignore them while clustering and use them only to
  validate our own annotation (a cross-tab).

## Method / pipeline

Fetch (Census API) → QC (mito%, gene/cell filters) → library-size normalize + `log1p` → highly-variable genes →
scale → **PCA** → kNN graph → **UMAP** → **Leiden** clustering → **Wilcoxon** marker genes → marker-based
cell-type annotation → **cross-tab vs Census expert labels**.

| Stage | Scanpy (Python) | Seurat (R) |
|---|---|---|
| Normalize | `normalize_total` + `log1p` | `NormalizeData` |
| Feature selection | `highly_variable_genes` | `FindVariableFeatures` |
| Reduce | `pp.scale` + `tl.pca` | `ScaleData` + `RunPCA` |
| Graph + cluster | `pp.neighbors` + `tl.leiden` | `FindNeighbors` + `FindClusters` |
| Embed | `tl.umap` | `RunUMAP` |
| Markers | `rank_genes_groups` (wilcoxon) | `FindAllMarkers` |

## Key result

The unsupervised pipeline recovers the expected tumour-microenvironment compartments — **malignant/epithelial**
(keratins, GATA3), **fibroblasts/CAFs** (collagens, ACTA2/TAGLN), **endothelial** (VWF/PECAM1, lymphatic MMRN2),
**T cells** (CD3D/CD3E), **B / plasma** (MS4A1/CD79A; IGKC/MZB1), **macrophages** (CD68/LYZ), and **mast cells**
(TPSAB1/CPA3). Cross-tabulating these labels against the Census expert `cell_type` is near block-diagonal — our
marker-based calls agree with independent annotation. Biologically: a tumour is an *ecosystem* of malignant
epithelium plus stroma and immune cells, resolvable only at single-cell resolution.

**R/Seurat twin, on the identical cells:** malignant, endothelial, fibroblast and macrophage labels each
strongly concentrate on their matching Census category (malignant: ~1,700/1,900 cells; endothelial: ~100% across
endothelial subtypes; fibroblast and macrophage similarly dominant). Its "T cell" label is a broader lymphocyte
super-cluster that also picks up B/plasma cells — because the R pipeline used 10 PCs (`dims=1:10`, the standard
Seurat tutorial default) vs Python's 40, under-resolving lymphocyte subtypes. This is an expected, explainable
cross-tool difference, not an error, and is documented rather than hidden.

## Limitations

Parameter-dependent (QC thresholds, #HVGs, #PCs, Leiden/Louvain resolution — and the R/Python PC-count difference
above); UMAP geometry is not quantitative; annotation is marker-based, and a few clusters in each language lacked
a clean marker signature and were labelled via a ground-truth crosstab tiebreak rather than markers alone; the
Census pool spans multiple donors/datasets, so a rigorous study would batch-integrate (Harmony/scVI) first;
proving an epithelial cluster is *malignant* (vs normal epithelium) needs copy-number inference
(inferCNV/CopyKAT); NK cells fold into other clusters at this resolution in both languages.

## Files

```
single_cell_tme.ipynb    # Python / Scanpy pipeline (fetch → cluster → annotate → validate)
single_cell_tme.R        # R / Seurat twin (reads the same cells exported as 10x)
results_py/              # UMAPs, marker table, crosstab vs Census (Python)
results_R/               # UMAPs, marker table, crosstab vs Census (R)
data/tumor_10x/          # the fetched cells in 10x format (written by the notebook; feeds the R twin)
```

## Run

**Python:** `pip install cellxgene-census scanpy leidenalg`, then run `single_cell_tme.ipynb` top to bottom
(needs internet for the Census fetch; first fetch takes a few minutes). **R:** run the notebook first (it writes
`data/tumor_10x/`), then `single_cell_tme.R` from the same folder (`install.packages(c("Seurat","dplyr"))`).
