# Reconstructing the bundled real-data examples

Run this from the package source root.

    Rscript inst/scripts/rebuild_example_data.R

The script requires the container packages listed in DESCRIPTION, Matrix,
scran, and SeuratObject. SeuratObject is used only during preparation to read
the published PBMC source object. The saved examples are Bioconductor
containers. Loading the examples and running the vignette require no Seurat
objects or downloads.

The script writes microbiome_tse.rds, umi_sce.rds, nonumi_sce.rds, and
franzosa_example_mae.rds for explicit assignment with readRDS(). It also retains
the combined tweedieverse_data.RData bundle for existing users. These files
contain the same objects and have the same per-object licences.

The optional arguments are an output directory, a cached HuMMANet RDS, and a
cache directory containing the four other source files named below.

    Rscript inst/scripts/rebuild_example_data.R /tmp/examples /path/to/franzosae_2019.rds /path/to/cache

Downloads are confined to this manual preparation script. The published PBMC
object is approximately 274 MB, while the packaged subset is compact.
Every downloaded file is checked against its recorded MD5. The GitHub sources
are pinned to revisions in the script. Subset selection is deterministic and
uses category labels and count variance, without fitting association models.

## Microbiome baseline example

The workflow follows the original [TweedieLabs microbiome tutorial](https://github.com/himelmallick/TweedieLabs/blob/master/R/Tweedieverse_microbiome.Rmd).
The source is the pair of HMP2 example tables distributed in MaAsLin2 at
revision 39eda3be7d801f0eb594eb7b09b26fad7d5c5875.

- HMP2_taxonomy.tsv MD5 is 78231446ddedd5c080cea53bbf37ac38.
- HMP2_metadata.tsv MD5 is 30dce926858e6feda6fbde8a149fd2b4.
- Both files are read from inst/extdata in the pinned
  [MaAsLin2 repository](https://github.com/biobakery/Maaslin2/tree/39eda3be7d801f0eb594eb7b09b26fad7d5c5875/inst/extdata).
- The original study is Lloyd-Price et al. (2019),
  [Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases](https://doi.org/10.1038/s41586-019-1237-9).

Order metadata by collection and sample ID, take the first sample per subject,
require collection < 3, and retain complete age, antibiotic, and diagnosis data.
Take the first 16 eligible samples in each of nonIBD, CD, and UC. Retain the 20
features with highest abundance variance, breaking ties by feature name. Preserve
the original relative-abundance values. Wrap the resulting 20-by-48 matrix
and its metadata in microbiome_tse, a TreeSummarizedExperiment.

The source package declares MIT + file LICENSE, copyright 2018 Harvard School
of Public Health. Its notice is retained in inst/COPYRIGHTS.

## PBMC3k UMI example

The source object is the public processed PBMC3k object used by the
[TweedieLabs UMI tutorial](https://github.com/himelmallick/TweedieLabs/blob/master/R/Tweedieverse_scRNASeq_UMI.Rmd),
available from the [published download](https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=1).
Its MD5 is ca57284952b2a17553cee5b87a71e7a7.

Read the RNA assay's raw counts and published cell identities. Combine Naive
CD4 T and Memory CD4 T as CD4, and retain CD8 T as CD8. Select the first 30
barcodes in sorted order from each group. Calculate library sizes using all
13,714 source genes in these cells, then retain the 20 genes with highest
count variance. Store raw counts and CellType/library_size metadata in umi_sce,
a SingleCellExperiment. The vignette uses these full-library sizes as offsets.

The counts originate from [10x Genomics PBMC3k](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-0-0),
released under CC BY 4.0. Cell labels and upstream quality control come from
the [Satija Lab PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

## Brain non-UMI example

The source is brain_scRNAseq.rda from SC2P revision
bdbb526cb81d598da99b63f55cc64c826863b954, MD5
15562abf4745bc7aa8147f878b210fe9. It supplies the count matrix Y and sample
annotations design used by the
[TweedieLabs non-UMI tutorial](https://github.com/himelmallick/TweedieLabs/blob/master/R/Tweedieverse_scRNASeq_nonUMI.Rmd).
The original study is Darmanis et al. (2015),
[A survey of human brain transcriptome diversity at the single cell level](https://doi.org/10.1073/pnas.1507125112),
GEO accession [GSE67835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835).

Assign design row names as count column names and remove unused factor levels.
Select the first 30 sorted cell IDs from each of astrocytes and oligodendrocytes.
Compute SCRAN factors using all 10,483 genes in the retained cells, with cell
classes as clusters and pool sizes 10 and 20. Retain the 20 genes with highest
count variance. Store raw read counts, original design variables, and
source_size_factor in nonumi_sce, a SingleCellExperiment.

SC2P declares GPL without specifying a version. Its dataset licence is recorded
separately from the package's MIT code licence. A copy of GNU GPL version 2 is
included in inst/extdata/SC2P-GPL-2.txt. Credit SC2P contributors Zhijin (Jean)
Wu, Jamie (Yi) Zhang, and Hao Wu, and the original Darmanis et al. study.

## Franzosa et al. data

The verified source is the HuMMANet v1 archive by Shivangi Verma, Nalin Arora,
Himel Mallick, and Tarini Shankar Ghosh (2026).

- Archive DOI and version-specific record
  https://doi.org/10.5281/zenodo.22053943
- Source file
  https://zenodo.org/records/22053943/files/franzosae_2019.rds?download=1
- Published file MD5
  9ab747a9c36003e9b76f4d9d58d67e77
- Data licence declared by the archive
  Creative Commons Attribution 4.0 International, CC BY 4.0
- Original study
  Franzosa EA et al. (2019), "Gut microbiome structure and metabolic activity
  in inflammatory bowel disease", Nature Microbiology 4, 293-305.
  https://doi.org/10.1038/s41564-018-0306-4

The archive contains processed study tables. Its upstream curation and
metabolite harmonization are described by HuMMANet
(https://doi.org/10.64898/2026.08.24.746727). This script begins with those
processed tables.

The reconstruction script performs the following steps.

1. Use sample_id as the row names of MetadataProfile and the assay tables.
2. Remove sample_id from the assay values and transpose to
   features-by-samples matrices.
3. Align both assays to the 155 metadata samples. Preserve abundance values
   without additional normalization, transformation, or imputation.
4. Create microbial row annotations from feature names and match metabolite
   annotations to metaboliteIdentifierMapping by Query_Name.
5. Wrap the 711 microbial features and 339 metabolites in
   SummarizedExperiment objects and a MultiAssayExperiment.
6. Save this object inside the list FranzosaE_2019_data as
   FransozaE_2019_MAE.RData. The historical filename spelling is retained.
7. Select the explicitly listed 16 samples and four features per layer for the
   vignette. Keep diseaseCat_clean and Age in primary colData, while retaining
   the layer annotations and metadata. Save franzosa_example_mae alongside the
   three real-data containers in tweedieverse_data.RData.

The reconstructed assay values, annotations, and metadata were checked against
both existing bundled RData files. This recovers a reproducible construction
from a verified public source. The Git history does not record the original
download command used in 2026.

See inst/extdata/DATA_LICENSES.md for the licence of each bundled object.
