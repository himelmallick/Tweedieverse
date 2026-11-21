# Reconstructing the bundled real-data examples

Run this from the package source root.

    Rscript inst/scripts/rebuild_example_data.R

The script requires the container packages listed in DESCRIPTION, Matrix,
scran, and a separately installed copy of SeuratObject. SeuratObject is a
manual preparation prerequisite, not a package dependency, and is used to read
the published PBMC source object. The saved examples are Bioconductor
containers. Loading the examples and running the vignette require no Seurat
objects or downloads.

The script writes microbiome_tse.rds, umi_sce.rds, nonumi_sce.rds, and
franzosa_example_mae.rds for explicit assignment with readRDS(). It also retains
the combined tweedieverse_data.RData bundle for existing users. These files
contain the same objects and have the same per-object licences.

The optional arguments are an output directory, a cached HuMMANet RDS, and a
cache directory containing the five other source files named below.

    Rscript inst/scripts/rebuild_example_data.R /tmp/examples /path/to/franzosae_2019.rds /path/to/cache

Downloads are confined to this manual preparation script. The published PBMC
object is approximately 274 MB. The Tabula Muris source archive is approximately
304 MB; only its brain count table is extracted. The packaged subsets are compact.
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

The replacement source is the original Tabula Muris Smart-seq2 FACS dataset,
Figshare version 8, [DOI 10.6084/m9.figshare.5829687.v8](https://doi.org/10.6084/m9.figshare.5829687.v8).
The record explicitly declares CC BY 4.0 for the deposited data and annotations.
A copy of its metadata is in sources/tabula-muris-5829687-v8.json.

The pinned files are:

- FACS.zip, Figshare file 10700143, MD5
  ef6e44ca2170155052591b5ab9d5d778.
- annotations_facs.csv, Figshare file 13088129, MD5
  8ba7b81d91fe569cea07c384c2d34912.

Download URLs are https://ndownloader.figshare.com/files/10700143 and
https://ndownloader.figshare.com/files/13088129. The script verifies both MD5s.
Extract only FACS/Brain_Non-Myeloid-counts.csv from the archive. The source
record describes these as gene-count tables from Smart-seq2, a full-length
non-UMI protocol; these are not counts from the project's droplet assay.

Filter annotations to Brain_Non-Myeloid and the cell ontology classes
astrocyte and oligodendrocyte. Among mice with at least 30 cells in each class,
choose the first sorted mouse ID (3_10_M), then the first 30 sorted cell IDs
per class. Read those 60 columns from the count table, retaining all source
genes. Compute SCRAN factors with the cell classes as clusters and pool sizes
10 and 20. Select the 20 genes with highest count variance using the common
variable_features() helper. Store unchanged source counts, the selected
annotations, a celltype factor, source_size_factor, and source/licence metadata
in nonumi_sce. This is an illustrative cell-type comparison within one mouse.

Credit Tabula Muris Consortium, James Webber, Joshua Batson, and Angela Pisco
(2018), the authors listed on the Figshare record, and cite the original study:
Tabula Muris Consortium et al. (2018), "Single-cell transcriptomics of 20 mouse
organs creates a Tabula Muris", https://doi.org/10.1038/s41586-018-0590-4.
Redistribution retains attribution, the CC BY 4.0 licence link, and a description
of subsetting and added size factors. No SC2P data are used in this example.

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
5. Wrap the original 711 microbial features and 339 metabolites in
   SummarizedExperiment objects for selecting the vignette subset.
6. Select the explicitly listed 16 samples and four features per layer. Keep
   diseaseCat_clean and Age in primary colData, while retaining layer
   annotations and metadata, and combine the subsets in a MultiAssayExperiment.
7. Save franzosa_example_mae.rds and include the same object alongside the
   three single-assay containers in tweedieverse_data.RData. The unused full
   MultiAssayExperiment bundle is no longer distributed or generated.

The reconstructed assay values, annotations, and metadata were checked against
the bundled vignette objects. This provides a reproducible construction from a
verified public source. The Git history does not record the original download
command used in 2026.

See inst/extdata/DATA_LICENSES.md for the licence of each bundled object.
