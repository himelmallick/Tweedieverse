# Bundled data provenance and licences

## tweedieverse_data.RData

This file contains objects with the following sources and licences. Each object
is also supplied in its own file (microbiome_tse.rds, umi_sce.rds,
nonumi_sce.rds, and franzosa_example_mae.rds), under the same licence.

| Object | Source | Licence |
| --- | --- | --- |
| microbiome_tse | HMP2 baseline data distributed with MaAsLin2 | MIT, copyright 2018 Harvard School of Public Health |
| umi_sce | 10x PBMC3k raw counts with published Satija Lab cell labels | CC BY 4.0 for the 10x counts |
| nonumi_sce | Tabula Muris Smart-seq2 FACS brain counts, Figshare v8 | CC BY 4.0 |
| franzosa_example_mae | A subset of the HuMMANet v1 FranzosaE_2019 data | CC BY 4.0 |

The Franzosa vignette object retains 16 samples and four features from each
assay. It preserves the source assay values. Container conversion, annotation
matching, and subsetting are documented in inst/scripts/README.md and
inst/scripts/rebuild_example_data.R.

## Attribution for the HuMMANet-derived example

Credit Shivangi Verma, Nalin Arora, Himel Mallick, and Tarini Shankar Ghosh
(2026), HuMMANet v1, https://doi.org/10.5281/zenodo.22053943.
The record metadata at https://zenodo.org/api/records/22053943 declares
the licence identifier cc-by-4.0.

CC BY 4.0 terms are available at
https://creativecommons.org/licenses/by/4.0/.

Cite the original study by Franzosa EA et al. (2019), "Gut microbiome structure
and metabolic activity in inflammatory bowel disease", Nature Microbiology 4,
293-305, https://doi.org/10.1038/s41564-018-0306-4.

## Other source notices

The HMP2 subset comes from the MIT-licensed MaAsLin2 distribution, copyright
2018 Harvard School of Public Health. Source and permission details are retained
in inst/COPYRIGHTS. Credit the original Lloyd-Price et al. (2019) HMP2 study,
https://doi.org/10.1038/s41586-019-1237-9.

The PBMC counts come from 10x Genomics' "3k PBMCs from a Healthy Donor" dataset,
which declares CC BY 4.0. Published cell identities and quality control are
attributed to the Satija Lab PBMC3k tutorial and the processed object linked in
inst/scripts/README.md. The package subset retains raw counts, selects CD4/CD8
cells, and calculates full-library sizes before gene selection.

The non-UMI example is a subset of the original Tabula Muris Smart-seq2 FACS
count tables and annotations, Figshare version 8:
https://doi.org/10.6084/m9.figshare.5829687.v8.
The original record explicitly declares CC BY 4.0 for these deposited files.
The source metadata snapshot is in
inst/scripts/sources/tabula-muris-5829687-v8.json; its licence URL is
https://creativecommons.org/licenses/by/4.0/.

Credit the record's authors: Tabula Muris Consortium, James Webber, Joshua Batson,
and Angela Pisco (2018). Cite Tabula Muris Consortium et al. (2018), "Single-cell
transcriptomics of 20 mouse organs creates a Tabula Muris",
https://doi.org/10.1038/s41586-018-0590-4.
The package subset retains 30 astrocytes and 30 oligodendrocytes from mouse
3_10_M, selects 20 genes by count variance, and adds SCRAN factors calculated
before gene selection. Raw source count values are unchanged. The cells come
from the FACS Smart-seq2 assay, not the droplet assay. The counts and annotations
retain CC BY 4.0, separately from the package code's MIT licence.

The package's MIT licence applies to its code. The bundled data objects retain
the source-specific licences listed above. Full source URLs, revision IDs,
checksums, and exact selection procedures are in inst/scripts/README.md and
inst/scripts/rebuild_example_data.R.
