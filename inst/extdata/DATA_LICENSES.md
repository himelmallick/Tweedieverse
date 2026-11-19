# Bundled data provenance and licences

## tweedieverse_data.RData

This file contains objects with the following sources and licences. Each object
is also supplied in its own file (microbiome_tse.rds, umi_sce.rds,
nonumi_sce.rds, and franzosa_example_mae.rds), under the same licence.

| Object | Source | Licence |
| --- | --- | --- |
| microbiome_tse | HMP2 baseline data distributed with MaAsLin2 | MIT, copyright 2018 Harvard School of Public Health |
| umi_sce | 10x PBMC3k raw counts with published Satija Lab cell labels | CC BY 4.0 for the 10x counts |
| nonumi_sce | SC2P brain count data, GEO GSE67835 | GPL, as declared by SC2P |
| franzosa_example_mae | A subset of the HuMMANet v1 FranzosaE_2019 data | CC BY 4.0 |

## FransozaE_2019_MAE.RData

This file contains FranzosaE_2019_data, a list holding one MultiAssayExperiment
with 155 samples, 711 microbial features, and 339 metabolites. It is reconstructed
from the HuMMANet v1 data archive, distributed under CC BY 4.0.

The vignette subset retains 16 samples and four features from each assay.
Both files preserve the source assay values. Container conversion, annotation
matching, and subsetting are documented in inst/scripts/README.md and
inst/scripts/rebuild_example_data.R.

## Attribution for the HuMMANet-derived objects

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

The non-UMI counts are distributed in the GPL-licensed SC2P package by Zhijin
(Jean) Wu, Jamie (Yi) Zhang, and Hao Wu. SC2P's DESCRIPTION declares GPL without
a version number. GNU GPL version 2 is included in SC2P-GPL-2.txt. Credit
Darmanis et al. (2015), https://doi.org/10.1073/pnas.1507125112, and GSE67835.
The package subset selects astrocytes and oligodendrocytes and calculates SCRAN
factors before selecting genes.

The package's MIT licence applies to its code. The bundled data objects retain
the source-specific licences listed above. Full source URLs, revision IDs,
checksums, and exact selection procedures are in inst/scripts/README.md and
inst/scripts/rebuild_example_data.R.
