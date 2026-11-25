# Exploratory R scripts

These files preserve the interactive analysis history behind the thesis.
Several use personal absolute paths, depend on objects created earlier in an
R session, or contain alternative analyses. The reusable model-fitting entry
point accepts an input CSV and an output path.
[Model script](../fit_radicality_models.R) · [Run the analysis](../../README.md#run-the-analysis)

## File guide

| File | Purpose |
|---|---|
| [BD.r](BD.r) | Load and inspect the study datasets and their filtered versions |
| [make_radicalisation.r](make_radicalisation.r) | Calculate amino-acid radicality measures |
| [make_models.r](make_models.r) | Explore group filtering and model specifications |
| [grouped_regressions.R](grouped_regressions.R) | Compare regressions across annotated site groups |
| [filter_senescent_sites.r](filter_senescent_sites.r) | Prepare senescence-related site subsets |
| [filter_exclude.r](filter_exclude.r) | Exclude sites associated with laboratory passage |
| [passage_phenotype_states.r](passage_phenotype_states.r) | Reconstruct clinical/laboratory states on internal tree nodes |
| [zamena.r](zamena.r) | Compare radicality at passage-associated sites and other sites |
| [deep.r](deep.r) | Inspect distance distributions in senescent-site subsets |
| [hipervar.r](hipervar.r) | Explore repeatedly observed groups and their distance/radicality variation |
| [work_with_data_from_tree.r](work_with_data_from_tree.r) | Compare extraction outputs and distance definitions |

## Known issues in the historical analysis

`make_models.r` contains a filter that compares `anc_aa` with a vector of site
numbers. The reusable model script applies the minimum group-size threshold
to the appropriate grouping variables. A comment in the historical script
identifies the affected block.

These details matter when tracing the original tables and comparing them with
new runs. The results guide describes the reported findings and model-specific
sample sizes. [Results](../../results/README.md)
