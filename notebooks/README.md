# Research notebooks

These notebooks document the development of tree traversal, distance
calculation, and laboratory-passage categorization. They retain their original
cell order, outputs, and study-specific paths. Running them requires the
corresponding local inputs and, in some cases, state from earlier cells.

| Notebook | Purpose |
|---|---|
| [get_data_from_tree_aa.ipynb](get_data_from_tree_aa.ipynb) | Extract successive amino-acid substitutions from the H3N2 tree and ancestral sequences |
| [get_data_from_three.ipynb](get_data_from_three.ipynb) | Earlier extraction workflow, with H1N1 inputs in the retained version |
| [dist_old.ipynb](dist_old.ipynb) | Inspect individual nodes, substitutions, and distance calculations |
| [passage_phenotype_categories.ipynb](passage_phenotype_categories.ipynb) | Convert passage metadata into clinical/laboratory categories for the passage-association analysis |

The passage-category notebook prepares the phenotype table used upstream of
EpiStat. It requires restricted strain metadata. The association test was run
outside this repository.

For extraction on new data, the command-line script provides input validation
and a defined CSV output. [Extraction script](../scripts/extract_substitution_pairs.py) ·
[Input format](../data/README.md) · [Run the analysis](../README.md#run-the-analysis)
