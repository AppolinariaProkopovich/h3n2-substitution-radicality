# Local data layout

No study data are distributed with this repository. Use this directory only
for inputs obtained under the appropriate access terms.

## Suggested layout

```text
data/
  raw/
    tree.newick
    ancestral_sequences.fasta
    strain_metadata.tsv
  processed/
    substitution_pairs.csv
    substitution_pairs_radicality.csv
```

`extract_substitution_pairs.py` expects a Newick tree in which every leaf and
internal node is named, plus an in-frame nucleotide FASTA containing one
aligned sequence for every tree node. Text after whitespace or `/` in a FASTA
header is removed when matching sequence identifiers to node names.

The extraction step writes the following columns.

| Column | Meaning |
|---|---|
| `site` | One-based amino-acid position |
| `ancestor_node` | Node at which the parental amino acid appeared |
| `grandanc_codon`, `grandanc_aa` | State before the first substitution |
| `anc_codon`, `anc_aa` | State produced by the first substitution |
| `der_codon`, `der_aa` | State produced by the next substitution |
| `distance` | Distance from the first substitution node to the midpoint of the next substitution branch |
| `distance_grand` | Equivalent distance from the pre-first-substitution node |

The radicality-annotation step adds two columns.

| Column | Meaning |
|---|---|
| `gr_rad_anc` | Grantham distance from `anc_aa` to `der_aa` |
| `miyata_rad` | Miyata distance from `anc_aa` to `der_aa` |

`distance` is expressed in the branch-length units of the input tree: expected
substitutions per site for the FastTree trees of D1 and D2, and years for the
TreeTime tree of D3.

Subset labels for clustered substitutions, senescent sites, and sites
associated with laboratory passage come from upstream analyses and are not
inferred by the reusable scripts in this repository. The passage-associated
sites are the 9 concordant sites and 15 R-sites listed in the results guide.
[Passage-associated sites](../results/README.md#laboratory-passage-association-analysis)

## Starting from an existing table

The R model script requires `distance`, `site`, `anc_aa`, `gr_rad_anc`, and
`miyata_rad`. Distance and radicality columns must be numeric. Site and
parental amino acid are treated as categorical variables.

Prepare a separate table for each dataset and site subset before fitting.
The script fits all three model specifications to each supplied table.
[Model fitting](../README.md#4-fit-the-models)

Raw inputs and processed tables are excluded from version control.
[Data policy](../DATA_POLICY.md)
