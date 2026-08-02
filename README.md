# Radicality of amino-acid substitutions in influenza A/H3N2

This repository accompanies a 2026 MSc thesis on whether successive amino-acid
substitutions in influenza A/H3N2 hemagglutinin become more or less chemically
radical with evolutionary distance. Radicality is measured with the Grantham
and Miyata scales and tested with linear and linear mixed models.

**Main result:** the negative trend seen in simple regressions is not a general
property of natural H3N2 evolution. It is driven by differences among protein
sites and by laboratory-passaged sequences. In two clinical-only datasets,
clustered substitutions show a reproducible positive trend after site-level
heterogeneity is included, although the effect becomes non-significant in the
more demanding model that nests parental amino acid within site.

The full thesis is available in Russian:
[`thesis/Prokopovich_2026_MSc_thesis_ru.pdf`](thesis/Prokopovich_2026_MSc_thesis_ru.pdf).

## Research question

For two successive nonsynonymous changes at one HA site, the analysis records:

- the amino acid before the first change (`grandanc_aa`);
- the amino acid produced by the first change (`anc_aa`);
- the amino acid produced by the next change (`der_aa`);
- the phylogenetic distance between the two changes;
- the Grantham and Miyata distance from `anc_aa` to `der_aa`.

The fixed effect of phylogenetic distance is compared across three models:

1. a simple linear regression;
2. a mixed model with site as a random intercept;
3. a mixed model with parental amino acid nested within site.

The third specification tests whether an apparent temporal trend remains after
accounting for different baseline radicalities among `site x parental amino
acid` groups.

## Data

Three HA datasets were analysed in the thesis:

| Dataset | Sequences | Tree / ancestral reconstruction | Substitution pairs | Sites |
|---|---|---|---:|---:|
| D1 | Clinical and laboratory-passaged | FastTree / MEGA | 2,492 | 125 |
| D2 | Clinical only | FastTree / MEGA | 1,293 | 53 |
| D3 | Clinical only, 2004 to September 2025 | TreeTime | 2,744 | 95 |

D1 was also used for a phylogenetic GWAS of laboratory passage. It identified
9 concordant sites, where changes coincide with transition to laboratory
culture, and 15 R-sites, where more than 60% of changes occur on laboratory
branches. Removing those 24 sites weakened but did not eliminate the artefact,
suggesting that laboratory adaptation has a distributed effect.

The sequence data originate from GISAID and cannot be redistributed here.
See [DATA_POLICY.md](DATA_POLICY.md) and [data/README.md](data/README.md) for
the local input layout and table schema.

## Results

| Analysis | D1: clinical + laboratory | D2: clinical | D3: clinical |
|---|---|---|---|
| All substitutions | Negative trend remains with site effect; disappears in nested model | No negative trend; weak positive Miyata effect with site | Simple Miyata trend disappears with site effect |
| Clustered substitutions | Negative simple-model trend disappears with random effects | Positive trend remains with site effect | Positive trend remains with site effect |
| Senescent sites | Negative trend disappears in nested model | No significant trend | Simple negative trend disappears with site effect |

The clustered-substitution signal is the only result reproduced in both
clinical datasets and with both radicality metrics. Its loss of significance in
the nested model leaves two plausible interpretations: residual confounding by
parental amino-acid context, or loss of power in small within-site groups.

## Workflow

The reusable part of the repository implements the downstream workflow:

```text
named Newick tree + ancestral nucleotide FASTA
                        |
                        v
           consecutive substitution pairs
                        |
                        v
             Grantham + Miyata distances
                        |
                        v
       linear model / site LMM / nested LMM
```

Tree construction, ancestral reconstruction, clustered-site calls, senescent-
site calls, and phylogenetic GWAS were performed upstream. The exploratory
notebooks and R sessions document those study-specific steps, but they still
refer to private inputs and are retained as provenance rather than presented as
a portable pipeline.

## Repository layout

```text
R/fit_radicality_models.R        reusable fixed and mixed-model analysis
R/exploratory/                   original interactive R research scripts
notebooks/                       original extraction and GWAS notebooks
scripts/extract_substitution_pairs.py
                                 tree-to-table extraction
scripts/add_radicality.py        Grantham and Miyata annotation
scripts/make_itol_*.py           iTOL annotation utilities
data/README.md                   local input contract
results/                         retained figure and generated-output policy
tests/                           tests for reusable Python code
thesis/                          submitted thesis PDF
```

## Setup

Python 3.10 or later and R 4.3 or later are recommended. With
[uv](https://docs.astral.sh/uv/):

```bash
uv sync --extra notebooks --extra test
Rscript scripts/install_r_dependencies.R
```

The equivalent Python installation with pip is:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Running the reusable analysis

Extract pairs of successive amino-acid changes. Every leaf and internal tree
node must be named and have an aligned nucleotide sequence in the FASTA file.

```bash
uv run python scripts/extract_substitution_pairs.py \
  data/raw/tree.newick \
  data/raw/ancestral_sequences.fasta \
  data/processed/substitution_pairs.csv
```

Add parent-to-child radicality values from AAindex records `GRAR740104`
(Grantham) and `MIYT790101` (Miyata):

```bash
uv run python scripts/add_radicality.py \
  data/processed/substitution_pairs.csv \
  data/processed/substitution_pairs_radicality.csv
```

Fit the three thesis model specifications for both metrics:

```bash
Rscript R/fit_radicality_models.R \
  data/processed/substitution_pairs_radicality.csv \
  results/tables/model_summary.csv
```

Run the model script separately for the full, clustered, and senescent tables
from each dataset. It writes one tidy row per response and model specification.

## Tests

```bash
uv run --extra test python -m pytest
```

The tests cover FASTA normalization, codon translation, extraction and distance
calculation on a small named tree, and AAindex radicality annotation. The
study-level numerical results cannot be rerun without the restricted inputs.

## License

The software is released under the [MIT License](LICENSE). The thesis PDF is
copyrighted by Appolinaria Prokopovich and is not covered by the software
license.
