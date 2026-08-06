# Research notebooks

These notebooks document the original development of the tree traversal,
distance comparison, and laboratory-passage categorization steps. Their cell
order, outputs, and study-specific paths are preserved as research provenance.

For new data, use `scripts/extract_substitution_pairs.py` rather than running
the extraction notebooks cell by cell.

`passage_phenotype_categories.ipynb` turns GISAID passage strings into the
binary clinical/laboratory phenotype used by the EpiStat association analysis.
It only prepares that phenotype table; the association test itself was run
outside this repository. It needs restricted strain metadata and cannot run
without separately authorised inputs.
