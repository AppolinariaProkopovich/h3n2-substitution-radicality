# Exploratory R scripts

These files preserve the interactive analysis history behind the thesis. They
include dataset comparisons, filtering experiments, grouped regressions,
laboratory-site checks, and alternative model specifications.

They are not the supported entry point: several scripts contain personal
absolute paths, rely on objects created earlier in an interactive R session, or
contain abandoned alternatives. Use `R/fit_radicality_models.R` for the clean
downstream model interface. The historical files are retained so that manual
decisions and intermediate checks are auditable.

Known divergence from the thesis: `filter_exclude.r`, which produced the
`..._filtered_concordant_r_site*` tables, lists 13 R-sites, while the thesis
reports 15 (sites 468 and 558 are missing from the script). The tables built by
that script therefore exclude 22 sites rather than the 24 the thesis describes.
The vector is left as it was run; a comment in the script marks the difference.
