#!/usr/bin/env python3
"""Add Grantham and Miyata radicality to a substitution-pair table."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import pandas as pd

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def annotate_radicality(data: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with parent-to-child Grantham and Miyata distances."""
    required = {"anc_aa", "der_aa"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")

    invalid = ~data["anc_aa"].isin(STANDARD_AMINO_ACIDS) | ~data["der_aa"].isin(
        STANDARD_AMINO_ACIDS
    )
    if invalid.any():
        row_numbers = ", ".join(str(index) for index in data.index[invalid][:10])
        raise ValueError(f"Non-standard amino acid in input row(s): {row_numbers}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        from aaindex import aaindex2
        grantham = aaindex2["GRAR740104"].matrix
        miyata = aaindex2["MIYT790101"].matrix
    output = data.copy()
    if output.empty:
        output["gr_rad_anc"] = pd.Series(dtype=float)
        output["miyata_rad"] = pd.Series(dtype=float)
        return output

    pairs = zip(output["anc_aa"], output["der_aa"], strict=True)
    values = [(grantham[a][b], miyata[a][b]) for a, b in pairs]
    output[["gr_rad_anc", "miyata_rad"]] = values
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input substitution-pair CSV")
    parser.add_argument("output", type=Path, help="Output annotated CSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = annotate_radicality(pd.read_csv(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(f"Wrote {len(output)} annotated substitutions to {args.output}")


if __name__ == "__main__":
    main()
