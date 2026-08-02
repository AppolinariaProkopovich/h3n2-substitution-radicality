#!/usr/bin/env python3
"""Extract consecutive amino-acid substitutions from an ancestral tree."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

try:
    from scripts.sequence_utils import read_fasta, translate_sequence
except ModuleNotFoundError:  # Direct execution from the scripts directory.
    from sequence_utils import read_fasta, translate_sequence

OUTPUT_COLUMNS = [
    "site",
    "ancestor_node",
    "grandanc_codon",
    "grandanc_aa",
    "anc_codon",
    "anc_aa",
    "der_codon",
    "der_aa",
    "distance",
    "distance_grand",
]


def validate_inputs(tree: Tree, sequences: dict[str, str]) -> int:
    """Validate tree/FASTA agreement and return the number of codon sites."""
    unnamed = [clade for clade in tree.find_clades() if not clade.name]
    if unnamed:
        raise ValueError(
            f"Tree contains {len(unnamed)} unnamed node(s); ancestral nodes must be named"
        )

    tree_names = {clade.name for clade in tree.find_clades()}
    missing = sorted(tree_names - sequences.keys())
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"FASTA lacks sequences for {len(missing)} tree node(s): {preview}")

    lengths = {len(sequences[name]) for name in tree_names}
    if len(lengths) != 1:
        raise ValueError(f"Aligned tree sequences have inconsistent lengths: {sorted(lengths)}")
    sequence_length = lengths.pop()
    if sequence_length % 3:
        raise ValueError(f"Alignment length {sequence_length} is not divisible by three")
    return sequence_length // 3


def extract_substitution_pairs(
    tree: Tree, sequences: dict[str, str]
) -> list[dict[str, Any]]:
    """Find pairs of successive nonsynonymous changes along root-to-tip paths.

    Distance is measured from the descendant node of the first substitution to
    the midpoint of the branch carrying the second substitution, matching the
    final extraction notebook used for the thesis.
    """
    n_sites = validate_inputs(tree, sequences)
    tree_names = {clade.name for clade in tree.find_clades()}
    amino_acids = {
        name: translate_sequence(sequences[name]) for name in tree_names
    }
    rows: list[dict[str, Any]] = []
    parent_map = {
        child: parent for parent in tree.find_clades() for child in parent.clades
    }

    def traverse(clade: Clade, last_events: dict[int, dict[str, Any]]) -> None:
        parent = parent_map.get(clade)

        if parent is not None:
            parent_name = parent.name
            node_name = clade.name
            parent_aa = amino_acids[parent_name]
            node_aa = amino_acids[node_name]

            for site_index in range(n_sites):
                before = parent_aa[site_index]
                after = node_aa[site_index]
                if before in {"X", "*"} or after in {"X", "*"} or before == after:
                    continue

                previous = last_events.get(site_index)
                if previous is not None:
                    start = site_index * 3
                    branch_length = clade.branch_length or 0.0
                    rows.append(
                        {
                            "site": site_index + 1,
                            "ancestor_node": previous["anc_node"].name,
                            "grandanc_codon": sequences[previous["grandanc_node"].name][start : start + 3],
                            "grandanc_aa": previous["grandanc_aa"],
                            "anc_codon": sequences[previous["anc_node"].name][start : start + 3],
                            "anc_aa": previous["anc_aa"],
                            "der_codon": sequences[node_name][start : start + 3],
                            "der_aa": after,
                            "distance": tree.distance(previous["anc_node"], clade)
                            - branch_length / 2,
                            "distance_grand": tree.distance(previous["grandanc_node"], clade)
                            - branch_length / 2,
                        }
                    )

                last_events[site_index] = {
                    "grandanc_aa": before,
                    "grandanc_node": parent,
                    "anc_aa": after,
                    "anc_node": clade,
                }

        for child in clade.clades:
            traverse(child, last_events.copy())

    traverse(tree.root, {})
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path, help="Named ancestral tree in Newick format")
    parser.add_argument("fasta", type=Path, help="Aligned nucleotide sequences for every tree node")
    parser.add_argument("output", type=Path, help="Output CSV path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tree = Phylo.read(args.tree, "newick")
    sequences = read_fasta(args.fasta)
    rows = extract_substitution_pairs(tree, sequences)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(args.output, index=False)
    print(f"Wrote {len(rows)} substitution pairs to {args.output}")


if __name__ == "__main__":
    main()
