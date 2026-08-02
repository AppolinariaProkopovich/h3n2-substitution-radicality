#!/usr/bin/env python3
"""Report leaf and internal-node counts for a Newick tree."""

import argparse
from pathlib import Path

from Bio import Phylo


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path, help="Tree in Newick format")
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    leaves = tree.count_terminals()
    total = sum(1 for _ in tree.find_clades())
    print(f"Leaves: {leaves}")
    print(f"Internal nodes: {total - leaves}")


if __name__ == "__main__":
    main()
