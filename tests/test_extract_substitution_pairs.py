from io import StringIO

import pytest
from Bio import Phylo

from scripts.extract_substitution_pairs import extract_substitution_pairs


def test_extracts_successive_changes_and_midpoint_distance() -> None:
    tree = Phylo.read(StringIO("(((leaf:1)node2:1)node1:1)root;"), "newick")
    sequences = {
        "root": "AAA",   # K
        "node1": "GAA",  # E: first substitution
        "node2": "GAC",  # D: second substitution
        "leaf": "GAC",
    }

    rows = extract_substitution_pairs(tree, sequences)

    assert len(rows) == 1
    assert rows[0] == {
        "site": 1,
        "ancestor_node": "node1",
        "grandanc_codon": "AAA",
        "grandanc_aa": "K",
        "anc_codon": "GAA",
        "anc_aa": "E",
        "der_codon": "GAC",
        "der_aa": "D",
        "distance": pytest.approx(0.5),
        "distance_grand": pytest.approx(1.5),
    }


def test_requires_sequence_for_every_named_node() -> None:
    tree = Phylo.read(StringIO("(leaf:1)root;"), "newick")

    with pytest.raises(ValueError, match="FASTA lacks sequences"):
        extract_substitution_pairs(tree, {"root": "AAA"})
