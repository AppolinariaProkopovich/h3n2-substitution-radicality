from pathlib import Path

import pytest

from scripts.sequence_utils import read_fasta, translate_codon, translate_sequence


def test_read_fasta_normalizes_names_and_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(">node1/comment metadata\naug\n>node2\nGAC\n", encoding="utf-8")

    assert read_fasta(fasta) == {"node1": "ATG", "node2": "GAC"}


def test_translation_handles_standard_and_ambiguous_codons() -> None:
    assert translate_codon("GAA") == "E"
    assert translate_codon("NNN") == "X"
    assert translate_sequence("ATGGAC") == ("M", "D")


def test_translation_rejects_out_of_frame_sequence() -> None:
    with pytest.raises(ValueError, match="divisible by three"):
        translate_sequence("ATGG")
