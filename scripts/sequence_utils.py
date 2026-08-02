"""Shared FASTA and translation helpers."""

from pathlib import Path

CODON_TO_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C",
    "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R",
    "CGG": "R", "AGA": "R", "AGG": "R", "AGT": "S",
    "AGC": "S", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G", "TAA": "*", "TAG": "*", "TGA": "*",
}


def normalize_sequence_name(header: str) -> str:
    """Match the node-name convention used by the reconstructed trees."""
    return header.split()[0].split("/")[0]


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file into an uppercase name-to-sequence mapping."""
    sequences: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []

    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks).upper().replace("U", "T")
                name = normalize_sequence_name(line[1:])
                if name in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                chunks = []
            elif name is None:
                raise ValueError("FASTA sequence encountered before the first header")
            else:
                chunks.append(line)

    if name is not None:
        sequences[name] = "".join(chunks).upper().replace("U", "T")
    if not sequences:
        raise ValueError(f"No sequences found in {path}")
    return sequences


def translate_codon(codon: str) -> str:
    """Translate a codon, returning X for gaps or ambiguous nucleotides."""
    normalized = codon.upper().replace("U", "T")
    if len(normalized) != 3 or any(base not in "ATCG" for base in normalized):
        return "X"
    return CODON_TO_AA.get(normalized, "X")


def translate_sequence(sequence: str) -> tuple[str, ...]:
    """Translate an in-frame nucleotide sequence into amino-acid symbols."""
    if len(sequence) % 3:
        raise ValueError("Nucleotide sequence length is not divisible by three")
    return tuple(translate_codon(sequence[i : i + 3]) for i in range(0, len(sequence), 3))
