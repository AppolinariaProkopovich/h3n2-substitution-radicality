"""
Генерация файла аннотации iTOL (DATASET_COLORSTRIP)
для одного аминокислотного сайта по ДНК-фасте.

Использование:
    python make_itol_site_annotation.py <fasta> <site_number> [output]

Пример:
    python make_itol_site_annotation.py ancestral_sequences.fasta 78 itol_site78.txt
"""

import sys

# ── Настройки ──────────────────────────────────────────────
SITE = 78  # аминокислотный сайт (1-based), можно переопределить аргументом

# Цвета для 20 стандартных аминокислот (+ X для неизвестных)
AA_COLORS = {
    "A": "#f44336", "R": "#e91e63", "N": "#9c27b0", "D": "#673ab7",
    "C": "#3f51b5", "Q": "#2196f3", "E": "#03a9f4", "G": "#00bcd4",
    "H": "#009688", "I": "#4caf50", "L": "#8bc34a", "K": "#cddc39",
    "M": "#ffeb3b", "F": "#ffc107", "P": "#ff9800", "S": "#ff5722",
    "T": "#795548", "W": "#607d8b", "Y": "#9e9e9e", "V": "#827717",
    "X": "#bdbdbd", "*": "#000000",
}

CODON2AA = {
    "TTT":"F","TTC":"F",
    "TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I",
    "ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y",
    "CAT":"H","CAC":"H",
    "CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N",
    "AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D",
    "GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C",
    "TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "AGT":"S","AGC":"S",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    "TAA":"*","TAG":"*","TGA":"*",
}


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split("/")[0]  # убираем суффикс после "/"
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def get_aa_at_site(nuc_seq, site_1based):
    """Извлекает аминокислоту в указанном сайте (1-based) из нуклеотидной последовательности."""
    start = (site_1based - 1) * 3
    codon = nuc_seq[start:start + 3].upper().replace("U", "T")
    if len(codon) < 3 or any(x not in "ATCG" for x in codon):
        return "X"
    return CODON2AA.get(codon, "X")


def make_itol_colorstrip(seqs, site, output_path):
    """Генерирует файл DATASET_COLORSTRIP для iTOL."""
    # Собираем данные: имя → аминокислота
    data = {}
    for name, seq in seqs.items():
        aa = get_aa_at_site(seq, site)
        data[name] = aa

    # Какие аминокислоты реально встречаются
    unique_aas = sorted(set(data.values()))
    print(f"Сайт {site}: найдено {len(data)} последовательностей")
    print(f"Уникальные аминокислоты: {', '.join(unique_aas)}")
    for aa in unique_aas:
        count = sum(1 for v in data.values() if v == aa)
        print(f"  {aa}: {count}")

    # Формируем легенду только для встречающихся аминокислот
    legend_labels = unique_aas
    legend_colors = [AA_COLORS.get(aa, "#bdbdbd") for aa in unique_aas]

    with open(output_path, "w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\tSite_{site}\n")
        f.write("COLOR\t#000000\n")
        f.write("STRIP_WIDTH\t30\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("MARGIN\t5\n\n")

        # Легенда
        tab = "\t"
        shapes_str = tab.join(["1"] * len(legend_labels))
        colors_str = tab.join(legend_colors)
        labels_str = tab.join(legend_labels)
        f.write(f"LEGEND_TITLE\tAmino acids at site {site}\n")
        f.write(f"LEGEND_SHAPES\t{shapes_str}\n")
        f.write(f"LEGEND_COLORS\t{colors_str}\n")
        f.write(f"LEGEND_LABELS\t{labels_str}\n\n")

        # Данные
        f.write("DATA\n")
        for name, aa in sorted(data.items()):
            color = AA_COLORS.get(aa, "#bdbdbd")
            f.write(f"{name}\t{color}\t{aa}\n")

    print(f"\nФайл сохранён: {output_path}")
    print("Перетащите его на дерево в iTOL для отображения.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Использование: python make_itol_site_annotation.py <fasta> [site] [output]")
        print("  fasta  — путь к выровненной ДНК-фасте (с предковыми последовательностями)")
        print("  site   — номер аминокислотного сайта, по умолчанию 78")
        print("  output — имя выходного файла, по умолчанию itol_site_<N>.txt")
        sys.exit(1)

    fasta_path = sys.argv[1]
    site = int(sys.argv[2]) if len(sys.argv) > 2 else SITE
    output = sys.argv[3] if len(sys.argv) > 3 else f"itol_site_{site}.txt"

    seqs = read_fasta(fasta_path)
    print(f"Прочитано последовательностей: {len(seqs)}")

    make_itol_colorstrip(seqs, site, output)
