"""
Генерация файла аннотации iTOL (TREE_COLORS)
для подсветки рёбер, на которых произошла аминокислотная замена в заданном сайте.

Также генерирует файл COLLAPSE для схлопывания клад без замен.

Использование:
    python make_itol_branch_colors.py <fasta> <newick> [site] [output_prefix]

Пример:
    python make_itol_branch_colors.py ancestral_sequences.fasta treefile.newick 78
"""

import sys
from pathlib import Path

# ── Настройки ──────────────────────────────────────────────
SITE =78
COLOR_SUBST = "#e41a1c"   # красный — ветка с заменой
COLOR_DEFAULT = "#cccccc"  # серый — без замены
WIDTH_SUBST = 4            # толщина ветки с заменой
WIDTH_DEFAULT = 1          # толщина обычной ветки

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
                name = line[1:].split("/")[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def get_aa_at_site(nuc_seq, site_1based):
    start = (site_1based - 1) * 3
    codon = nuc_seq[start:start + 3].upper().replace("U", "T")
    if len(codon) < 3 or any(x not in "ATCG" for x in codon):
        return "X"
    return CODON2AA.get(codon, "X")


def main():
    if len(sys.argv) < 3:
        print("Использование: python make_itol_branch_colors.py <fasta> <newick> [site] [output_prefix]")
        sys.exit(1)

    fasta_path = sys.argv[1]
    tree_path = sys.argv[2]
    site = int(sys.argv[3]) if len(sys.argv) > 3 else SITE
    prefix = sys.argv[4] if len(sys.argv) > 4 else f"itol_site_{site}"

    # ── Читаем данные ──
    from ete3 import Tree
    seqs = read_fasta(fasta_path)
    print(f"Прочитано последовательностей: {len(seqs)}")

    # Получаем аминокислоту для каждого имени
    aa_map = {}
    for name, seq in seqs.items():
        aa_map[name] = get_aa_at_site(seq, site)

    t = Tree(tree_path, format=1)

    # ── Находим рёбра с заменами ──
    subst_edges = []   # (node_name, parent_aa, child_aa)
    no_subst_nodes = []

    for node in t.traverse("preorder"):
        if node.is_root():
            continue
        parent = node.up
        if not parent.name or not node.name:
            continue
        p_aa = aa_map.get(parent.name)
        n_aa = aa_map.get(node.name)
        if p_aa is None or n_aa is None:
            continue
        if p_aa in ("X", "*") or n_aa in ("X", "*"):
            continue

        if p_aa != n_aa:
            subst_edges.append((node.name, p_aa, n_aa))
        else:
            no_subst_nodes.append(node.name)

    print(f"\nСайт {site}:")
    print(f"  Рёбер с заменой: {len(subst_edges)}")
    print(f"  Рёбер без замены: {len(no_subst_nodes)}")
    print(f"\nЗамены:")
    for name, p_aa, n_aa in subst_edges:
        print(f"  {name}: {p_aa} -> {n_aa}")

    # ── Файл 1: TREE_COLORS ──
    colors_file = f"{prefix}_branches.txt"
    with open(colors_file, "w") as f:
        f.write("TREE_COLORS\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"# Подсветка замен в аминокислотном сайте {site}\n")
        f.write(f"# Красные ветви = замена, серые = без замены\n\n")
        f.write("DATA\n")

        # Серые ветки (без замены)
        for name in no_subst_nodes:
            f.write(f"{name}\tbranch\t{COLOR_DEFAULT}\tnormal\t{WIDTH_DEFAULT}\n")

        # Красные ветки (с заменой)
        for name, p_aa, n_aa in subst_edges:
            f.write(f"{name}\tbranch\t{COLOR_SUBST}\tnormal\t{WIDTH_SUBST}\n")

        # Подписи на ветках с заменами
        for name, p_aa, n_aa in subst_edges:
            f.write(f"{name}\tlabel\t{COLOR_SUBST}\tbold\t2\n")

    print(f"\nСохранён: {colors_file}")

    # ── Файл 2: DATASET_TEXT — подписи замен на ветках ──
    labels_file = f"{prefix}_labels.txt"
    with open(labels_file, "w") as f:
        f.write("DATASET_TEXT\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\tSubstitutions at site {site}\n")
        f.write("COLOR\t#e41a1c\n")
        f.write("ALIGN_TO_TREE\t1\n")
        f.write("SIZE_FACTOR\t0.8\n\n")
        f.write("DATA\n")
        for name, p_aa, n_aa in subst_edges:
            f.write(f"{name}\t{p_aa}->{n_aa}\t-1\t{COLOR_SUBST}\tbold\t1\n")

    print(f"Сохранён: {labels_file}")

    # ── Файл 3: COLLAPSE — список клад для схлопывания ──
    # Один обход postorder (снизу вверх): O(n)
    # Для каждого узла считаем: есть ли замена в поддереве и сколько листьев.
    subst_node_names = set(name for name, _, _ in subst_edges)

    has_subst = {}   # node id -> bool
    leaf_count = {}  # node id -> int

    for node in t.traverse("postorder"):
        nid = id(node)
        if node.is_leaf():
            has_subst[nid] = node.name in subst_node_names
            leaf_count[nid] = 1
        else:
            child_has = any(has_subst[id(c)] for c in node.children)
            has_subst[nid] = child_has or (node.name in subst_node_names)
            leaf_count[nid] = sum(leaf_count[id(c)] for c in node.children)

    # Обход preorder: берём самые крупные клады без замен.
    # Если узел без замен — добавляем его и НЕ спускаемся в потомков
    # (так автоматически убираются вложенные).
    final_collapse = []
    stack = [t]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            continue
        nid = id(node)
        if not has_subst[nid] and not node.is_root() and leaf_count[nid] >= 3:
            # Вся клада без замен — collapse, не спускаемся глубже
            if node.name:
                final_collapse.append((node.name, leaf_count[nid]))
        else:
            # Есть замены внутри — спускаемся дальше
            stack.extend(node.children)

    collapse_file = f"{prefix}_collapse.txt"
    with open(collapse_file, "w") as f:
        f.write("COLLAPSE\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATA\n")
        for name, n_leaves in final_collapse:
            f.write(f"{name}\n")

    print(f"Сохранён: {collapse_file}")
    print(f"  Клад для схлопывания: {len(final_collapse)}")
    for name, n_leaves in sorted(final_collapse, key=lambda x: -x[1])[:10]:
        print(f"    {name} ({n_leaves} листьев)")

    print(f"\n=== Инструкция ===")
    print(f"1. Загрузите дерево в iTOL")
    print(f"2. Перетащите {colors_file} — покрасит ветви")
    print(f"3. Перетащите {labels_file} — добавит подписи замен")
    print(f"4. Перетащите {collapse_file} — схлопнет клады без замен")


if __name__ == "__main__":
    main()
