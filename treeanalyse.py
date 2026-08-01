from ete3 import Tree


t = Tree("/home/ria/HSE/lab/new_three/data/treefile_without_outliers.newick", format=1)
leaves = len(t.get_leaves())
internal = len([n for n in t.traverse() if not n.is_leaf()])
print(f"Листьев: {leaves}, внутренних узлов: {internal}")