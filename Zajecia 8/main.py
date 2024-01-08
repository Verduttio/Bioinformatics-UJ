from Bio import Phylo


if __name__ == '__main__':
    tree = Phylo.read("phylotree.txt", "newick")
    print(tree)
    Phylo.draw_ascii(tree)
    Phylo.draw(tree)

    terminals = tree.get_terminals()
    for terminal in terminals:
        if terminal.name != "Homo":
            print("Dystans od Homo Sapiens do", terminal.name, "wynosi", tree.distance("Homo", terminal.name))

