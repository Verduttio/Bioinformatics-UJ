from Bio import SeqIO


def parse_fasta_to_dna_array(filename):
    dnas = []
    for record in SeqIO.parse(filename, "fasta"):
        dnas.append(str(record.seq))
    return dnas


def parse_fasta_to_dna_name_array(filename):
    names = []
    for record in SeqIO.parse(filename, "fasta"):
        names.append(str(record.name))
    return names


def get_3_letter_suffix(string):
    return string[len(string)-3:]


def get_3_letter_prefix(string):
    return string[:3]


def get_adjacency_list(dna_structures, dna_names):
    adjacency_list = []
    for dna in dna_structures:
        dna_suffix = get_3_letter_suffix(dna)
        dna_name = dna_names[dna_structures.index(dna)]
        for dna_neighbour_candidate in dna_structures:
            if dna_neighbour_candidate != dna:
                dna_neighbour_candidate_prefix = get_3_letter_prefix(dna_neighbour_candidate)
                if dna_suffix == dna_neighbour_candidate_prefix:
                    dna_neighbour_candidate_name = dna_names[dna_structures.index(dna_neighbour_candidate)]
                    adjacency_list.append((dna_name, dna_neighbour_candidate_name))
    return adjacency_list


def print_formatted_adjacency_list(adjacency_list):
    for name1, name2 in adjacency_list:
        print(name1, name2)


if __name__ == '__main__':
    dnas = parse_fasta_to_dna_array("input.txt")
    dna_name_array = parse_fasta_to_dna_name_array("input.txt")

    adjacency_list = get_adjacency_list(dnas, dna_name_array)
    print_formatted_adjacency_list(adjacency_list)

