from Bio import SeqIO

codon_aa_table = {"UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
                  "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
                  "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
                  "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
                  "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
                  "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
                  "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
                  "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
                  "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
                  "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
                  "UAA": "STOP", "CAA": "Q", "AAA": "K", "GAA": "E",
                  "UAG": "STOP", "CAG": "Q", "AAG": "K", "GAG": "E",
                  "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
                  "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
                  "UGA": "STOP", "CGA": "R", "AGA": "R", "GGA": "G",
                  "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
                  }


def parse_fasta_to_dna_array(filename):
    dnas = []
    for record in SeqIO.parse(filename, "fasta"):
        dnas.append(str(record.seq))
    return dnas


def reverse_complement(section):
    reversed_section = section[::-1]
    complements = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    replaced_chars = [complements.get(char) for char in reversed_section]
    return ''.join(replaced_chars)


def dna_to_rna(dna):
    rna = dna.replace("T", "U")
    return rna


def cumulative_sum(input_list):
    output_list = [input_list[0]]
    for i in range(1, len(input_list)):
        output_list.append(output_list[-1] + input_list[i] + 3)

    return output_list


def get_all_start_codon_indexes(rna):
    indexes = []
    while rna.find("AUG") != -1:
        index = rna.find("AUG")
        indexes.append(index)
        rna = rna[index+3:]

    return cumulative_sum(indexes)


def get_codons(rna):
    codon_len = 3
    return [rna[i: i + codon_len] for i in range(0, len(rna), codon_len)]


def is_orf(proteins):
    return proteins.find("STOP") != -1


def cut_to_stop_indicator(proteins):
    return proteins[:proteins.find("STOP")]


def proteins_code(rna):
    codons = get_codons(rna)
    protein = [codon_aa_table.get(codon, "") for codon in codons]
    protein_str = ''.join(protein)

    return protein_str


def rnas_by_start_codon(rna):
    start_codon_indexes = get_all_start_codon_indexes(rna)
    rnas = []
    for index in start_codon_indexes:
        rnas.append(rna[index:])

    return rnas


def get_orfs(rna):
    orfs = []
    rnas_formatted = rnas_by_start_codon(rna)
    for rna_formatted in rnas_formatted:
        proteins = proteins_code(rna_formatted)
        if is_orf(proteins):
            orfs.append(cut_to_stop_indicator(proteins))

    return orfs


def remove_duplicates(collection):
    return list(set(collection))


def display_formatted_result(collection):
    for element in collection:
        print(element)


if __name__ == '__main__':

    dnas = parse_fasta_to_dna_array("input.txt")
    dna = dnas[0]
    dna_reverse_complement = reverse_complement(dna)

    rna_from_dna = dna_to_rna(dna)
    rna_from_dna_reversed_complement = dna_to_rna(dna_reverse_complement)

    all_orfs = get_orfs(rna_from_dna) + get_orfs(rna_from_dna_reversed_complement)
    all_orfs = remove_duplicates(all_orfs)

    display_formatted_result(all_orfs)
