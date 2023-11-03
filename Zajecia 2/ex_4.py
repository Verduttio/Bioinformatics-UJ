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


def count_gc_content(dna):
    dna_length = len(dna)
    letter_gc_counter = 0

    for letter in dna:
        if letter == "G" or letter == "C":
            letter_gc_counter += 1

    return (letter_gc_counter / dna_length) * 100


def dna_array_to_gc_content_array(dnas):
    return list(map(lambda dna: count_gc_content(dna), dnas))


if __name__ == '__main__':
    dnas = parse_fasta_to_dna_array("input.txt")
    dna_name_array = parse_fasta_to_dna_name_array("input.txt")

    gc_array = dna_array_to_gc_content_array(dnas)
    max_gc = max(gc_array)
    index_of_max_gc = gc_array.index(max_gc)
    dna_name_of_max_gc = dna_name_array[index_of_max_gc]

    print(dna_name_of_max_gc)
    print(max_gc)

