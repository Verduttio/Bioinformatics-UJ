from Bio import SeqIO


def parse_fasta_to_dna_array(filename):
    dnas = []
    for record in SeqIO.parse(filename, "fasta"):
        dnas.append(str(record.seq))
    return dnas


def count_transitions(dna1, dna2):
    counter = 0
    for symbol1, symbol2 in zip(dna1, dna2):
        if is_transition(symbol1, symbol2):
            counter += 1
    return counter


def count_transversions(dna1, dna2):
    counter = 0
    for symbol1, symbol2 in zip(dna1, dna2):
        if is_transversion(symbol1, symbol2):
            counter += 1
    return counter


def is_transition(symbol1, symbol2):
    if (symbol1 == "A" and symbol2 == "G") or (symbol1 == "G" and symbol2 == "A"):
        return True
    elif (symbol1 == "C" and symbol2 == "T") or (symbol1 == "T" and symbol2 == "C"):
        return True
    else:
        return False


def count_ratio(dna1, dna2):
    return count_transitions(dna1, dna2)/count_transversions(dna1, dna2)


def is_transversion(symbol1, symbol2):
    return (symbol1 != symbol2) and (not is_transition(symbol1, symbol2))


if __name__ == '__main__':
    dnas = parse_fasta_to_dna_array("input.txt")
    dna1 = dnas[0]
    dna2 = dnas[1]
    print("Ratio: ", count_ratio(dna1, dna2))
