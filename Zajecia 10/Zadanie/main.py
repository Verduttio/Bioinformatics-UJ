from Bio import SeqIO
import pylcs
import random


def generate_random_reads(sequence, fragment_length=200, expected_coverage=5):
    number_of_draws = len(sequence) * expected_coverage // fragment_length
    reads = []
    random.seed()
    for _ in range(number_of_draws):
        start_index = random.randint(0, len(sequence) - fragment_length)
        reads.append(sequence[start_index: start_index + fragment_length])
    return reads


def sequences_to_fasta_file(sequences, file_name):
    with open(file_name, "w") as handle:
        for i, sequence in enumerate(sequences):
            handle.write(">read_{}\n".format(i))
            handle.write(sequence + "\n")


def read_sequences_from_fasta_file(file_name):
    sequences = []
    with open(file_name, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    return sequences


def count_greatest_prefix_suffix_len_per_each_one(sequences):
    # output: [[first sequence with second, first with third, ..., first with last],
    #          [second sequence with third, second with fourth, ..., second with last],
    #          [third with fourth, third with fifth, ..., third with last],
    #          ...,
    #          [last-1 with last]
    #         ]

    greatest_prefix_suffix_len_per_each_one = []
    for i in range(len(sequences) - 1):
        greatest_prefix_suffix_len_per_one = prefix_suffix_match(sequences[i], sequences[i + 1:])
        greatest_prefix_suffix_len_per_each_one.append(greatest_prefix_suffix_len_per_one)
    return greatest_prefix_suffix_len_per_each_one


def find_prefix_suffix_length(s1, s2):
    length = 0
    max_length = min(len(s1), len(s2))

    for i in range(1, max_length + 1):
        if s1[:i] == s2[-i:]:
            length = i

    if length == 0:
        s1, s2 = s2, s1
        for i in range(1, max_length + 1):
            if s1[:i] == s2[-i:]:
                length = i

    return length


def prefix_suffix_match(main_str, str_array):
    result = []

    for s in str_array:
        prefix_suffix_length = max(find_prefix_suffix_length(main_str, s), find_prefix_suffix_length(s, main_str))
        result.append(prefix_suffix_length)

    return result


def find_indexes_of_maximum_value_in_2d_list(list_2d):
    max_value = max(map(max, list_2d))
    if max_value == 0:
        return None, None

    for i, row in enumerate(list_2d):
        for j, value in enumerate(row):
            if value == max_value:
                return i, j


def parse_indexes_into_sequence_ids(i, j):
    return i, j+i+1


def find_indexes_of_two_sequences_with_greatest_prefix_suffix_len(sequences):
    greatest_prefix_suffix_len_per_each_one = count_greatest_prefix_suffix_len_per_each_one(sequences)
    i, j = find_indexes_of_maximum_value_in_2d_list(greatest_prefix_suffix_len_per_each_one)
    if i is None:
        return None, None
    return parse_indexes_into_sequence_ids(i, j)


def connect_two_sequences_through_common_sequence(seq1, seq2):
    res = pylcs.lcs_string_idx(seq1, seq2)
    common_sequence = ''.join([seq2[i] for i in res if i != -1])

    seq1_start_index_of_common_seq = seq1.find(common_sequence)

    # seq1 should be beginning sequence
    # [seq1_common_seq2]
    if seq1_start_index_of_common_seq == 0:
        seq1, seq2 = seq2, seq1

    return seq1 + seq2[len(common_sequence):]


def assembly_sequence_from_reads(sequences):
    while len(sequences) > 1:
        print("Log: {} sequences left".format(len(sequences)))
        seq1_index, seq2_index = find_indexes_of_two_sequences_with_greatest_prefix_suffix_len(sequences)
        if seq1_index is None:
            break
        seq1, seq2 = sequences[seq1_index], sequences[seq2_index]
        connected_sequence = connect_two_sequences_through_common_sequence(seq1, seq2)
        sequences.append(connected_sequence)
        sequences.remove(seq1)
        sequences.remove(seq2)
    return sequences


def assembly_sequence_from_file(file_name):
    sequences = read_sequences_from_fasta_file(file_name)
    return assembly_sequence_from_reads(sequences)


if __name__ == '__main__':
    y_chromosome = read_sequences_from_fasta_file("seq.fasta")[0]
    random_reads = generate_random_reads(y_chromosome, 200, 5)
    sequences_to_fasta_file(random_reads, "reads.fasta")

    result_sequences = assembly_sequence_from_file("reads.fasta")

    print("Liczba kontigów: {}".format(len(result_sequences)))
    print("Długości kontigów: {}".format([len(x) for x in result_sequences]))
    print("Sumaryczna długość kontigów: {}".format(sum([len(x) for x in result_sequences])))
    print("Długość chromosomu Y: {}".format(len(y_chromosome)))

    if len(result_sequences) == 1:
        lcs_len = pylcs.lcs(y_chromosome, result_sequences[0])
        print("Długość najdłuższego podciągu: {}".format(lcs_len))
        print("Podobieństwo: {}".format(lcs_len / len(y_chromosome)))

    # Zwrócona sekwencja nie jest identyczna jak bazowa.
    # Jednak jest do niej bardzo podobna.
    # Podobieństwo było różne i w zależności od rozlosowania zbliżało się do 98%.
