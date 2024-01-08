import numpy as np
from Bio.Align import substitution_matrices
from Bio import Entrez, SeqIO


def which_is_max(first, second, third):
    values = [first, second, third]
    max_value = max(values)
    max_indices = [i + 1 for i, value in enumerate(values) if value == max_value]
    return max_indices


def coordinates_of_max_value(which_max_value, current_i, current_j):
    if which_max_value == 1:
        return current_i - 1, current_j - 1
    elif which_max_value == 2:
        return current_i - 1, current_j
    else:
        return current_i, current_j - 1


def algorithm_with_blosum62(gap, s1, s2):
    arr = np.zeros((len(s2) + 1, len(s1) + 1))
    my_parents_arr = np.empty((len(s2) + 1, len(s1) + 1), dtype=object)

    # fill first row and column
    for i in range(1, len(arr[0])):
        arr[0][i] = i * gap
        my_parents_arr[0][i] = [(0, i - 1)]
    for i in range(1, len(arr)):
        arr[i][0] = i * gap
        my_parents_arr[i][0] = [(i - 1, 0)]

    for i in range(1, len(arr)):
        for j in range(1, len(arr[0])):
            match_or_not_value = get_value_from_blosum62(s1[j-1], s2[i-1])
            first = arr[i - 1][j - 1] + match_or_not_value
            second = arr[i - 1][j] + gap
            third = arr[i][j - 1] + gap

            arr[i][j] = max(first, second, third)
            which_max_value = which_is_max(first, second, third)

            my_parents_arr[i][j] = []
            for max_value in which_max_value:
                my_parents_arr[i][j].append(coordinates_of_max_value(max_value, i, j))

    return arr, my_parents_arr


def read_paths_coordinates(matrix):
    rows = len(matrix)
    cols = len(matrix[0])

    def get_paths(row, col):
        if row == 0 and col == 0:
            return [[]]

        paths = []
        for parent in matrix[row][col]:
            r, c = parent
            for path in get_paths(r, c):
                paths.append(path + [(row, col)])

        return paths

    return get_paths(rows - 1, cols - 1)


def get_optimal_match(s1, s2, path):
    optimal_match_s1 = ""
    optimal_match_s2 = ""
    last_i = 0
    last_j = 0
    for coordinates in path:
        coordinate_i = coordinates[0]
        coordinate_j = coordinates[1]

        if last_i != coordinate_i:
            optimal_match_s2 += s2[coordinate_i-1] + " "
            last_i = coordinate_i
        else:
            optimal_match_s2 += "- "

        if last_j != coordinate_j:
            optimal_match_s1 += s1[coordinate_j-1] + " "
            last_j = coordinate_j
        else:
            optimal_match_s1 += "- "

    return optimal_match_s1 + "\n" + optimal_match_s2


def get_optimal_matches(s1, s2, paths):
    optimal_matches = []
    for path in paths:
        optimal_match = get_optimal_match(s1, s2, path)
        optimal_matches.append(optimal_match)
    return optimal_matches


def swap_sequences_first_should_be_longer(s1, s2):
    if len(s1) < len(s2):
        return s2, s1
    else:
        return s1, s2


def get_value_from_blosum62(nucleotide1, nucleotide2):
    blosum62 = substitution_matrices.load("BLOSUM62")
    return blosum62[nucleotide1][nucleotide2]


def optimal_matches_algorithm_blosum62(seq1, seq2, gap):
    seq1, seq2 = swap_sequences_first_should_be_longer(seq1, seq2)

    output_array, my_parents = algorithm_with_blosum62(gap, seq1, seq2)
    paths_coordinates = read_paths_coordinates(my_parents)
    optimal_matches = get_optimal_matches(seq1, seq2, paths_coordinates)
    return optimal_matches


def print_optimal_matches(optimal_matches):
    for match in optimal_matches:
        print(match)
        print()


def optimal_matches_algorithm(seq1, seq2, match, mismatch, gap):
    seq1, seq2 = swap_sequences_first_should_be_longer(seq1, seq2)

    output_array, my_parents = algorithm(match, mismatch, gap, seq1, seq2)
    paths_coordinates = read_paths_coordinates(my_parents)
    optimal_matches = get_optimal_matches(seq1, seq2, paths_coordinates)
    return optimal_matches


def compare_two_sequences(seq1, seq2, gap):
    optimal_matches_arr = optimal_matches_algorithm(seq1, seq2, match=1, mismatch=-1, gap=gap)
    print_optimal_matches(optimal_matches_arr)


### Optional ###
def check_match(seq1, seq2, pos1, pos2):
    return seq1[pos1] == seq2[pos2]


def algorithm(match, mismatch, gap, s1, s2):
    arr = np.zeros((len(s2) + 1, len(s1) + 1))
    my_parents_arr = np.empty((len(s2) + 1, len(s1) + 1), dtype=object)

    # fill first row and column
    for i in range(1, len(arr[0])):
        arr[0][i] = i * gap
        my_parents_arr[0][i] = [(0, i - 1)]
    for i in range(1, len(arr)):
        arr[i][0] = i * gap
        my_parents_arr[i][0] = [(i - 1, 0)]

    for i in range(1, len(arr)):
        for j in range(1, len(arr[0])):
            match_or_not_value = match if check_match(s1, s2, pos1=j-1, pos2=i-1) else mismatch
            first = arr[i - 1][j - 1] + match_or_not_value
            second = arr[i - 1][j] + gap
            third = arr[i][j - 1] + gap

            arr[i][j] = max(first, second, third)
            which_max_value = which_is_max(first, second, third)

            my_parents_arr[i][j] = []
            for max_value in which_max_value:
                my_parents_arr[i][j].append(coordinates_of_max_value(max_value, i, j))

    return arr, my_parents_arr


def get_arrays_with_only_path_elements(arr, paths):
    arrays = []
    for path in paths:
        array_with_path = get_array_with_only_path_elements(arr, path)
        arrays.append(array_with_path)
    return arrays


def get_array_with_only_path_elements(arr, path):
    arr_copy = np.copy(arr)
    for i in range(len(arr)):
        for j in range(len(arr[0])):
            if (i, j) not in path:
                arr_copy[i][j] = np.nan
    return arr_copy


def get_gi_40886941_sequence():
    Entrez.email = "bszwaja20@gmail.com"
    handle = Entrez.efetch(db="protein", id="40886941", rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank").seq


def get_gi_34849618_sequence():
    Entrez.email = "bszwaja20@gmail.com"
    handle = Entrez.efetch(db="protein", id="34849618", rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank").seq


if __name__ == '__main__':
    seq_40886941 = get_gi_40886941_sequence()
    seq_34849618 = get_gi_34849618_sequence()

    print("Optymalne dopasowanie dla sekwencji hemoglobiny beta człowieka i szczura:")
    compare_two_sequences(seq_40886941, seq_34849618, gap=-7)
