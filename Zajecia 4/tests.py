import numpy as np
import unittest
from Bio.Align import substitution_matrices
from main import which_is_max
from main import algorithm
from main import algorithm_with_blosum62
from main import read_paths_coordinates
from main import get_arrays_with_only_path_elements
from main import get_optimal_matches
from main import optimal_matches_algorithm, optimal_matches_algorithm_blosum62
from main import get_value_from_blosum62
from main import get_gi_40886941_sequence
from main import get_gi_34849618_sequence
from main import print_optimal_matches


class TestSum(unittest.TestCase):
    def test_which_is_max(self):
        self.assertEqual(which_is_max(1, 2, 3), [3])
        self.assertEqual(which_is_max(1, 2, 1), [2])
        self.assertEqual(which_is_max(1, 1, 3), [3])
        self.assertEqual(which_is_max(1, 2, 2), [2, 3])
        self.assertEqual(which_is_max(3, 2, 1), [1])
        self.assertEqual(which_is_max(1, 3, 2), [2])
        self.assertEqual(which_is_max(1, 2, 3), [3])
        self.assertEqual(which_is_max(1, 1, 1), [1, 2, 3])
        self.assertEqual(which_is_max(1, 1, 2), [3])

    def test_output_array(self):
        s1 = "TGCTCGTA"
        s2 = "TTCATA"

        match = 5
        mismatch = -2
        gap = -6

        output_array, my_parents = algorithm(match, mismatch, gap, s1, s2)
        output_array = output_array.tolist()
        assert_output_array = [[0., -6., -12., -18., -24., -30., -36., -42., -48.],
                               [-6., 5., -1., -7., -13., -19., -25., -31., -37.],
                               [-12., -1., 3., -3., -2., -8., -14., -20., -26.],
                               [-18., -7., -3., 8., 2., 3., -3., -9., -15.],
                               [-24., -13., -9., 2., 6., 0., 1., -5., -4.],
                               [-30., -19., -15., -4., 7., 4., -2., 6., 0.],
                               [-36., -25., -21., -10., 1., 5., 2., 0., 11.]]
        self.assertEqual(assert_output_array, output_array)

    def test_one_path_coordinates(self):
        my_parents = [
            [None, list([(0, 0)]), list([(0, 1)]), list([(0, 2)]), list([(0, 3)]), list([(0, 4)]), list([(0, 5)]),
             list([(0, 6)]), list([(0, 7)])],
            [list([(0, 0)]), list([(0, 0)]), list([(1, 1)]), list([(1, 2)]),
             list([(0, 3), (1, 3)]), list([(1, 4)]), list([(1, 5)]),
             list([(0, 6), (1, 6)]), list([(1, 7)])],
            [list([(1, 0)]), list([(1, 0), (1, 1)]), list([(1, 1)]), list([(1, 2), (2, 2)]),
             list([(1, 3)]), list([(2, 4)]), list([(2, 5)]), list([(1, 6), (2, 6)]),
             list([(2, 7)])],
            [list([(2, 0)]), list([(2, 1)]), list([(2, 1), (2, 2)]), list([(2, 2)]),
             list([(3, 3)]), list([(2, 4)]), list([(3, 5)]), list([(3, 6)]),
             list([(3, 7)])],
            [list([(3, 0)]), list([(3, 1)]), list([(3, 1), (3, 2)]), list([(3, 3)]),
             list([(3, 3)]), list([(3, 4), (4, 4)]), list([(3, 5)]),
             list([(3, 6), (4, 6)]), list([(3, 7)])],
            [list([(4, 0)]), list([(4, 0), (4, 1)]), list([(4, 1), (4, 2)]), list([(4, 3)]),
             list([(4, 3)]), list([(4, 4)]), list([(4, 5), (5, 5)]), list([(4, 6)]),
             list([(5, 7)])],
            [list([(5, 0)]), list([(5, 1)]), list([(5, 1), (5, 2)]), list([(5, 3)]),
             list([(5, 4)]), list([(5, 4)]), list([(5, 5)]), list([(5, 7)]),
             list([(5, 7)])]]

        paths_coordinates = read_paths_coordinates(my_parents)
        # print(paths_coordinates)
        assert_paths_coordinates = [[(1, 1), (1, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 7), (6, 8)]]
        self.assertEqual(assert_paths_coordinates, paths_coordinates)

    def test_two_paths_coordinates(self):
        my_parents = [[None, list([(0, 0)]), list([(0, 1)]), list([(0, 2)]), list([(0, 3)])],
                      [list([(0, 0)]), list([(0, 0)]), list([(0, 1), (1, 1)]), list([(1, 2)]),
                       list([(1, 3)])],
                      [list([(1, 0)]),
                       list([(1, 1)]),
                       list([(1, 1)]),
                       list([(1, 2)]),
                       list([(2, 3)])],
                      [list([(2, 0)]), list([(2, 1)]), list([(2, 1), (2, 2)]), list([(2, 2)]),
                       list([(2, 3)])]]

        paths_coordinates = read_paths_coordinates(my_parents)
        # print(paths_coordinates)
        assert_paths_coordinates = [[(0, 1), (1, 2), (2, 3), (3, 4)], [(1, 1), (1, 2), (2, 3), (3, 4)]]
        self.assertEqual(assert_paths_coordinates, paths_coordinates)

    def test_my_parents(self):
        s1 = "TGCTCGTA"
        s2 = "TTCATA"

        match = 5
        mismatch = -2
        gap = -6

        output_array, my_parents = algorithm(match, mismatch, gap, s1, s2)
        my_parents = my_parents.tolist()
        assert_my_parents = [
            [None, list([(0, 0)]), list([(0, 1)]), list([(0, 2)]), list([(0, 3)]), list([(0, 4)]), list([(0, 5)]),
             list([(0, 6)]), list([(0, 7)])],
            [list([(0, 0)]), list([(0, 0)]), list([(1, 1)]), list([(1, 2)]),
             list([(0, 3), (1, 3)]), list([(1, 4)]), list([(1, 5)]),
             list([(0, 6), (1, 6)]), list([(1, 7)])],
            [list([(1, 0)]), list([(1, 0), (1, 1)]), list([(1, 1)]), list([(1, 2), (2, 2)]),
             list([(1, 3)]), list([(2, 4)]), list([(2, 5)]), list([(1, 6), (2, 6)]),
             list([(2, 7)])],
            [list([(2, 0)]), list([(2, 1)]), list([(2, 1), (2, 2)]), list([(2, 2)]),
             list([(3, 3)]), list([(2, 4)]), list([(3, 5)]), list([(3, 6)]),
             list([(3, 7)])],
            [list([(3, 0)]), list([(3, 1)]), list([(3, 1), (3, 2)]), list([(3, 3)]),
             list([(3, 3)]), list([(3, 4), (4, 4)]), list([(3, 5)]),
             list([(3, 6), (4, 6)]), list([(3, 7)])],
            [list([(4, 0)]), list([(4, 0), (4, 1)]), list([(4, 1), (4, 2)]), list([(4, 3)]),
             list([(4, 3)]), list([(4, 4)]), list([(4, 5), (5, 5)]), list([(4, 6)]),
             list([(5, 7)])],
            [list([(5, 0)]), list([(5, 1)]), list([(5, 1), (5, 2)]), list([(5, 3)]),
             list([(5, 4)]), list([(5, 4)]), list([(5, 5)]), list([(5, 7)]),
             list([(5, 7)])]]
        self.assertEqual(assert_my_parents, my_parents)

    def test_get_one_optimal_matches(self):
        seq1 = "TGCTCGTA"
        seq2 = "TTCATA"

        paths = [[(1, 1), (1, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 7), (6, 8)]]
        optimal_matches = get_optimal_matches(seq1, seq2, paths)

        assert_optimal_matches = ["T G C T C G T A \nT - - T C A T A "]
        self.assertEqual(assert_optimal_matches, optimal_matches)

    def test_get_multiple_opimal_matches(self):
        seq1 = "AAGT"
        seq2 = "AGT"

        paths = [[(0, 1), (1, 2), (2, 3), (3, 4)], [(1, 1), (1, 2), (2, 3), (3, 4)]]
        optimal_matches = get_optimal_matches(seq1, seq2, paths)

        assert_optimal_matches = ["A A G T \n- A G T ", "A A G T \nA - G T "]
        self.assertEqual(assert_optimal_matches, optimal_matches)

    def get_arrays_with_only_path_elements(self):
        array = [[0., -6., -12., -18., -24., -30., -36., -42., -48.],
                 [-6., 5., -1., -7., -13., -19., -25., -31., -37.],
                 [-12., -1., 3., -3., -2., -8., -14., -20., -26.],
                 [-18., -7., -3., 8., 2., 3., -3., -9., -15.],
                 [-24., -13., -9., 2., 6., 0., 1., -5., -4.],
                 [-30., -19., -15., -4., 7., 4., -2., 6., 0.],
                 [-36., -25., -21., -10., 1., 5., 2., 0., 11.]]
        paths_coordinates = [[(6, 8), (5, 7), (4, 6), (3, 5), (2, 4), (1, 3), (1, 2), (1, 1), (0, 0)]]

        arrays_with_path = get_arrays_with_only_path_elements(array, paths_coordinates)

        assert_array = np.empty((len(array), len(array[0])))
        assert_array[:] = np.nan
        assert_array[0][0] = 0
        assert_array[1][1] = 5
        assert_array[1][2] = -1
        assert_array[1][3] = -7
        assert_array[2][4] = -2
        assert_array[3][5] = 3
        assert_array[4][6] = 1
        assert_array[5][7] = 6
        assert_array[6][8] = 11

        assert_array = np.array([assert_array])
        self.assertEqual(np.array_equal(arrays_with_path, assert_array, equal_nan=True), True)

    def test_full_algorithm_one_match(self):
        sequence1 = "TGCTCGTA"  # columns
        sequence2 = "TTCATA"  # rows

        optimal_matches_arr = optimal_matches_algorithm(sequence1, sequence2, match=5, mismatch=-2, gap=-6)

        assert_optimal_matches = ["T G C T C G T A \nT - - T C A T A "]
        self.assertEqual(assert_optimal_matches, optimal_matches_arr)

    def test_full_algorithm_two_matches_no_blosum62(self):
        sequence1 = "AAGT"  # columns
        sequence2 = "AGT"  # rows

        optimal_matches_arr = optimal_matches_algorithm(sequence1, sequence2, match=5, mismatch=-2, gap=-6)

        assert_optimal_matches = ["A A G T \n- A G T ", "A A G T \nA - G T "]
        self.assertEqual(assert_optimal_matches, optimal_matches_arr)

    def test_blosum62(self):
        blosum62 = substitution_matrices.load("BLOSUM62")
        print(blosum62)
        self.assertEqual(get_value_from_blosum62("T", "G"), get_value_from_blosum62("G", "T"))

    def test_algorithm_with_blosum62(self):
        seq_40886941 = get_gi_40886941_sequence()
        seq_34849618 = get_gi_34849618_sequence()

        optimal_matches_arr = optimal_matches_algorithm_blosum62(seq_40886941, seq_34849618, gap=-7)
        assert_optimal_matches = ["M V H L T P E E K S A V T A L W G K V N V D E V G G E A L G R L L V V Y P W T Q R L F E S F G D L F T P D A V M G N P K V K A H G K K V L G A F S D G P A H L D N L K G T F A T L S E L H C D K L H V D P E N F R L L G N V L V C V L A H H F G K E F T P P V Q A A Y Q K V V A G V A N A L A H K Y H \nM V H L T D A E K A A V N G L W G K V N P D D V G G E A L G R L L V V Y P W T Q R Y F D S F G D L S S A S A I M G N P K V K A H G K K V I N A F N D G L K H L D N L K G T F A H L S E L H C D K L H V D P E N F R L L G N M I V I V L G H H L G K E F T P C A Q A A F Q K V V A G V A S A L A H K Y H "]
        self.assertEqual(assert_optimal_matches, optimal_matches_arr)

    def test_algorithm_with_blosum_long_assertion(self):
        seq1 = "MLAVLPEKREMTECHLSDEEIRKLNRDLRILIATNGTLTRILNVLANDEIVVEIVKQQIQDAAPEMDGCDHSSIGRVLRRDIVLKGRRSGIPFVAAESFIAIDLLPPEIVASLLETHRPIGEVMAASCIETFKEEAKVWAGESPAWLELDRRRNLPPKVVGRQYRVIAEGRPVIIITEYFLRSVFEDNSREEPIRHQRSVGTSARSGRSICT"
        seq2 = "MTNRTLSREEIRKLDRDLRILVATNGTLTRVLNVVANEEIVVDIINQQLLDVAPKIPELENLKIGRILQRDILLKGQKSGILFVAAESLIVIDLLPTAITTYLTKTHHPIGEIMAASRIETYKEDAQVWIGDLPCWLADYGYWDLPKRAVGRRYRIIAGGQPVIITTEYFLRSVFQDTPREELDRCQYSNDIDTRSGDRFVLHGRVFKNL"
        optimal_matches_arr = optimal_matches_algorithm_blosum62(seq1, seq2, gap=-4)
        assert_match = "M L A V L P E K R E M T E C H L S D E E I R K L N R D L R I L I A T N G T L T R I L N V L A N D E I V V E I V K Q Q I Q D A A P E M D G C D H S S I G R V L R R D I V L K G R R S G I P F V A A E S F I A I D L L P P E I V A S L L E T H R P I G E V M A A S C I E T F K E E A K V W A G E S P A W L E L D R R R - N L P P K V V G R Q Y R V I A E G R P V I I I T E Y F L R S V F E D N S R E E P I - R H Q R S - - V G T - S A - R - - - S G R S I C T - \nM - - T - - N - R - - T - - - L S R E E I R K L D R D L R I L V A T N G T L T R V L N V V A N E E I V V D I I N Q Q L L D V A P K I P E L E N L K I G R I L Q R D I L L K G Q K S G I L F V A A E S L I V I D L L P T A I T T Y L T K T H H P I G E I M A A S R I E T Y K E D A Q V W I G D L P C W L A - D Y G Y W D L P K R A V G R R Y R I I A G G Q P V I I T T E Y F L R S V F Q D T P R E E - L D R C Q Y S N D I D T R S G D R F V L H G R V F K N L "
        self.assertEqual(True, assert_match in optimal_matches_arr)


if __name__ == '__main__':
    unittest.main()
    print("Everything passed")
