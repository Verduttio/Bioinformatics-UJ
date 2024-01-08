import unittest

from main import parse_indexes_into_sequence_ids, find_indexes_of_maximum_value_in_2d_list, connect_two_sequences_through_common_sequence


class Tests(unittest.TestCase):
    def test_parse_indexes_into_sequence_ids_1(self):
        self.assertEqual(parse_indexes_into_sequence_ids(0, 1), (0, 2))

    def test_parse_indexes_into_sequence_ids_2(self):
        self.assertEqual(parse_indexes_into_sequence_ids(3, 2), (3, 6))

    def test_parse_indexes_into_sequence_ids_3(self):
        self.assertEqual(parse_indexes_into_sequence_ids(2, 4), (2, 7))

    def test_find_indexes_of_maximum_value_in_2d_list(self):
        self.assertEqual(find_indexes_of_maximum_value_in_2d_list([[1, 2, 3, 4], [5, 60, 7], [8, 9]]), (1, 1))

    # def test_find_indexes_of_two_sequences_with_greatest_common_sequence_len(self):
    #     self.assertEqual(find_indexes_of_two_sequences_with_greatest_common_sequence_len(['AABA', 'CCC', 'ACABAD', 'AADA']), (0, 2))

    def test_connect_two_sequences_through_common_sequence_1(self):
        self.assertEqual('AAHMUYBAD', connect_two_sequences_through_common_sequence('AAHMUY', 'MUYBAD'))

    def test_connect_two_sequences_through_common_sequence_2(self):
        self.assertEqual('BADMUYGOOD', connect_two_sequences_through_common_sequence('MUYGOOD', 'BADMUY'))
