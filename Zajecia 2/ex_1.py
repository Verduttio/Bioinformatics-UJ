def hamming_distance(sequence_1st, sequence_2nd):
    distance = 0
    sequences = zip(sequence_1st, sequence_2nd)
    for (symbol_seq_1st, symbol_seq_2nd) in sequences:
        if symbol_seq_1st != symbol_seq_2nd:
            distance += 1
    return distance


if __name__ == '__main__':
    sequence_1st = "TGTTTGTAGAGTCAAACCGGGCATTGTGCTACAACTAGGGTTGTCGGGTCCTCCTGTTAAACCGTTAAAGCATCGAGTTTATTATCGGACGTGGAGAGTTGAAGCAGTGCAACACAGGCCTAGCTTTTGGGCACTCCACCGTACTGTATTCGGCTTCTCGACGTTGAAGTTAATGAGACGACGGTCAACTAAGTAACGAGGATTGCGAAAGCCATCTCGCATTGAGTGTCCGCACCCTGGCTAAGGGAAAGCAAACCGACGTTCTGTATACATTGCCGGCACTAGTAATGTTCGTAAAAACTGCTGCTCGCTGTCATTCTAGGTTGAAGTCATATCTCCGAGCGACCCGTACAGGACACATGTAGGCTTCACCTTGACGGTCACCGGATTGTCCTCAACCCCTCATGATAGATGAAATGACTCAAGCCCTCGCCGCCGGGCAAACCAATTATAGAGACGAGAACCGTACGTATCGACGCAGCCTCGCAGTAGTCCAGTATTTCAACACCTACCGCCCCACTGGACGTAAAACCAGTATGGTACGTGGCCTCGGTTTGCAGTCTTGTCAAACGAAAAGTGTGTCTTAGTTGACCTCGTAGACGCGAAAATGTACTCGACAACTGGTCGTGTTATGAGTGCAATGAGTTTTGTGACCCCGGCTCTGCACCGGGCCGGTTCTTAATATGTCGGTCCAGTCGTGGCCACAGTTACATCCGGGGGCAATATGTACATCGATAACTGGCTGGCGCCAGTCCTATGGACGACGGTTTTAGCTGCCGGTGCGTTACGCGTACGCACACCTCGACTGCCACTGGTGTCGGCGGGCTTGGCGACCTGCACATGTTCGATCTCCGTGGCCGAATGACCAAGTCCATAATCGTCGGACATACTCCCTCTCCTCGTTGCGGAACGAAGGATCTGGCCATGGCTTGGACCCGACCGTTGTGA"
    sequence_2nd = "TGCTTGCAGTGTTCAGTCCGGCCCTATGCTATCAGTAGGATTGGTCCTTAGTCGCGCTCTACGGTTGAATCATCGTCTTTATTATCTGACGAGGCGAGTGGCATCGGTGCAGCCTGGCCGTCATTCTTGAACTCCACGGAGAACGCGATTACGCTCGTAGCCGCAAAATTTAATGAGCTGGTAGATAGCTAAGTAATTAAGATGGTATAGGAAATCTCTTATTCTGGTGCCCGACGCAAGTTGTAAGTAACTAACTCGTCCTGAAGTCTATACTTCCTCCCATGTTGATTTGCGTTAACAGCGCCGCTATCTGCCCTTCGACGTCGACGCCATATGAACGAGCTACCCATGCATCCCTCGGCTAGGTATCATCTAGGCGGCCTCCCAAATGGCTTGAACGATTGATCCGACCCTCAAGGAAACAGGCTGACGTCCTGCGATGCACTAGGCATATAGGCTCGGACCATCTGTAGAGTGGCTATCCCCCACCGATCCATTGTACCAAATAATACAGCTCCAGCTGACTTGGCATTAACACGGCACAAGGACTGTGCAACTAGTCTGGTCACGTGACAGGTGTAGCTTGGTTGACCTTAAAGTCTGAGGAACGTACCCCTACACTGGGTTTGTTATGGGGGCCTTGAATTTTCAGAAACTAGCCCCGCTAGAGATAGATTCAATGCATTTCTGTTTTCTGATAGCCACAGCTACGTCGCATTAGGCGATGCACAAAGCATATAAGTTGACGCGCGCGCCATTCTCTACAGCCTTAGGAGGATTAATGTCAAGGTCTATTCAGAACAGACGGCCACTGGCGTGGTCACGCATTGACTCTTGCACATGTCCGATCGGCTTGGCGGTTTGGCCAGGGAGAAGATCGTCAGGCCCTAATTCTCACTCGCTTACGTAACCAAGGATCTCTCTATGGGAGTTACCGCAAGGTTATGA"
    print(hamming_distance(sequence_1st, sequence_2nd))

