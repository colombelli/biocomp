def select_test_set():
    test_set = input("Select (1 or 2) the test set...")
    sequences = False
    motif_lens = False

    if test_set == "1":
        # First test set
        s1 = "cccctgatagacgctatctggctatccacgtacataggtcctctgtgcgaatctatgcgtttccaaccat"
        s2 = "agtttactggtgtacatttgatacgtacgtacaccggcaacctgaaacaaacgctcagaaccagaagtgc"
        s3 = "aaaggagtccgtgcaccctctttcttcgtggctctggccaacgagggctgatgtataagacgaaaatttt"
        s4 = "agcccctccgatgtaagtcatagctgtaactattacctgccacccctattacatcttacgtacgtataca"
        s5 = "ctgggttatacaacgcgtcatggcggggtatgcgttttggtcgtcgtacgctcgatcgttaacgtaggtc"
        sequences = [s1,s2,s3,s4,s5]
        motif_lens = [8, 5, 3]


    elif test_set == "2":
        # Second test set
        s1 = "gtcacgcttctgcataccatcctgactactcgtggcgaatacggttcgtctcagaacattgacgagtaggacctccatgtacacgtgagttcgccagtagagggcagaactagaggcccgagctcgttacccagtatatgtactcggcacacactgggatataatactacacgggatactaatagtggcatatcacgccg"
        s2 = "atccctctaacaagttgttttgacggaccgtatttccaaatgtgctcggcttcagaaacaacctttctgccctctactggcgacgtcacaacgacgacaacagaccatatggagtggaccctactcatgtaattgagaccgtcgcatgtagttgatttatgtaaacatatggctctagtttcaggcccctgtaaaggtaa"
        s3 = "ttacataggttccttcacgtcactccttgtccgcgatatctcctcttacccttactaccaagcgtttcctgaaaggcaatgaaaagttgccatgcgctgtcgccagtagagggcagaataccaaggcgcttcagacaactgtcgctgttcgtgggtgggagggattgtatctataatataggatagttcgtatcgaaaaa"
        s4 = "ttatcgaccgccactttctcgccagtagagggcagaaccacaaagtgactccccgagcaatggctgacctactagttatccggcatcacatcggcacatatacgggcgagaccgagccctctccgtaaccaccagtcccactacttcacaggcatatcctgtatcaatgaaatcacaaacgttcgcatgaagataatcgt"
        s5 = "gacggcacattttaacggcccaggttggcagacaggaaatctacgatggtgctactgctttcccgagctctgccacgatgccacacagcacaattctgccctctactggcgatcacctcgaataaaaccgaatgcaagaccgagtaacagcggctggtaacatgcgggaggacgcgctttccgcaagtatattaataggt"
        s6 = "tgcatcataggttagtaagagttataaatcttcgatccctaagtgtggtgcactactcggtcgacctcgcattgacacaacgcgagagtcgccagtagagggcagaagccggcacttttgacctcttctatagaaggtagaccgtgagatcgcgcccgaaggggcccgacggtctccaaggtggaacgtattaggtaatc"
        s7 = "gccgtgtataggcctccgatcgtgcggttctgccctctactggcgaaaggggcatttgctattccaatcgcatagattaccaaataaaaaacgaaagaaggccgtccttgcaaagcttagtccttaaactgagatgcttggcgaccggccataagctccactcgcttgagcacatcaccaagaatcaaagtagcaaaccc"
        s8 = "tgtccgctctgcagacgtccgggatctacgttggtgttctctctagtaacagtacggcagttctttttcggtccgaagcgaatcccacccgccaaggttacataagcattatctgaaggcaaccatacgaactctcattggctcgccagtagagggcagaaggacatcgtgtcatagcacatgcccacagaggagattcg"
        s9 = "ccggtctcaatagccgaacgaggatcgactggtaggcgtgtcgggtgtgtgtggaccggctttggaagaaccacacttctctggccctcaattggccaaaagtcatcttaaggatcggttggccggcagcccctcgtgactacgataaccccaggttctgccctctactggcgaccttgacagagcacttacccactgta"
        s10 = "aaaagagtagtgatgagttagaaagaatttaagggacatcctcttgatttggacggctatccccaggaaatcgtaggcgggggtgcacatggatatctttaggtattaattcccccattccctcttcgttctgccctctactggcgaatgttgctcgcaatactacagcctcctcaatacaggtagggtattttacatat"
        sequences = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10]
        motif_lens = [5, 3]
    
    return sequences, motif_lens




if __name__ == "__main__":
    
    sequences, motif_lens = select_test_set() 
    while not sequences:
        sequences, motif_lens = select_test_set() 


    for motif_len in motif_lens:
        find_motifs(sequences, motif_len)



def find_motifs(sequences, motif_len):

    matrices =  build_first_matrices(sequences[0], motif_len)



def build_first_matrices(sequence, motif_len):

    matrices = []
    num_matrices = len(sequence) - motif_len + 1

    for i in range(num_matrices):
        l_mer = get_l_mer_by_iteration(sequence, motif_len, i)
        matrix = build_matrix([l_mer])
        matrices.append(matrix)
    
    return matrices



def get_l_mer_by_iteration(sequence, motif_len, i):
    start_idx = motif_len * i
    end_idx = start_idx + motif_len
    return sequence[start_idx:end_idx]



def build_matrix(l_mers):

    num_cols = len(l_mers[0])
    matrix = []

    for base in ['a','c','g','t']:
        row = []
        for position in range(num_cols):
            matches = 0
            for l_mer in l_mers:
                if l_mer[position] == base:
                    matches += 1
            row.append(matches)

    return (l_mers, matrix)
             



    