from copy import deepcopy
from math import log2

def select_test_set():
    test_set = input("Select (1 or 2) the test set... ")
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




def find_motifs(sequences, motif_len):

    matrices =  build_first_matrices(sequences[0], motif_len)
    for sequence in sequences[1:]:
        matrices =  update_matrices(sequence, matrices, motif_len)

    best_matrix = select_matrix(matrices)
    return best_matrix
        



def build_first_matrices(sequence, motif_len):

    matrices = []
    num_matrices = len(sequence) - motif_len + 1

    for i in range(num_matrices):
        l_mer = get_l_mer_by_iteration(sequence, motif_len, i)
        matrix = build_matrix([l_mer])
        matrices.append(matrix)
    
    return matrices



def get_l_mer_by_iteration(sequence, motif_len, i):
    start_idx = i
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
        matrix.append(row)

    I = compute_I(matrix, l_mers)

    return (l_mers, matrix, I)
             

def compute_I(matrix, l_mers):
    
    num_sequences = len(l_mers)
    pa = get_genomic_frequency(l_mers, 'a')
    pc = get_genomic_frequency(l_mers, 'c')
    pg = get_genomic_frequency(l_mers, 'g')
    pt = get_genomic_frequency(l_mers, 't')
    frequencies = [pa, pc, pg, pt]

    I = 0

    for idx, row in enumerate(matrix):
        for col in row:
            N_bi = col
            if N_bi == 0:
                continue

            I += (N_bi / num_sequences) * \
                log2((N_bi / num_sequences) / frequencies[idx])

    return I


def get_genomic_frequency(l_mers, base):

    num_bases = len(l_mers) * len(l_mers[0])
    base_occurrences = 0
    for l_mer in l_mers:
        base_occurrences += l_mer.count(base)
    
    return base_occurrences/num_bases



def update_matrices(sequence, matrices, motif_len):

    # For each matrix
    num_matrices = len(matrices)
    updated_matrices = []

    for matrix in matrices:
        current_matrices = []
        
        for i in range(num_matrices):
            l_mers = deepcopy(matrix[0])
            new_l_mer = get_l_mer_by_iteration(sequence, motif_len, i)
            l_mers.append(new_l_mer)
            matrix_temp = build_matrix(l_mers)
            current_matrices.append(matrix_temp)
        
        best_matrix = select_matrix(current_matrices)
        updated_matrices.append(best_matrix)
    
    return updated_matrices



def select_matrix(matrices):

    best_score = 0
    best_matrix = []

    for matrix in matrices:
        if matrix[2] > best_score:
            best_score = matrix[2]
            best_matrix = matrix
    
    return best_matrix


def print_matrix(matrix):
    print("Matrix:")
    for row in matrix:
        print(row)


def print_patterns(patterns):
    print("Patterns:")
    for pattern in patterns:
        print(pattern)



if __name__ == "__main__":
    
    sequences, motif_lens = select_test_set() 
    while not sequences:
        sequences, motif_lens = select_test_set() 

    print("\n\n#######################################\n")
    for motif_len in motif_lens:
        
        matrix_tuple = find_motifs(sequences, motif_len)
        print("\nMotif length considered:", motif_len)
        print("\n")
        print_patterns(matrix_tuple[0])
        print("\n")
        print_matrix(matrix_tuple[1])
        print("\nI score:", matrix_tuple[2])
        print("\n\n#######################################\n")



