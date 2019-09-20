import numpy as np
import matplotlib.pyplot as plt


class Aligner:

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.pMatrix = np.zeros(
            (len(self.seq1)+1, len(self.seq2)+1), dtype=int)
        self.mutations = 0

    # align() simply iterates trhough the matrix putting 1 where there's a match and 0 otherwise
    def align(self):

        for i, j in np.ndindex(self.pMatrix.shape):
            if i == 0:
                continue
            if j == 0:
                continue

            if (self.seq1[i-1] == self.seq2[j-1]):
                self.pMatrix[i][j] = 1
            else:
                self.pMatrix[i][j] = 0

                # also checks if it's in a diagonal to look for the mutation
                if i == j:
                    self.mutations += 1
                    print("Mutation found!\nOriginal amino acid: " +
                          self.seq1[i-1] + "\nChanged to: " + self.seq2[j-1] + "\nAt position: " + str(i) + "\n\n")
                # post-run note: this doesn't work because the sequence has a gap and there's a point where
                # most of the sequences are different and i currently don't know how to detect another imperfect
                # diagonal algorithmically while/after building the dot matrix

file1 = open("hemoglobins/pkp2.txt", "r")
pkp2 = file1.read()
file1.close()

file1 = open("hemoglobins/pkp2-mutated.txt", "r")
pkp2m = file1.read()
file1.close()


aligner = Aligner(pkp2, pkp2m)
aligner.align()
print("Number of mutations found:", aligner.mutations)

# plots and saves the matrix as a .png file
plt.figure(figsize=(15, 15))
plt.matshow(aligner.pMatrix, cmap='gray', fignum=1)
plt.savefig('dot_matrix.png')
