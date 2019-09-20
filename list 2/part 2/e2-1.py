import numpy as np


class Aligner:

    def __init__(self, seq1, seq2, gapPenalty, missPenalty, matchScore):
        self.seqAlignments = []
        self.seq1 = seq1
        self.seq2 = seq2
        self.gapPenalty = gapPenalty
        self.missPenalty = missPenalty
        self.matchScore = matchScore
        self.alignMatrix = np.zeros(
            (len(self.seq1)+1, len(self.seq2)+1), dtype=int)
        self.traceBackMatrix = np.zeros(
            (len(self.seq1)+1, len(self.seq2)+1), dtype='U4')
        # now we have another stop criteria in the traceback matrix: when the score is equals to 0
        self.indexToTrace = {
            0: "d",
            1: "l",
            2: "u",
            3: "f"
        }
        self.finalScore = 0
        self.identity = 0

    
    # get the best possible value according to the recursion formula given
    def getValue(self, i, j):
        if self.seq1[i-1] == self.seq2[j-1]:
            missOrMatch = self.matchScore
        else:
            missOrMatch = self.missPenalty

        possibleValues = [
            self.alignMatrix[i-1][j-1] + missOrMatch,
            self.alignMatrix[i][j-1] + self.gapPenalty,
            self.alignMatrix[i-1][j] + self.gapPenalty,
            0
        ]

        return max(possibleValues), possibleValues.index(max(possibleValues))


    # align the sequences building the matrixes    
    def align(self):

        for row in self.traceBackMatrix:
            row[0] = "f"

        for i in range(len(self.traceBackMatrix[0])):
            self.traceBackMatrix[0][i] = "f"
        self.traceBackMatrix[0][0] = "f"


        for i, j in np.ndindex(self.alignMatrix.shape):
            if i == 0:
                continue
            if j == 0:
                continue

            self.alignMatrix[i][j], index = self.getValue(i, j)
            self.traceBackMatrix[i][j] = self.indexToTrace[index]
        
        # find maximum value(s) from the matrix
        self.finalScore = np.amax(self.alignMatrix)
        # find index of maximum value(2) matrix
        self.maxIndexes = np.where(self.alignMatrix == self.finalScore)

        # and now we need to slice only the parts we are interested in
        self.sliceMatrix()

    
    
    
    def sliceMatrix(self):

        # "zips" the x(es) and y(s) corresponding values    
        listOfCoordinates = list(zip(self.maxIndexes[0], self.maxIndexes[1]))
        
        # make an alignment for each maximum value found (if more the one was found)
        for coord in listOfCoordinates:
            self.slicedAlignMatrix = self.alignMatrix[0:coord[0]+1, 0:coord[1]+1]
            self.makeAlignment()
        
        
    # make the textual alignment
    def makeAlignment(self):

        s1 = ''
        s2 = ''

        i = self.slicedAlignMatrix.shape[0] - 1
        j = self.slicedAlignMatrix.shape[1] - 1
        print("Indexes from score matrix: (%d, %d)" %(i, j))
        alignType = self.traceBackMatrix[i][j]

        while alignType != 'f':

            if alignType == 'd':
                s1 = self.seq1[i-1] + s1
                s2 = self.seq2[j-1] + s2
                i -= 1
                j -= 1

            if alignType == 'l':
                s1 = '-' + s1
                s2 = self.seq2[j-1] + s2
                j -= 1

            if alignType == 'u':
                s1 = self.seq1[i-1] + s1
                s2 = '-' + s2
                i -= 1

            alignType = self.traceBackMatrix[i][j]

        self.s1 = s1
        self.s2 = s2
        self.getIdentity()
        self.printResults()


    # get the number of matches/total number (identity)
    def getIdentity(self):

        ident = 0

        for i in range(len(self.s1)):
            if self.s1[i] == self.s2[i]:
                ident += 1

        totalPositions = max([len(self.seq1), len(self.seq2)])
        self.identity = ident/totalPositions

    

    def printResults(self):
        
        lastI = self.slicedAlignMatrix.shape[0] -1
        lastJ = self.slicedAlignMatrix.shape[1] -1
        
        firstI = lastI - len(self.s1)
        firstJ = lastJ - len(self.s2)
        
        
        # the printed part of the slicedAlignMatrix is only the used matrix part for making the local alignment
        print(self.slicedAlignMatrix[firstI:, firstJ:], '\n')
        print(self.s1)
        print(self.s2, '\n')
        print('Final Score:', self.finalScore)
        print('Identity:', self.identity)
        print('\n\n\n')



# Opens the human sequence
file1 = open("hemoglobins/human.txt","r")
human = file1.read()
file1.close()


# Opens the biomphalaria sequence
file1 = open("hemoglobins/biomphalaria.txt","r")
biomphalaria = file1.read()
file1.close()


# Assign the task given values 
gap = -2
match = 1
missmatch = -1


# Align process (prints inside the method's callings)
aligner = Aligner(human, biomphalaria, gap, missmatch, match)
aligner.align()
