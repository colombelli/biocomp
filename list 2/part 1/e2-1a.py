import numpy as np


class Aligner:

    def __init__(self, seq1, seq2, gapPenalty, missPenalty, matchScore):
        self.seq1 = seq1
        self.seq2 = seq2
        self.gapPenalty = gapPenalty
        self.missPenalty = missPenalty
        self.matchScore = matchScore
        self.alignMatrix = np.zeros(
            (len(self.seq1)+1, len(self.seq2)+1), dtype=int)
        self.traceBackMatrix = np.zeros(
            (len(self.seq1)+1, len(self.seq2)+1), dtype='U4')
        self.indexToTrace = {
            0: "d",
            1: "l",
            2: "u"
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
            self.alignMatrix[i-1][j] + self.gapPenalty
        ]

        return max(possibleValues), possibleValues.index(max(possibleValues))

    # align the sequences building the matrixes
    def align(self):

        for row in self.traceBackMatrix:
            row[0] = "u"

        for i in range(len(self.traceBackMatrix[0])):
            self.traceBackMatrix[0][i] = "l"
        self.traceBackMatrix[0][0] = "f"

        accGap = 0
        for row in self.alignMatrix:
            row[0] = accGap
            accGap += self.gapPenalty

        accGap = 0
        for i in range(len(self.alignMatrix[0])):
            self.alignMatrix[0][i] = accGap
            accGap += self.gapPenalty

        for i, j in np.ndindex(self.alignMatrix.shape):
            if i == 0:
                continue
            if j == 0:
                continue

            self.alignMatrix[i][j], index = self.getValue(i, j)
            self.traceBackMatrix[i][j] = self.indexToTrace[index]

        self.finalScore = self.alignMatrix[len(self.seq1)][len(self.seq2)]
        self.makeAlignment()

    # make the textual alignment
    def makeAlignment(self):

        s1 = ''
        s2 = ''

        i = len(self.seq1)
        j = len(self.seq2)
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

    # get the number of matches/total number (identity)
    def getIdentity(self):

        ident = 0

        for i in range(len(self.s1)):
            if self.s1[i] == self.s2[i]:
                ident += 1

        totalPositions = max([len(self.seq1), len(self.seq2)])
        self.identity = ident/totalPositions

    def printResults(self):
        print(self.alignMatrix, '\n')
        print(self.s1)
        print(self.s2, '\n')
        print('Final Score:', self.finalScore)
        print('Identity:', self.identity)


# open the human sequence
file1 = open("hemoglobins/human.txt", "r")
human = file1.read()
file1.close()

animals = {}

animalList = ["chicken", "cow", "deer", "horse", "pig", "trout", "wolf"]
animals = {}

# open each animal sequence
for animal in animalList:
    file1 = open("hemoglobins/"+animal+".txt", "r")
    animalSequence = file1.read()
    file1.close()
    animals[animal] = animalSequence


# default values given by the exercise
gap = -4
match = 5
missmatch = -3

scores = {}


# compare all animal hemoglobin sequences against the human one
for animal in animalList:

    aligner = Aligner(human, animals[animal], gap, missmatch, match)
    aligner.align()
    scores["human vs "+animal] = (aligner.finalScore, aligner.identity)

    print("Results for human vs "+animal+":\n")
    aligner.printResults()
    print("\n")


# get the best score obtained by the comparsions using identity as a tie breaker
bestVal = (0, 0)
bestKey = ''
draws = []

for key, value in scores.items():

    if value[0] > bestVal[0]:
        bestVal = value
        bestKey = key
        draws = []

    elif value[0] == bestVal[0]:
        if value[1] > bestVal[1]:
            bestVal = value
            bestKey = key
            draws = []

        elif value[1] == bestVal[1]:
            draws.append(key)

if draws == []:
    print("\nThe best score was achieved when "+bestKey +
          " were compared. \nObtained score:", scores[bestKey][0], "\nObtained identity:", scores[bestKey][1])

else:
    print("\nMore than one sequence aligment produced the same final score and identity.\n")
    print(bestKey, "- (score, identity):", scores[bestKey])

    for d in draws:
        print(d, "- (score, identity):", scores[d])
