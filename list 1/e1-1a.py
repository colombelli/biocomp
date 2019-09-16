from Bio import SeqIO
records = list(SeqIO.parse("sequence.fasta", "fasta"))

seq = str(records[0].seq[0:])

searchingSeq = 'CAATTGAATAATTG'
nucleotides = ['A', 'C', 'T', 'G']

possibleSeq = {}

# generates all possible sequences
for idx in range(len(searchingSeq)):
    for n in nucleotides:
        newSeq = list(searchingSeq)
        newSeq[idx] = n
        strNewSeq = ''.join(newSeq)
        if not(strNewSeq in possibleSeq):
            possibleSeq[strNewSeq] = [idx, 0]

# counts the occurence of every possible sequence
for s in possibleSeq:
    occurrences = seq.count(s)
    idx, oc = possibleSeq[s]
    possibleSeq[s] = [idx, occurrences]

newFinalDict = {}

# isolates in another dictionary only the sequences appearing one time in the chromosome 7
for s in possibleSeq:
    idx, oc = possibleSeq[s]
    if oc == 1:
        newFinalDict[s] = [idx+1, oc] 

# print the unique sequences found, the ith nucleotide that has been modified, 
# the original given nucleotide with the changed one,
# and the starting index of the subsequence in the chromosome 7 sequence
print("Original  given subsequence:", searchingSeq)
print("-----------------------------------------------------------------\n")

for key, val in newFinalDict.items():
    givenNucleotide = searchingSeq[val[0]-1]
    nucleotideChanged = key[val[0]-1]
    print("Unique sequence:", key) 
    print("ith position that changed:", val[0])
    print("Starting index position to find the subsequence in the original chromosome 7 sequence:", seq.find(key))
    print("Original (mutated) given nucleotide:", givenNucleotide)
    print("Changed (theoretical unique and non-mutated) nucleotide:", nucleotideChanged)
    print("-----------------------------------------------------------------\n")