from Bio import SeqIO
records = list(SeqIO.parse("sequence.fasta", "fasta"))

seq = str(records[0].seq[0:])


# function that counts palindrome of given size 
def countPalindromes(paSize):
    i = 0
    j = paSize
    
    seqLen = len(seq)
    counted = 0

    while j <= seqLen:
        currentSeq = seq[i:j]
        if(currentSeq==currentSeq[::-1]):
            counted += 1

        i += 1
        j += 1

    print("Palindromes of size %d:" %(paSize), counted)

countPalindromes(11)
countPalindromes(9)