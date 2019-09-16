from Bio import SeqIO
records = list(SeqIO.parse("sequence.fasta", "fasta"))

seq = str(records[0].seq[0:])
nucleotides = ['A', 'C', 'T', 'G']

finalCount = {}

# Count nucleotides and other possible non-nucleotide characters 
for n in seq:
    if n in finalCount:
        finalCount[n] += 1
    else:
        finalCount[n] = 1

# print the final count result
for n in finalCount:
    if n in nucleotides:
        print(n, "- amount found:", finalCount[n])
    else:
        print("Non-nucleotide found:", n, "- amount found:", finalCount[n])