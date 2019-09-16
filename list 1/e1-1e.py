from Bio import SeqIO
records = list(SeqIO.parse("sequence.fasta", "fasta"))

seq = str(records[0].seq[0:])

# Switch case of the complementary DNA
compRib = ''
for n in seq:
    if n == 'A':
        compRib += 'T'
    elif n == 'T':
        compRib += 'A'
    elif n == 'C':
        compRib += 'G'
    elif n == 'G':
        compRib += 'C'
    else:
        compRib += n

# Print of the first 100 nucleotides
print("The first 100 given nucleotides:\n")
print(seq[0:100], "\n")
print("The first 100 nucleotides of the cDNA:\n")
print(compRib[0:100])