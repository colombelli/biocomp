from Bio import SeqIO
records = list(SeqIO.parse("sequence.fasta", "fasta"))

seq = str(records[0].seq[0:])

seqOcur = {}

i = 0
j = 37  
seqLen = len(seq)

# Build a dictionary with key meaning the sequence and value meaning the amount of times that it appears
while j <= seqLen:
    
    if (seq[i:j] in seqOcur):
        seqOcur[seq[i:j]] += 1
    else:
        seqOcur[seq[i:j]] = 1
    
    i += 1
    j += 1

# Print the first 100 different sequences with their respective count
print(len(seqOcur), "different subsequences of size 37 were found. \n\n")
limit = 1
print("The first 100 sequences with their respective count:")
for key, value in seqOcur.items():
    if limit == 101:
        break
    print("Sequence:", key, "Count:", value) 
    limit += 1