from toolkit import * 

with open('sequence.fasta', 'r') as file:
    
    fastaRead = file.readlines()

# print(fastaRead)
    
fastaLabel = ""
fastaSeq = ""

# Append first line from fasta text to the fastalabel variable
# Append any other line that contains sequence and concatonate it together.

for line in fastaRead:
    if line[0] == ">":
        fastaLabel += line

    else:
        fastaSeq += line

fastaSeq= fastaSeq.replace('\n','') #removes the newline operator from text
print(fastaLabel)       
# print(fastaSeq)


print(countNucs(fastaSeq))
print(gcContent(fastaSeq))
print(translation_total(fastaSeq))
print(count_kmer(fastaSeq,))

#Writing out a txt file of the protein output
with open("protein_output.txt", 'w') as f:
    for pro_seq in meaningful_pro_seq(fastaSeq):
        f.write(pro_seq + "\n")
 