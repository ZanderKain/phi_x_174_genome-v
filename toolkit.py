from collections import Counter

def extract_sequences_no_biopython(fasta_file, headers_to_extract=None):
    """
    Extracts sequences from a FASTA file without using Biopython.

    Args:
        fasta_file (str): Path to the FASTA file.
        headers_to_extract (list, optional): A list of headers to extract. 
                                            If None, extracts all sequences. 
                                            Defaults to None.

    Returns:
        dict: A dictionary where keys are the headers and values are the 
              corresponding sequences.  Returns an empty dictionary if there's an error.
    """

    extracted_sequences = {}
    current_header = None
    current_sequence = ""

    try:
        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()  # Remove leading/trailing whitespace

                if line.startswith(">"):  # Header line
                    if current_header: # Save previous sequence if it exists
                        if headers_to_extract is None or current_header in headers_to_extract:
                            extracted_sequences[current_header] = current_sequence
                    current_header = line[1:].split()[0]  # Extract header (remove '>')
                    current_sequence = ""  # Reset sequence for the new header

                else:  # Sequence line
                    current_sequence += line  # Append to the current sequence

            # Save the last sequence (important!)
            if current_header:
                if headers_to_extract is None or current_header in headers_to_extract:
                    extracted_sequences[current_header] = current_sequence

    except FileNotFoundError:
        print(f"Error: FASTA file '{fasta_file}' not found.")
        return {}
    except Exception as e:
        print(f"An error occurred: {e}")
        return {}

    return extracted_sequences



nucleotides = {"A": 0,
               "T": 0,
               "C": 0,
               "G": 0}

DNA_codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}

def countNucs(seq): #Count Nucleotides and seperates by bases from string sequences
    
    total = 0 #Will return total bases in entire sequences if only nucleotides
    notBases = []
    
    for nucs in seq:
        if nucs in nucleotides:
            nucleotides[nucs] += 1
        else:
            notBases.append(nucs) #If bases are not ATCG, appends to a list then.
            
    for key, value in nucleotides.items(): #Checks through list of tuples and adds only the values to the total variable
        total += value
    
        
    return (f'Nucleotides: {nucleotides} \nTotal Nucleotide bases: {total}\nBases not ATCG: {notBases}')       

def gcContent(seq): #Calculates the GC % within DNA sequence
    gCount = 0
    cCount = 0
    total = 0    
    
    for base in seq:    
        if base == 'G':
            gCount += 1
            total += 1
        elif base == 'C':
            cCount += 1
            total += 1
        
        elif base in nucleotides.keys():
            total += 1
        
        else:
            continue
               
    
            
    
    return (f'The GC% is: {((gCount + cCount)/total) * 100:.2f}%') 
        
def translation_total(seq,init_pos=0):
    return f'The total proteins in entire genome: {len([DNA_codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)])} '
    

def meaningful_pro_seq(seq):
    translated_pro = [DNA_codons[seq[pos:pos + 3]] for pos in range(0, len(seq) - 2, 3)]
    
    current_prot = []
    proteins = []
    
    for proseq in translated_pro:
        if proseq == "_":
            
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
                
        else:
            
            if proseq == "M":
                current_prot.append("")
                
            for i in range(len(current_prot)):
                current_prot[i] += proseq
    
    return proteins
                

def genKmers(seq, k):
    kmers = []
    for i in range(len(seq) - k + 1):
        kmers.append(seq[i:i + k])
    return kmers

def count_kmer(seq, k=3):
    kmers = genKmers(seq, k)
    kmer_counts = {}
    
    for kmer in kmers:
        if kmer not in kmer_counts:
            kmer_counts[kmer] = 0
        kmer_counts[kmer] += 1
        
    kmer_c = Counter(kmer_counts)
    total_count = sum(kmer_c.values())
    
    
    return Counter({
    kmer: (count, f'{count / total_count:.3f}') for kmer, count in kmer_counts.items()
})
   
        
        
    
   
    
    