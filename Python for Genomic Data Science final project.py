# Final project for the Python for Genomic Data Science course
# David - 220619
# Building functions for DNA fasta file reading, ORF finding, and repeat ID'ing

def fasta_readv2(file):
    # This function reads all sequences in a fasta file and assigns them in a dictionary
    fasta_file = open(file,"r")
    lines = fasta_file.readlines()
    last = lines[-1]
    fasta_names = []
    fasta_seqs = []
    seq = ""
    for line in lines:
        if line[0] == ">":
            fasta_seqs.append(seq)
            seq = ""
            fasta_names.append(line[1:len(line)-1])
        else:
            seq = seq + line[0:len(line)-1]
        if line == last:
            fasta_seqs.append(seq)
            fasta_seqs.pop(0)
    fasta_file.close()
    return dict(zip(fasta_names, fasta_seqs))
fasta_files_dict = fasta_readv2(file = ".\dna.example.fasta")
len(fasta_files_dict)

# Lengths of fasta entries
fasta_lengths = []
for key in fasta_files_dict:
    fasta_lengths.append(len(fasta_files_dict[key]))
max(fasta_lengths)
min(fasta_lengths)

# Question 3: ORFs
def open_reading_frames(DNA, strand = 0):
    # strand = 0 means all reading frame
    # strand = 0 means only reading frame 1
    # and so on
    if strand == 0:
        start = 0
        by_n = 1

    if strand == 1:
        start = 0
        by_n = 3

    if strand == 2:
        start = 1
        by_n = 3

    if strand == 3:
        start = 2
        by_n = 3

    start_codon = "ATG"
    stop_codon = ["TAA", "TAG", "TGA"]
    ORF = {}
    ORFnumber = 1
    for pos in range(start,len(DNA),by_n):
        if DNA[pos:pos+3] == "ATG":
            seq = ""
            for sub_pos in range(pos,len(DNA),3):
                codon = DNA[sub_pos:sub_pos+3]
                seq = seq + codon
                if codon in stop_codon:
                    ORF["orf_" + str(ORFnumber) + "_pos_" + str(pos+1)] = seq
                    ORFnumber += 1
                    break
    return ORF
open_reading_frames("ATGAGGTGACATGACCGTAACAAGCCTTATGATTAGCTAA",strand = 1)
orfs = open_reading_frames("ATGAGGTGACATGACCGTAACAAGCCTTATGATTAGCTAA",strand = 2)

# Count the lengths
orf_lengths = []
for key in orfs:
    orf_lengths.append(len(orfs[key]))

# Function 4: identifying DNA repeats
def dna_repeats(DNA,n):
    nmers = []
    dict_nmers = {}
    for pos in range(0,len(DNA)):
        nmer_n = DNA[pos:pos+n]
        if len(nmer_n) < n: 
            break
        if nmer_n in nmers:
            dict_nmers[nmer_n] = nmers.count(nmer_n) + 1
        nmers.append(nmer_n)
    return dict_nmers
dna_repeats("ATGGTGCTGTGCGTGCATGG",4)

# FINAL QUIZ
# Questions 1,2 and 3
# How many records are in the multi-FASTA file?
# What is the length of the longest sequence in the file?
# What is the length of the shortest sequence in the file?
quiz_fasta = fasta_readv2(".\dna2.fasta")
len(quiz_fasta)
fasta_lengths = []
for key in quiz_fasta:
    fasta_lengths.append(len(quiz_fasta[key]))
max(fasta_lengths)
min(fasta_lengths)

# Question 4 What is the length of the longest ORF appearing in reading frame 2 of any of the sequences?
def find_longest_orf(fasta_dict, frame = 0):
    longest_orf ={}
    for key in fasta_dict:
        orfs = open_reading_frames(fasta_dict[key], strand = frame)
        orf_lengths = []
        for key2 in orfs:
            orf_lengths.append(len(orfs[key2]))
        if orf_lengths:
            longest_orf[key] = max(orf_lengths)
    return longest_orf

orf_lengths = find_longest_orf(quiz_fasta,frame = 2)
max_key = max(orf_lengths, key=orf_lengths.get)
max_key
orf_lengths[max_key]

# Question 5 - What is the starting position of the longest ORF in reading frame 3 in any of the sequences? The position should indicate the character number where the ORF begins. For instance, the following ORF:
orf_lengths = find_longest_orf(quiz_fasta,frame = 3)
max_key = max(orf_lengths, key=orf_lengths.get)
orf_lengths[max_key]

orfs_in_max3 = open_reading_frames(quiz_fasta[max_key],3)

orf_len = []
or_name = []
for key in orfs_in_max3:
    orf_len.append(len(orfs_in_max3[key]))
    or_name.append(key)
or_name[3]



# Question 6 - What is the length of the longest ORF appearing in any sequence and in any forward reading frame?
orf_lengths = find_longest_orf(quiz_fasta,frame = 0)
max_key = max(orf_lengths, key=orf_lengths.get)
orf_lengths[max_key]

# Question 7 - What is the length of the longest forward ORF that appears in the sequence with the identifier  gi|142022655|gb|EQ086233.1|16?
id = "gi|142022655|gb|EQ086233.1|16 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence"
q7_orfs = open_reading_frames(quiz_fasta[id])
q7_orf_lengths = {}
for key in q7_orfs:
    q7_orf_lengths[key] = len(q7_orfs[key])
max_key = max(q7_orf_lengths, key=q7_orf_lengths.get)
q7_orf_lengths[max_key]

# Question 8 - Find the most frequently occurring repeat of length 6 in all sequences. How many times does it occur in all?
repeats_dict = {}
for key in quiz_fasta:
    repeats = dna_repeats(quiz_fasta[key],6)
    for repeat_name in repeats:
        if not repeat_name in repeats_dict:
            repeats_dict[repeat_name] = repeats[repeat_name]
        else:
            repeats_dict[repeat_name] = repeats_dict[repeat_name] + repeats[repeat_name]
len(repeats_dict)
# We have 2392 distinct 6-mer DNA repeats
max_key = max(repeats_dict, key=repeats_dict.get)
repeats_dict[max_key]

# Question 9 - Find all repeats of length 12 in the input file. Let's use Max to specify the number of copies
# of the most frequent repeat of length 12.  How many different 12-base sequences occur Max times?
repeats_dict = {}
for key in quiz_fasta:
    repeats = dna_repeats(quiz_fasta[key],12)
    for repeat_name in repeats:
        if not repeat_name in repeats_dict:
            repeats_dict[repeat_name] = repeats[repeat_name]
        else:
            repeats_dict[repeat_name] = repeats_dict[repeat_name] + repeats[repeat_name]
len(repeats_dict)
# We have 2392 distinct 6-mer DNA repeats
max_key = max(repeats_dict, key=repeats_dict.get)
repeats_dict[max_key]

# Question 10 - Which one of the following repeats of length 7 has a maximum number of occurrences?
repeats_dict = {}
for key in quiz_fasta:
    repeats = dna_repeats(quiz_fasta[key],7)
    for repeat_name in repeats:
        if not repeat_name in repeats_dict:
            repeats_dict[repeat_name] = repeats[repeat_name]
        else:
            repeats_dict[repeat_name] = repeats_dict[repeat_name] + repeats[repeat_name]
len(repeats_dict)
repeats_dict["CGCGCCG"]
repeats_dict["TGCGCGC"]
repeats_dict["CATCGCC"]
repeats_dict["GCGGCCG"]
