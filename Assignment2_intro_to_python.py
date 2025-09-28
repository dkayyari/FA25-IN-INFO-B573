
# dna_analysis.py
# B573 Assignment 2 - DNA Sequence Analysis
# Author: Deeksha Kayyari
# Date: 09/28/2025

# -------------------------------
# Step 1: Read the DNA sequence
# -------------------------------
def read_sequence(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence

seq = read_sequence('chr1_GL383518v1_alt.fa')

#  Print specific letters from the sequence
print("\n🔹 10th letter of sequence:", seq[9])
print("🔹 758th letter of sequence:", seq[757])

# -----------------------------------------------
# Step 2: Create the reverse complement sequence
# -----------------------------------------------
def reverse_complement(seq):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c'
    }
    rev_comp = ''.join([complement.get(base, base) for base in reversed(seq)])
    return rev_comp

rev_seq = reverse_complement(seq)

# Print specific letters from reverse complement
print("\n🔹 79th letter of reverse complement:", rev_seq[78])
print("🔹 500th to 800th letters of reverse complement:\n", rev_seq[499:800])

# ------------------------------------------------------------
# Step 3: Create nested dictionary of nucleotide frequencies
# ------------------------------------------------------------
def count_by_kilobase(seq):
    my_dict = {}
    for i in range(0, len(seq), 1000):
        chunk = seq[i:i+1000]
        base_counts = {}
        for base in set(chunk):
            base_counts[base] = chunk.count(base)
        my_dict[i] = base_counts
    return my_dict

my_dict = count_by_kilobase(seq)

#  Print example: counts for kilobase starting at position 5000
print("\n🔹 Counts for kilobase starting at position 5000:")
print("   ", my_dict.get(5000, {}))

# ------------------------------------------------------------------
# Step 4a: Count A, C, G, T in first 1000 base pairs (case-insensitive)
# ------------------------------------------------------------------
first_kb = my_dict.get(0, {})
first_list = [
    first_kb.get('A', 0) + first_kb.get('a', 0),
    first_kb.get('C', 0) + first_kb.get('c', 0),
    first_kb.get('G', 0) + first_kb.get('g', 0),
    first_kb.get('T', 0) + first_kb.get('t', 0)
]
print("\n🔹 First 1000 bp counts [A, C, G, T]:", first_list)

# ---------------------------------------------------
# Step 4b: Repeat for all kilobases and store results
# ---------------------------------------------------
all_lists = []
for kb in sorted(my_dict.keys()):
    counts = my_dict[kb]
    all_lists.append([
        counts.get('A', 0) + counts.get('a', 0),
        counts.get('C', 0) + counts.get('c', 0),
        counts.get('G', 0) + counts.get('g', 0),
        counts.get('T', 0) + counts.get('t', 0)
    ])

print("\n🔹 Total number of kilobases:", len(all_lists))

# ---------------------------------------------------
# Step 4c: Print nucleotide counts and sum per kilobase
# ---------------------------------------------------
print("\n🔹 Nucleotide counts per kilobase:")
for i, lst in enumerate(all_lists):
    print(f"   • Kilobase {i*1000}: A={lst[0]}, C={lst[1]}, G={lst[2]}, T={lst[3]}")

sums = [sum(lst) for lst in all_lists]
print("\n🔹 Sums of each kilobase:")
print("   ", sums)

# ---------------------------------------------------
# Step 5: Analyze discrepancies in kilobase sums
# ---------------------------------------------------
print("\n🔹 Kilobases with unexpected sums:")
for i, total in enumerate(sums):
    if total != 1000:
        print(f"   ⚠️ Kilobase {i*1000} has unexpected sum: {total}")

# ---------------------------------------------------
# Explanation for Differences in Expected Results
# ---------------------------------------------------
# Each kilobase (chunk of 1000 base pairs) is expected to contain exactly 1000 nucleotides.
# However, some kilobases may have a total count less than 1000. This can happen due to:
#
# 1️. Incomplete Final Chunk:
#     - If the total sequence length isn't a perfect multiple of 1000, the last kilobase will be shorter.
#
# 2️. Non-standard Characters:
#     - DNA sequences may include ambiguous bases like 'N', or symbols not part of A, C, G, T.
#     - These characters are counted but may not be included in the A/C/G/T summary.
#
# 3️. Case Sensitivity:
#     - Although both uppercase and lowercase letters are handled, inconsistencies in input formatting
#       could still affect counts if not properly normalized.
#
# 4️. File Formatting Issues:
#     - Unexpected line breaks, spaces, or hidden characters in the FASTA file could lead to miscounts.
#
# These factors explain why some kilobase sums may differ from the expected value of 1000.

