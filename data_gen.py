import random
import numpy as np
from Bio import SeqIO

def random_string(length, alphabet="ATCG"):
    #Generates a random DNA string of specified length from the given alphabet
    return ''.join(random.choices(alphabet, k=length))

def perturbed_string(original, num_changes, alphabet="ATCG"):
    #applies random insertions/deletions/substitutions to the original string
    result = list(original)
    changes = 0
    while changes < num_changes:
        change_type = random.choice(['insert', 'delete', 'substitute'])
        pos = random.randint(0, len(result) - 1)
        if change_type == 'insert':
            result.insert(pos, random.choice(alphabet))
        elif change_type == 'delete' and len(result) > 1:
            result.pop(pos)
        elif change_type == 'substitute':
            result[pos] = random.choice(alphabet) 
        changes += 1
    
    return ''.join(result)

def generate_test_strings(length, num_changes, num_pairs=5, alphabet="ATCG"):
    #returned all generated test strings
    pairs = []
    for i in range(num_pairs):
        original = random_string(length, alphabet)
        perturbed = perturbed_string(original, num_changes, alphabet)
        pairs.append((original, perturbed))
    return pairs

def generate_genomic_strings(fna_file, sequence_length, num_changes, num_pairs):
    # Load a reference genome
    pairs = []
    with open(fna_file, "r") as file:
        # Parse the FASTA file and load sequences into memory
        records = list(SeqIO.parse(file, "fasta"))
        
        for i in range(num_pairs):
            # Select a random chromosome or scaffold
            record = random.choice(records)
            chrom_length = len(record.seq)
            
            # Randomly select a subsequence from the chosen record
            start_pos = random.randint(0, chrom_length - sequence_length)
            original = str(record.seq[start_pos:start_pos + sequence_length]).upper()
            
            # Generate a perturbed version of the subsequence
            perturbed = perturbed_string(str(original), num_changes)
            pairs.append((str(original), perturbed))
    return pairs

# Parameters for testing
#fna_file = "/Users/olivialanglois/Library/CloudStorage/OneDrive-JohnsHopkins/sequences/drosophilia_genome_chromosome_X.fna"
length = 10  # Length of each string
num_changes = 5  # Number of insertions/deletions/substitutions
num_pairs = 2  # Number of test string pairs to generate

# Generate random string pairs for benchmarking
test_pairs = generate_test_strings(length, num_changes, num_pairs)

#generate genome string pairs for benchmarking
#test_pairs = generate_genomic_strings(fna_file, length, num_changes, num_pairs)


for i in range(num_pairs):
    test_file_name = f"testfile{i+1}.txt"
    with open(test_file_name, 'w') as outfile:
        outfile.write(f"{test_pairs[i][0]}\n{test_pairs[i][1]}")
