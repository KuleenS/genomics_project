import random
import numpy as np

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

# Parameters for testing
length = 1000  # Length of each string
num_changes = 50  # Number of insertions/deletions/substitutions
num_pairs = 2  # Number of test string pairs to generate

# Generate random string pairs for benchmarking
random_test_pairs = generate_test_strings(length, num_changes, num_pairs)
for i in range(num_pairs):
    test_file_name = f"testfile{i+1}.txt"
    with open(test_file_name, 'w') as outfile:
        outfile.write(f"{random_test_pairs[i][0]}\n{random_test_pairs[i][1]}")
