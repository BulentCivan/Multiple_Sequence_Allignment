import numpy as np
from scipy.cluster._optimal_leaf_ordering import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import pandas as pd

def read_Blosum62():
    Blosum62_file_path = 'Blosum62.txt'
    scoring_matrix = {}
    with open(Blosum62_file_path, 'r') as file:
        amino_acids = file.readline().split()
        for line in file:
            values = line.split()
            amino_acid_row = values[0]
            for i, value in enumerate(values[1:]):
                amino_acid_col = amino_acids[i]
                scoring_matrix[(amino_acid_row, amino_acid_col)] = int(value)
    return scoring_matrix

def read_fasta_file(file_path):
    sequence_dict = {}
    current_sequence_id = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence_id = line[1:]
                sequence_dict[current_sequence_id] = ''
            else:
                sequence_dict[current_sequence_id] += line
    return sequence_dict

def getting_gap_penalty():
    while True:
        try:
            gap_penalty = float(input("Please enter the gap penalty: "))
            break
        except ValueError:
            print("Error: Please enter a valid number!")
    return gap_penalty


similarity_table = {}
def pairwise_sequence_alignment(input,blosum62, gap_penalty):
    sequence_names = list(input.keys())

    for k in range(len(sequence_names)):
        for l in range(k + 1, len(sequence_names)):

            seq1_length = len(input[(sequence_names[k])])
            seq2_length = len(input[(sequence_names[l])])
            seq1=input[(sequence_names[k])]
            seq2=input[(sequence_names[l])]
            score_matrix = np.zeros((seq1_length+1, seq2_length+1), dtype=int)
            traceback_matrix = np.zeros((seq1_length+1, seq2_length+1), dtype=int)
            for x in range(seq1_length+1):
                score_matrix[x][0] = x * gap_penalty
            for y in range(seq2_length+1):
                score_matrix[0][y] = y * gap_penalty
            newseq1=""
            newseq2=""


            for i in range(1,seq1_length+1):
                for j in range (1,seq2_length+1):
                    matched = score_matrix[i - 1][j - 1] + blosum62[seq1[i-1],seq2[j-1]]
                    deleted = score_matrix[i - 1][j] + gap_penalty
                    inserted = score_matrix[i][j - 1] + gap_penalty
                    score_matrix[i][j] = max(matched, deleted, inserted)



            i, j = seq1_length, seq2_length
            exact_matches=0
            allign_length=0
            while i > 0 or j > 0:
                if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + blosum62[seq1[i-1],seq2[j-1]]:
                    newseq1 = seq1[i - 1] + newseq1
                    newseq2 = seq2[j - 1] + newseq2
                    i -= 1
                    j -= 1
                elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
                    newseq1 = seq1[i - 1] + newseq1
                    newseq2 = '-' + newseq2
                    i -= 1
                else:
                    newseq1 = '-' + newseq1
                    newseq2 = seq2[j - 1] + newseq2
                    j -= 1

            for i in range(0,len(newseq1)):
                if(newseq1[i] == newseq2[i]):
                    exact_matches +=1
                allign_length += 1

            similarity_score=exact_matches/allign_length
            similarity_table[(sequence_names[k],sequence_names[l])] = f"{similarity_score:.2f}"

            #print(newseq1)
            #print(newseq2)

            #print(similarity_table)
            #print( score_matrix)

def guide_tree(similarity_table):
    entities = list(set(entity for pair in similarity_table.keys() for entity in pair))
    entity_to_index = {entity: i for i, entity in enumerate(entities)}

    num_entities = len(entities)
    distance_matrix = np.zeros((num_entities, num_entities))

    for pair, similarity_score in similarity_table.items():
        entity1, entity2 = pair
        index1, index2 = entity_to_index[entity1], entity_to_index[entity2]
        distance = 1.0 - float(similarity_score)

        distance_matrix[index1, index2] = distance
        distance_matrix[index2, index1] = distance

    # Convert the distance matrix to a Pandas DataFrame
    distance_df = pd.DataFrame(distance_matrix, index=entities, columns=entities)

    # Print the distance matrix with entity names
    print("Distance Matrix:")
    print(distance_df)



def main():
    my_hash_table = read_Blosum62()
    input1=read_fasta_file("Input.txt")
    gap_penalty = getting_gap_penalty()

    pairwise_sequence_alignment(input1,my_hash_table, gap_penalty)

    guide_tree(similarity_table)

if __name__ == "__main__":
    main()