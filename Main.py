import numpy as np
from scipy.cluster._optimal_leaf_ordering import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import pandas as pd
import json


similarity_table = {}

class GuideTreeNode:
    def __init__(self, name):
        self.name = name
        self.left = None
        self.right = None
        self.distance = 0.0

        
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

            print(allign_length)
            print(exact_matches)
            similarity_score=exact_matches/allign_length
            similarity_table[(sequence_names[k],sequence_names[l])] = f"{similarity_score:.2f}"

            print(newseq1)
            print(newseq2)

            print(similarity_table)

            print( score_matrix)

    return similarity_table

order_list = ["temp"]

def order_list_names(row,column):
    new_elements = []

    for k in order_list:
        if k not in row:
            new_elements.append(row)
        if k not in column:
            new_elements.append(column)

    order_list.extend(new_elements)

    print(order_list)

def guide_tree(similarity_table):
    entities = list(set(entity for pair in similarity_table.keys() for entity in pair))
    entity_to_index = {entity: i for i, entity in enumerate(entities)}

    num_entities = len(entities)
    distance_matrix = np.ones((num_entities, num_entities))

    for pair, similarity_score in similarity_table.items():
        entity1, entity2 = pair
        index1, index2 = entity_to_index[entity1], entity_to_index[entity2]
        distance = 1.0 - float(similarity_score)

        distance_matrix[index1, index2] = distance
        distance_matrix[index2, index1] = distance

    # Convert the distance matrix to a Pandas DataFrame
    distance_df = pd.DataFrame(distance_matrix, index=entities, columns=entities)
    mergedtable = distance_df
    count = 1
    merged_named_list = []
    while not mergedtable.columns.size < 1 :  # Continue until the DataFrame is empty

        min_index = mergedtable.stack().idxmin()
        row, column = min_index
        print("Minimum Distance Names:", row, column)
        print("Guide Tree: ")
        print(mergedtable)
        order_list_names(row,column)

        for i in range(0, len(mergedtable)):
            new_entity = row + column
            print()
            if mergedtable.columns[i] != row and mergedtable.columns[i] != column:

                temp = ((mergedtable.iloc[i][row] + mergedtable.iloc[i][column])/2).round(3)
                print(mergedtable.columns[i], mergedtable.iloc[i][row].round(3), "/2", "+",
                      mergedtable.iloc[i][column].round(3), "/2", "=", temp)

                mergedtable.at[new_entity, mergedtable.columns[i]] = temp
                mergedtable.at[mergedtable.columns[i], new_entity] = temp
                mergedtable.at[new_entity, new_entity] = 1
        
        merged_named = str(count), '. Step: Merged', row, 'and', column, 'to', new_entity
        merged_named_list.append(merged_named)

        mergedtable = mergedtable.drop([row, column], axis=1)
        mergedtable = mergedtable.drop([row, column], axis=0)
        count += 1


    merged_named_list = [' '.join(map(str, x)) for x in merged_named_list]
    # print every list element in a new line
    for i in range(0, len(merged_named_list)):
        print(merged_named_list[i])

    return mergedtable

# use upgma method to my distance dict to  matrix  using scipy
def UPGMA(similarity_table):
    entities = list(set(entity for pair in similarity_table.keys() for entity in pair))
    entity_to_index = {entity: i for i, entity in enumerate(entities)}

    num_entities = len(entities)
    distance_matrix = np.ones((num_entities, num_entities))

    for pair, similarity_score in similarity_table.items():
        entity1, entity2 = pair
        index1, index2 = entity_to_index[entity1], entity_to_index[entity2]
        distance = 1.0 - float(similarity_score)

        distance_matrix[index1, index2] = distance
        distance_matrix[index2, index1] = distance

    # Convert the distance matrix to a Pandas DataFrame
    distance_df = pd.DataFrame(distance_matrix, index=entities, columns=entities)
    print(distance_df)
    # Convert the DataFrame to a condensed distance matrix
    condensed_distance_matrix = squareform(distance_df, checks=False)

    # Calculate the linkage matrix
    linkage_matrix = linkage(condensed_distance_matrix, "average")

    # Convert the linkage matrix to a dendrogram
    dendrogram(linkage_matrix, labels=distance_df.index, orientation="left", color_threshold=0, above_threshold_color='black')

    # Show the plotted dendrogram
    plt.show()

# multiple sequence alignment

    

def main():
    my_hash_table = read_Blosum62()

    input1=read_fasta_file("Input.txt")
    gap_penalty = getting_gap_penalty()
    
    similarity_table = pairwise_sequence_alignment(input1,my_hash_table, gap_penalty)

    alper=guide_tree(similarity_table)






if __name__ == "__main__":
    main()