import numpy as np

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

def pairwise_sequence_alignment(seq1, seq2,blosum62, gap_penalty):
    seq1_length = len(seq1)
    seq2_length = len(seq2)

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



    print(newseq1)
    print(newseq2)  

    return score_matrix

def main():
    my_hash_table = read_Blosum62()
    #gap_penalty= getting_gap_penalty()
    input1=read_fasta_file("Input2.txt")
    gap_penalty = getting_gap_penalty()
    values_iterator = iter(input1.values())
    first_value = next(values_iterator, None)
    second_value = next(values_iterator, None)

    #print(pairwise_sequence_alignment(first_value,second_value,my_hash_table, gap_penalty))


    print(pairwise_sequence_alignment("CQLMKTERPRPNTFVIRCLQWTTVIERTFHVDSPDEREEWMRAIQMVANSLKQRGPGEDA","MKTERPRPNTFIIRCLQWTTVIERTFHVETPEEREEWTTAIQTVADGLKKQEEEE",my_hash_table, gap_penalty))

if __name__ == "__main__":
    main()