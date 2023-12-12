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

def main():

    my_hash_table = read_Blosum62()
    #gap_penalty= getting_gap_penalty()
    input1=read_fasta_file("Input2.txt")

if __name__ == "__main__":
    main()