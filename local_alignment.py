import numpy as np

#Converts the string to list for the code.
def convert_str_to_list(original_seq1, original_seq2):
    seq1 = list(original_seq1)
    seq1.insert(0, "-")

    seq2 = list(original_seq2)
    seq2.insert(0, "-")

    return seq1, seq2

#The compute_alignment function takes in two nucleotide sequences and scores for matches, mismatches, and gaps to compute the optimal local alignment graph with the given parameters.
def compute_alignment(original_seq1, original_seq2, match, mismatch, gap):
    seq1, seq2 = convert_str_to_list(original_seq1, original_seq2)

    l = np.zeros((len(seq2), len(seq1)), dtype=int)

    for i, row in enumerate(seq2):
        for j, col in enumerate(seq1):
            if (i == 0):
                l[i][j] = j * 0  
            elif (j == 0):
                l[i][j] = i * 0 
            else:
                row_score = l[i][j-1] + gap# Compute the score if there is a gap on the row sequence
                col_score = l[i-1][j] + gap# Compute the score if there is a gap on the column sequence

                if row_score < 0: #If row score is below 0, set row score as 0. 
                    row_score = 0
                else:
                    pass

                if col_score < 0: #If column score is below 0, set column score as 0. 
                    col_score = 0
                else:
                    pass

                # Compute the score if there is a match/mismatch

                diagonal_score = l[i-1][j-1] 

                if (row == col):
                    diagonal_score += match 
                else:
                    diagonal_score += mismatch

                if diagonal_score < 0:
                    diagonal_score = 0
                else:
                    pass

                l[i][j] = max(row_score, col_score, diagonal_score)
    return l

#Prints the matrix on the backend. 
def print_matrix(l, original_seq1, original_seq2):

    seq1, seq2 = convert_str_to_list(original_seq1, original_seq2)

    print('  |', end='')
    for n in seq1:
        print(' {} |'.format(n), end='')
    print()
    print('--' + '+---' * (l.shape[0]) + '+')
    for i in range(l.shape[1]):
        print(seq2[i] + ' ', end='')
        for n in l[i]:
            if n > 9 or n < -9:
                print('|{}'.format(n), end='')
            elif n < 0:
                print('|{} '.format(n), end='')
            else:
                print('| {} '.format(n), end='')
        print('|', end='')
        print('\n--' + '+---' * (l.shape[0]) + '+')
        
#The get_alignment function takes the local alignment matrix and sequences and then traces back the optimal path taken by the algorithm. It also links the optimal path to the alignment sequences.

def get_alignment(l, original_seq1, original_seq2):
    
    seq1, seq2 = convert_str_to_list(original_seq1, original_seq2)

    index_positions = list(np.argwhere(l == np.max(l)))

    indices = []
    output_path = []
    output_align1 = []
    output_align2 = []

    for pos in index_positions:
        path = []

        scores = [0, 0, 0]

        aligned_seq1 = []
        aligned_seq2 = []

        max_value = np.max(l)

        path.append(max_value)

        x = pos[0]
        y = pos[1]

        while (l[x][y] >= 1):
            scores[0] = l[x][y-1] # Backtracked score from the same row
            scores[1] = l[x-1][y] # Backtracked score from the same column
            scores[2] = l[x-1][y-1] # Backtracked score from the diagonal

            if (seq2[x] == seq1[y]):
                path.append(scores[2])
                if (x != 0 and y != 0):
                    indices.append((x, y))
                    aligned_seq1.append(seq1[y])
                    aligned_seq2.append(seq2[x])
                x -= 1
                y -= 1
            else:
                path.append(max(scores))

                if np.argmax(scores) == 0:
                    if (len(aligned_seq1) == 0 and len(aligned_seq2) == 0):
                        aligned_seq1.append(seq1[x])
                        aligned_seq2.append('-')
                        indices.append((x, y))
                        y -= 1
                    else:
                        aligned_seq1.append(seq1[x+1])
                        aligned_seq2.append('-')
                        indices.append((x, y))
                        y -= 1
                elif np.argmax(scores) == 1:
                    if (len(aligned_seq1) == 0 and len(aligned_seq2) == 0):
                        aligned_seq1.append('-')
                        aligned_seq2.append(seq2[y])
                        indices.append((x, y))
                        x -= 1
                    else:
                        aligned_seq1.append('-')
                        aligned_seq2.append(seq2[y])
                        indices.append((x, y))
                        x -= 1
                elif np.argmax(scores) == 2:
                    if (x != 1 and y != 1):
                        aligned_seq1.append(seq1[y])
                        aligned_seq2.append(seq2[x])
                        indices.append((x, y))
                    x -= 1
                    y -= 1   


        path.reverse()
        path.remove(0)
        output_path.append(path)
        aligned_seq1.reverse()
        output_align1.append(aligned_seq1)
        aligned_seq2.reverse()
        output_align2.append(aligned_seq2)
        indices.reverse()

    return output_path, indices, output_align1, output_align2

#print_alignment outputs the local alignment and its details. 

def print_alignment(l, output_path, output_align1, output_align2, original_seq1, original_seq2):
    print()
    print('Local Alignment Matrix:\n')
    print_matrix(l, original_seq1, original_seq2)
    print()

    for num in range(len(output_path)):
        if len(output_path) == 1:
            print("Only Optimal Path")
        elif len(output_path) > 1:
            print("Optimal Path", num+1)
        print(output_path[num])

    print()

    for m in range(len(output_align1)):
        if len(output_align1) == 1:
            print("Only Possibility")
        elif len(output_align1) > 1:
            print("Possibility", m+1)
        for n in output_align1[m]:
            print(n, end=' ')
        print()
        alignments = ''
        for s in range(len(output_align1[m])):
            if (output_align1[m][s] == output_align2[m][s]):
                alignments += '| '
            elif (output_align1[m][s] == '-' and output_align2[m][s] != '-'):
                alignments += '  '
            elif (output_align1[m][s] == '-' and output_align2[m][s] != '-'):
                alignments += '  '
            elif (output_align1[m][s] != output_align2[m][s]):
                alignments += '  '

        print(alignments)

        for n in output_align2[m]:
            print(n, end=' ')
        print()


if __name__ == '__main__':

    original_seq1 = 'GATTACATATACG'
    original_seq2 = "GTCGACGCTACGT"

    match = 2
    mismatch = -1
    gap = -2

    l = compute_alignment(original_seq1, original_seq2, match, mismatch, gap)

    output_path, indices, output_align1, output_align2 = get_alignment(l, original_seq1, original_seq2)

    print_alignment(l, output_path, output_align1, output_align2, original_seq1, original_seq2)
