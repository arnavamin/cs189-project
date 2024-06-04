import numpy as np

'''
The compute_alignment function takes in two nucleotide sequences and scores for matches, mismatches, and gaps to compute the
optimal global alignment graph with the given parameters.

Input:
seq1: Sequence 1
seq2: Sequence 2
match: Match score for alignment
mismatch: Mismatch score for alignment
gap: Gap score for alignment

Output:
g: An mxn matrix which holds all the computed values of the optimal global alignment
'''
def compute_alignment(seq1, seq2, match, mismatch, gap):
    # Initialize a matrix to hold the values for global alignment (dimensions are the lengths of the 2 given sequences)
    g = np.zeros((len(seq2), len(seq1)), dtype=int)

    # Compute global alignment
    for i, row in enumerate(seq2):
        for j, col in enumerate(seq1):
            if (i == 0):
                g[i][j] = j * gap   # Initialize the "gap row"
            elif (j == 0):
                g[i][j] = i * gap   # Initialize the "gap column"
            else:
                row_score = g[i][j-1] + gap # Compute the score if there is a gap on the column sequence
                col_score = g[i-1][j] + gap # Compute the score if there is a gap on the row sequence
                
                # Compute the score if there is a match/mismatch
                diagonal_score = g[i-1][j-1]    
                if (row == col):
                    diagonal_score += match 
                else:
                    diagonal_score += mismatch

                # Choose the max score from the three options and set it to the current cell in the matrix
                g[i][j] = max(row_score, col_score, diagonal_score)

    return g

'''
The get_alignment function takes the global alignment matrix and sequences and then traces back the optimal path taken by
the algorithm. It also links the optimal path to the alignment sequences.

Input:
g: mxn matrix which holds all the computed values of the optimal global alignment
seq1: sequence 1
seq2: sequence 2

Output:
path: The optimal global alignment path (max values at each step)
indices: The matrix indices of the optimal path (for tkinter application)
aligned_seq1: The aligned version of sequence 1 (adjusting for gaps)
aligned_seq2: The aligned version of sequence 2 (adjusting for gaps)

'''
def get_alignment(g, seq1, seq2):
    # Trace back to get optimal path
    path = []
    indices = []
    scores = [0, 0, 0] # Initialize an array to hold the computed scores while backtracking

    # Initialize the aligned sequences with the last nucleotide in the alignment
    # aligned_seq1 = [seq1[-1]]
    # aligned_seq2 = [seq2[-1]]
    aligned_seq1 = []
    aligned_seq2 = []

    # Set up counters for the indexes
    i = len(seq2)-1
    j = len(seq1)-1
    
    path.append(g[i][j]) # Append the final score to the path
    
    while (i > 0 or j > 0):
        scores[0] = g[i][j-1] # Backtracked score from the same row
        scores[1] = g[i-1][j] # Backtracked score from the same column
        scores[2] = g[i-1][j-1] # Backtracked score from the diagonal
        if (seq2[i] == seq1[j]):
            path.append(scores[2])
            if (i != 0 and j != 0):
                aligned_seq1.append(seq1[j])
                aligned_seq2.append(seq2[i])
                indices.append((i,j))
            i -= 1
            j -= 1
        else:
            path.append(max(scores))    # Append the max score to the path
        
            # Check to see which backtracked score was used and then append the associated nucleotide or gap to the aligned sequences
            if np.argmax(scores) == 0:
                if (len(aligned_seq1) == 0 and len(aligned_seq2) == 0):
                    aligned_seq1.append(seq1[i])
                    aligned_seq2.append('-')
                    indices.append((i,j))
                    j -= 1
                else:
                    aligned_seq1.append(seq1[j])
                    aligned_seq2.append('-')
                    indices.append((i,j))
                    j -= 1
            elif np.argmax(scores) == 1:
                if (len(aligned_seq1) == 0 and len(aligned_seq2) == 0):
                    aligned_seq1.append('-')
                    aligned_seq2.append(seq2[j])
                    indices.append((i,j))
                    i -= 1
                else:
                    aligned_seq1.append('-')
                    aligned_seq2.append(seq2[i])
                    indices.append((i,j))
                    i -= 1
            elif np.argmax(scores) == 2:
                if (i != 1 and j != 1):
                    aligned_seq1.append(seq1[j])
                    aligned_seq2.append(seq2[i])
                    indices.append((i,j))
                i -= 1
                j -= 1    

    # Reverse the path and aligned sequences to start from the beginning
    path.reverse()
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    indices.reverse()
    return path, indices, aligned_seq1, aligned_seq2

'''
Outputs the global alignment and its details.

Input:
g: mxn matrix which holds all the computed values of the optimal global alignment
path: The optimal global alignment path (max values at each step)
seq1: sequence 1
seq2: sequence 2
aligned_seq1: The aligned version of sequence 1 (adjusting for gaps)
aligned_seq2: The aligned version of sequence 2 (adjusting for gaps)
'''
def print_alignment(g, path, seq1, seq2, aligned_seq1, aligned_seq2):
    # Print the global alignment matrix
    print('Global Alignment Matrix:\n')
    print_matrix(g, seq1, seq2)
    print()

    # Print the optimal path (scores)
    print('Optimal path:', path)
    print()
    
    # Print aligned sequences
    print('Aligned Sequences:')

    for n in aligned_seq1:
        print(n, end=' ')
    print()

    alignments = ''
    for i in range(len(aligned_seq1)):
        if (aligned_seq1[i] == aligned_seq2[i]):
            alignments += '| '
        elif (aligned_seq1[i] == '-' and aligned_seq2[i] != '-'):
            alignments += '  '
        elif (aligned_seq2[i] == '-' and aligned_seq1[i] != '-'):
            alignments += '  '
        elif (aligned_seq1[i] != aligned_seq2[i]):
            alignments += '. '

    print(alignments)

    for n in aligned_seq2:
        print(n, end=' ')
    print()

def print_matrix(g, seq1, seq2):
    print('  |', end='')
    for n in seq1:
        print(' {} |'.format(n), end='')
    print()
    print('--' + '+---' * (g.shape[0]) + '+')
    for i in range(g.shape[1]):
        print(seq2[i] + ' ', end='')
        for n in g[i]:
            if n > 9 or n < -9:
                print('|{}'.format(n), end='')
            elif n < 0:
                print('|{} '.format(n), end='')
            else:
                print('| {} '.format(n), end='')
        print('|', end='')
        print('\n--' + '+---' * (g.shape[0]) + '+')

'''
Format the input sequence from a string to a list with padding.

Input:
seq: String sequence input.

Output:
formatted_seq: The formatted sequence list.
'''
def format_sequence(seq):
    formatted_seq = list(seq)
    formatted_seq.insert(0, '-')
    return formatted_seq

if __name__ == '__main__':

    #Initialize the sequences and scoring weights (will be changed for user input later)
    seq1 = 'GATTACATATACG'
    seq2 = 'GTCGACGCTACGT'

    seq1 = format_sequence(seq1)
    seq2 = format_sequence(seq2)

    match = 2
    mismatch = -1
    gap = -2

    # Compute the global alignment matrix
    g = compute_alignment(seq1, seq2, match, mismatch, gap)

    # Get the global alignment from the matrix
    path, indices, aligned_seq1, aligned_seq2 = get_alignment(g, seq1, seq2)

    # Print the global alignment and its relevant information
    print_alignment(g, path, seq1, seq2, aligned_seq1, aligned_seq2)
