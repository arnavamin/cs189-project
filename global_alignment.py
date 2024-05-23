import numpy as np

def run():
    seq1 = ['-', 'G', 'A', 'C', 'A']
    seq2 = ['-', 'A', 'C', 'A']
    path = []
    
    g = np.zeros((len(seq2), len(seq1)), dtype=int)

    # Alignment parameters
    match = 2
    mismatch = -1
    gap = -2

    # Compute global alignment
    for i, row in enumerate(seq2):
        for j, col in enumerate(seq1):
            if (i == 0):
                g[i][j] = j * gap
            elif (j == 0):
                g[i][j] = i * gap
            else:
                row_score = g[i][j-1] + gap
                col_score = g[i-1][j] + gap
                diagonal_score = g[i-1][j-1]
                
                if (row == col):
                    diagonal_score += match
                else:
                    diagonal_score += mismatch

                g[i][j] = max(row_score, col_score, diagonal_score)

    # Trace back to get optimal path
    aligned_seq1 = [seq1[-1]]
    aligned_seq2 = [seq2[-1]]
    seq1_counter = 1
    seq2_counter = 1

    i = len(seq2)-1
    j = len(seq1)-1
    path.append(g[i][j])
    scores = [0, 0, 0]
    while (i != 0 and j != 0):
        scores[0] = g[i][j-1]
        scores[1] = g[i-1][j]
        scores[2] = g[i-1][j-1]
        path.append(max(scores))
        if np.argmax(scores) == 0:
            j -= 1
            seq2_counter += 1
            aligned_seq1.append('-')
            if (seq2_counter >= len(seq2)):
                aligned_seq2.append('-')
            else:
                aligned_seq2.append(seq2[-seq2_counter])
        elif np.argmax(scores) == 1:
            i -= 1
            seq1_counter += 1
            if (seq1_counter >= len(seq1)):
                aligned_seq1.append('-')
            else:
                aligned_seq1.append(seq1[-seq1_counter])
            aligned_seq2.append('-')
        elif np.argmax(scores) == 2:
            i -= 1
            j -= 1
            seq1_counter += 1
            seq2_counter += 1
            if (seq1_counter >= len(seq1)):
                aligned_seq1.append('-')
            else:
                aligned_seq1.append(seq1[-seq1_counter])
            if (seq2_counter >= len(seq2)):
                aligned_seq2.append('-')
            else:
                aligned_seq2.append(seq2[-seq2_counter])
            
    
    
        


    print(g)
    path.reverse()
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    print('Optimal path:', path)
    print(aligned_seq1)
    print(aligned_seq2)

if __name__ == '__main__':
    run()
