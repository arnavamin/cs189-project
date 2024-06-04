import tkinter as tk
from tkinter import ttk
import numpy as np
import io

# Provided alignment functions
def convert_str_to_list(original_seq1, original_seq2):
    seq1 = list(original_seq1)
    seq1.insert(0, "-")
    seq2 = list(original_seq2)
    seq2.insert(0, "-")
    return seq1, seq2

def compute_alignment(original_seq1, original_seq2, match, mismatch, gap):
    seq1, seq2 = convert_str_to_list(original_seq1, original_seq2)
    l = np.zeros((len(seq2), len(seq1)), dtype=int)
    for i, row in enumerate(seq2):
        for j, col in enumerate(seq1):
            if i == 0:
                l[i][j] = j * 0  
            elif j == 0:
                l[i][j] = i * 0 
            else:
                row_score = l[i][j-1] + gap
                col_score = l[i-1][j] + gap
                if row_score < 0:
                    row_score = 0
                if col_score < 0:
                    col_score = 0
                diagonal_score = l[i-1][j-1] 
                if row == col:
                    diagonal_score += match 
                else:
                    diagonal_score += mismatch
                if diagonal_score < 0:
                    diagonal_score = 0
                l[i][j] = max(row_score, col_score, diagonal_score)
    return l

def get_alignment(l, original_seq1, original_seq2):
    seq1, seq2 = convert_str_to_list(original_seq1, original_seq2)
    index_positions = list(np.argwhere(l == np.max(l)))
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
        while l[x][y] >= 1:
            scores[0] = l[x][y-1] 
            scores[1] = l[x-1][y] 
            scores[2] = l[x-1][y-1] 
            if seq2[x] == seq1[y]:
                path.append(scores[2])
                if x != 0 and y != 0:
                    aligned_seq1.append(seq1[y])
                    aligned_seq2.append(seq2[x])
                x -= 1
                y -= 1
            else:
                path.append(max(scores))
                if np.argmax(scores) == 0:
                    aligned_seq1.append(seq1[y])
                    aligned_seq2.append('-')
                    y -= 1
                elif np.argmax(scores) == 1:
                    aligned_seq1.append('-')
                    aligned_seq2.append(seq2[x])
                    x -= 1
                elif np.argmax(scores) == 2:
                    if x != 1 and y != 1:
                        aligned_seq1.append(seq1[y])
                        aligned_seq2.append(seq2[x])
                    x -= 1
                    y -= 1   
        path.reverse()
        path.remove(0)
        output_path.append(path)
        aligned_seq1.reverse()
        output_align1.append(aligned_seq1)
        aligned_seq2.reverse()
        output_align2.append(aligned_seq2)
    return output_path, output_align1, output_align2

def print_matrix(l, original_seq1, original_seq2):
    seq1, seq2 = convert_str_to_list(original_seq1, original_seq2)
    output = io.StringIO()
    output.write('  |')
    for n in seq1:
        output.write(' {} |'.format(n))
    output.write('\n')
    output.write('--' + '+---' * len(seq1) + '+\n')
    for i in range(len(seq2)):
        output.write(seq2[i] + ' ')
        for n in l[i]:
            if n > 9 or n < -9:
                output.write('|{}'.format(n))
            elif n < 0:
                output.write('|{} '.format(n))
            else:
                output.write('| {} '.format(n))
        output.write('|\n')
        output.write('--' + '+---' * len(seq1) + '+\n')
    return output.getvalue()

def print_alignment(l, output_path, output_align1, output_align2, original_seq1, original_seq2):
    output = io.StringIO()
    output.write('\nLocal Alignment Matrix:\n\n')
    output.write(print_matrix(l, original_seq1, original_seq2))
    output.write('\n')
    for num in range(len(output_path)):
        if len(output_path) == 1:
            output.write("Only Optimal Path\n")
        elif len(output_path) > 1:
            output.write("Optimal Path {}\n".format(num+1))
        output.write(str(output_path[num]) + '\n')
    output.write('\n')
    for m in range(len(output_align1)):
        if len(output_align1) == 1:
            output.write("Only Possibility\n")
        elif len(output_align1) > 1:
            output.write("Possibility {}\n".format(m+1))
        for n in output_align1[m]:
            output.write(n + ' ')
        output.write('\n')
        alignments = ''
        for s in range(len(output_align1[m])):
            if output_align1[m][s] == output_align2[m][s]:
                alignments += '| '
            else:
                alignments += '  '
        output.write(alignments + '\n')
        for n in output_align2[m]:
            output.write(n + ' ')
        output.write('\n\n')
    return output.getvalue()

# GUI implementation
class AlignmentApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Local Alignment")
        self.geometry("800x600")
        self.create_widgets()
    
    def create_widgets(self):
        # Input frame
        input_frame = ttk.LabelFrame(self, text="Input Sequences")
        input_frame.pack(padx=10, pady=10, fill="x", expand="yes")

        ttk.Label(input_frame, text="Sequence 1:").grid(row=0, column=0, padx=5, pady=5)
        self.seq1_entry = ttk.Entry(input_frame, width=40)
        self.seq1_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="Sequence 2:").grid(row=1, column=0, padx=5, pady=5)
        self.seq2_entry = ttk.Entry(input_frame, width=40)
        self.seq2_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="Match Score:").grid(row=2, column=0, padx=5, pady=5)
        self.match_entry = ttk.Entry(input_frame, width=10)
        self.match_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(input_frame, text="Mismatch Score:").grid(row=3, column=0, padx=5, pady=5)
        self.mismatch_entry = ttk.Entry(input_frame, width=10)
        self.mismatch_entry.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(input_frame, text="Gap Penalty:").grid(row=4, column=0, padx=5, pady=5)
        self.gap_entry = ttk.Entry(input_frame, width=10)
        self.gap_entry.grid(row=4, column=1, padx=5, pady=5, sticky="w")

        ttk.Button(input_frame, text="Compute Alignment", command=self.compute_alignment).grid(row=5, column=0, columnspan=2, pady=10)

        # Output frame
        output_frame = ttk.LabelFrame(self, text="Results")
        output_frame.pack(padx=10, pady=10, fill="both", expand="yes")

        self.result_text = tk.Text(output_frame, wrap="word", height=20)
        self.result_text.pack(padx=10, pady=10, fill="both", expand="yes")

    def compute_alignment(self):
        seq1 = self.seq1_entry.get()
        seq2 = self.seq2_entry.get()
        match = int(self.match_entry.get())
        mismatch = int(self.mismatch_entry.get())
        gap = int(self.gap_entry.get())

        l = compute_alignment(seq1, seq2, match, mismatch, gap)
        output_path, output_align1, output_align2 = get_alignment(l, seq1, seq2)

        self.result_text.delete("1.0", tk.END)
        alignment_output = print_alignment(l, output_path, output_align1, output_align2, seq1, seq2)
        self.result_text.insert(tk.END, alignment_output)

if __name__ == "__main__":
    app = AlignmentApp()
    app.mainloop()
