import tkinter as tk
from tkinter import ttk
import numpy as np

import global_alignment as ga

# Function to display the alignment score matrix
def show_alignment():
    seq1 = ga.format_sequence(entry_seq1.get())
    seq2 = ga.format_sequence(entry_seq2.get())
    match = int(entry_match.get())
    mismatch = int(entry_mismatch.get())
    gap = int(entry_gap.get())
    
    # For testing purposes
    # seq1 = ga.format_sequence('GATTACATATACG')
    # seq2 = ga.format_sequence('GTCGACGCTACGT')
    # match = 2
    # mismatch = -1
    # gap = -2

    matrix = ga.compute_alignment(seq1, seq2, match, mismatch, gap)
    path, indices, aligned_seq1, aligned_seq2 = ga.get_alignment(matrix, seq1, seq2)

    # Clear previous content in result_frame
    for widget in result_frame.winfo_children():
        widget.destroy()

    # Create labels for column headers (seq1 characters)
    for j, char in enumerate(seq1):
        label = ttk.Label(result_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
        label.grid(row=0, column=j+1, sticky='nsew')

    # Create labels for row headers (seq2 characters)
    for i, char in enumerate(seq2):
        label = ttk.Label(result_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
        label.grid(row=i+1, column=0, sticky='nsew')

    # Create labels for the matrix values
    temp_label = tk.Label(result_frame)
    default_bg = temp_label.cget("bg")
    temp_label.destroy()
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            bg = 'blue' if (i, j) in indices else default_bg
            label = tk.Label(result_frame, text=str(int(matrix[i][j])), borderwidth=1, relief="solid", width=5, anchor='center', bg=bg)
            label.grid(row=i+1, column=j+1, sticky='nsew')

    # Configure frame to expand with window
    for i in range(len(matrix) + 1):
        result_frame.rowconfigure(i, weight=1)
    for j in range(len(matrix[0]) + 1):
        result_frame.columnconfigure(j, weight=1)

    

    ttk.Label(root, text="Aligned Sequences:").grid(row=7, column=0, padx=10, pady=10, sticky='e')
    aligned_frame = ttk.Frame(root)
    aligned_frame.grid(row=8, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    # Clear previous content in aligned_frame
    for widget in aligned_frame.winfo_children():
        widget.destroy()
    
    for j, char in enumerate(aligned_seq1):
        label = ttk.Label(aligned_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
        label.grid(row=0, column=j+1, sticky='nsew')
    
    for j, char in enumerate(aligned_seq2):
        label = ttk.Label(aligned_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
        label.grid(row=1, column=j+1, sticky='nsew')
        

# Create the main window
root = tk.Tk()
root.title("Global Alignment")

# Configure the main window grid to expand
root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)
root.grid_rowconfigure(6, weight=1)

# Input fields for sequences
ttk.Label(root, text="Sequence 1:").grid(row=0, column=0, padx=5, pady=5, sticky='e')
entry_seq1 = ttk.Entry(root)
entry_seq1.grid(row=0, column=1, padx=5, pady=5, sticky='ew')

ttk.Label(root, text="Sequence 2:").grid(row=1, column=0, padx=5, pady=5, sticky='e')
entry_seq2 = ttk.Entry(root)
entry_seq2.grid(row=1, column=1, padx=5, pady=5, sticky='ew')

# Input fields for scoring parameters
ttk.Label(root, text="Match Score:").grid(row=2, column=0, padx=5, pady=5, sticky='e')
entry_match = ttk.Entry(root)
entry_match.grid(row=2, column=1, padx=5, pady=5, sticky='ew')

ttk.Label(root, text="Mismatch Score:").grid(row=3, column=0, padx=5, pady=5, sticky='e')
entry_mismatch = ttk.Entry(root)
entry_mismatch.grid(row=3, column=1, padx=5, pady=5, sticky='ew')

ttk.Label(root, text="Gap Score:").grid(row=4, column=0, padx=5, pady=5, sticky='e')
entry_gap = ttk.Entry(root)
entry_gap.grid(row=4, column=1, padx=5, pady=5, sticky='ew')

# Button to trigger the alignment calculation
ttk.Button(root, text="Align", command=show_alignment).grid(row=5, column=0, columnspan=2, pady=10)

# Frame to display the result matrix
result_frame = ttk.Frame(root)
result_frame.grid(row=6, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')

root.mainloop()
