import tkinter as tk
from tkinter import ttk
import numpy as np

import global_alignment as ga
import local_alignment as la

# Function to display the alignment score matrix
def show_alignment(alignment_type):
    original_seq1 = entry_seq1.get()
    original_seq2 = entry_seq2.get()
    seq1 = ga.format_sequence(entry_seq1.get())
    seq2 = ga.format_sequence(entry_seq2.get())
    
    aligned_frame = ttk.Frame(root)
    aligned_frame.grid(row=9, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')

    try:
        match = int(entry_match.get())
        mismatch = int(entry_mismatch.get())
        gap = int(entry_gap.get())
    except ValueError:
        # Clear previous content in result_frame
        for widget in result_frame.winfo_children():
            widget.destroy()
        for widget in aligned_frame.winfo_children():
            widget.destroy()
        ttk.Label(root, text="Invalid score inputted. Please try again.").grid(row=8, column=0, padx=10, pady=10, sticky='e')
        return
    
    # For testing purposes
    # seq1 = ga.format_sequence('GATTACATATACG')
    # seq2 = ga.format_sequence('GTCGACGCTACGT')
    # match = 2
    # mismatch = -1
    # gap = -2
    if alignment_type == 'global':
        matrix = ga.compute_alignment(seq1, seq2, match, mismatch, gap)
        path, indices, aligned_seq1, aligned_seq2 = ga.get_alignment(matrix, seq1, seq2)
    elif alignment_type == 'local':
        matrix = la.compute_alignment(original_seq1, original_seq2, match, mismatch, gap)
        path, indices, aligned_seq1, aligned_seq2 = la.get_alignment(matrix, original_seq1, original_seq2)

    # Clear previous content in result_frame
    for widget in result_frame.winfo_children():
        widget.destroy()
    for widget in aligned_frame.winfo_children():
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

    

    ttk.Label(root, text="Aligned Sequences:").grid(row=8, column=0, padx=10, pady=10, sticky='e')
    # Clear previous content in aligned_frame
    for widget in aligned_frame.winfo_children():
        widget.destroy()
    if (type(aligned_seq1[0]) == list):
        for count, lst in enumerate(aligned_seq1):
            for j, char in enumerate(lst):
                label = ttk.Label(aligned_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
                label.grid(row=(count*2), column=j+1, sticky='nsew')
        for count, lst in enumerate(aligned_seq2):
            for j, char in enumerate(lst):
                label = ttk.Label(aligned_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
                label.grid(row=(count*2)+1, column=j+1, sticky='nsew')
    else:
        for j, char in enumerate(aligned_seq1):
            label = ttk.Label(aligned_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
            label.grid(row=0, column=j+1, sticky='nsew')
        
        for j, char in enumerate(aligned_seq2):
            label = ttk.Label(aligned_frame, text=char, borderwidth=1, relief="solid", width=5, anchor='center')
            label.grid(row=1, column=j+1, sticky='nsew')
        

# Create the main window
root = tk.Tk()
root.title("Alignment App")

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
ttk.Button(root, text="Global Align", command=lambda: show_alignment('global')).grid(row=5, column=0, columnspan=2, padx=10, pady=10)
ttk.Button(root, text="Local Align", command=lambda: show_alignment('local')).grid(row=6, column=0, columnspan=2, padx=10, pady=10)

# Frame to display the result matrix
result_frame = ttk.Frame(root)
result_frame.grid(row=7, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')

root.mainloop()
