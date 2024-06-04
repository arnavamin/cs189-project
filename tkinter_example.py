import tkinter as tk
from tkinter import ttk
import numpy as np

import global_alignment as ga

# Function to display the alignment score matrix
def show_alignment():
    seq1_original = 'GATTACATATACG'
    seq2_original = 'GTCGACGCTACGT'

    seq1 = ga.format_sequence(seq1_original)
    seq2 = ga.format_sequence(seq2_original)

    match = 2
    mismatch = -1
    gap = -2

    matrix = ga.compute_alignment(seq1, seq2, match, mismatch, gap)

    # Clear previous content in result_frame
    for widget in result_frame.winfo_children():
        widget.destroy()

    # Create Treeview to display matrix
    tree = ttk.Treeview(result_frame, show='headings')
    tree.grid(row=0, column=0, sticky='nsew')

    # # Configure Treeview columns
    # print(len(seq1), len(matrix))
    print(seq1)
    tree['columns'] = seq1
    for col in tree['columns']:
        tree.heading(col, text=col)
        tree.column(col, width=40, anchor='center')
    # for i, col in enumerate(matrix[0]):
    #     col_label = '-' if i == 0 else seq1_original[i-1]
    #     tree.insert('', 'end', values=[col_label])

    # Add rows to Treeview
    for i, row in enumerate(matrix):
        row_label = '-' if i == 0 else seq2_original[i-1]
        tree.insert('', 'end', values=[row_label] + [int(val) for val in row])

    # Configure frame to expand with window
    result_frame.rowconfigure(0, weight=1)
    result_frame.columnconfigure(0, weight=1)

    # Add horizontal and vertical scrollbars
    x_scroll = ttk.Scrollbar(result_frame, orient='horizontal', command=tree.xview)
    y_scroll = ttk.Scrollbar(result_frame, orient='vertical', command=tree.yview)
    tree.configure(xscrollcommand=x_scroll.set, yscrollcommand=y_scroll.set)
    x_scroll.grid(row=1, column=0, sticky='ew')
    y_scroll.grid(row=0, column=1, sticky='ns')

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
