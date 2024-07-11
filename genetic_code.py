import tkinter as tk
from tkinter import messagebox, Toplevel, Text, filedialog 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
from Bio import pairwise2, SeqIO, Phylo
from collections import Counter
from PIL import Image, ImageTk  
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio import pairwise2, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np

# Function to calculate GC content of a DNA sequence
def calculate_gc_content():
    dna_sequence = entry_dna.get()  # Get the DNA sequence from the entry 
    if dna_sequence:  # Check if a sequence is provided
        seq = Seq(dna_sequence)  # Create a Seq object
        gc_content = seq.count("G") + seq.count("C")  # Count G and C bases
        total_bases = len(seq)  # Get the total number of bases
        gc_percentage = (gc_content / total_bases) * 100  # Calculate GC percentage
        messagebox.showinfo("GC Content", f"The GC content is: {gc_percentage:.2f}%")  # Show the result
    else:
        messagebox.showerror("Error", "Please load a DNA sequence from a FASTA file.")  # Show an error message if no sequence is provided

# Function to find the reverse complement of a DNA sequence
def find_reverse_complement():
    dna_sequence = entry_dna.get()  
    if dna_sequence:  
        seq = Seq(dna_sequence)  
        reverse_complement_seq = seq.reverse_complement()  # Get the reverse complement
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])  # Open a save dialog
        if file_path:  # Check if a file path is provided
            with open(file_path, "w") as file:  # Save the reverse complement to the file
                file.write(f"The reverse complement sequence is: {reverse_complement_seq}")
                messagebox.showinfo("Reverse Complement", f"Reverse complement sequence saved to {file_path}")  # Show success message
    else:
        messagebox.showerror("Error", "Please load a DNA sequence from a FASTA file.")  # Show an error message if no sequence is provided

# Function to find codon usage in a DNA sequence
def find_codon_usage():
    dna_sequence = entry_dna.get()  
    if dna_sequence:  
        seq = Seq(dna_sequence)  
        codon_counts = Counter(seq[i:i+3] for i in range(0, len(seq)-2, 3))  # Count codons in the sequence
        codon_usage_str = '\n'.join([f'{codon}: {count}' for codon, count in codon_counts.items()])  # Create a string of codon counts
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])  # Open a save dialog
        if file_path:  # Check if a file path is provided
            with open(file_path, "w") as file:  # Save the codon usage to the file
                file.write(f"The codon usage is:\n{codon_usage_str}")
                messagebox.showinfo("Codon Usage", f"Codon usage saved to {file_path}")  # Show success message
    else:
        messagebox.showerror("Error", "Please load a DNA sequence from a FASTA file.")  # Show an error message if no sequence is provided


# Function to open the alignment window
def open_alignment_window():
    align_window = Toplevel(root)  # Create a new window for sequence alignment
    align_window.title("Sequence Alignment") # Set the title of the new window
    
    text_sequences = Text(align_window, height=10, width=50) # Create a text widget for inputting sequences
    text_sequences.pack()

    alignment_option = tk.StringVar()  # Create a string variable for the alignment option
    alignment_option.set("pairwise")  # Set the default alignment option to pairwise

    def perform_alignment():
        sequences = [seq.strip() for seq in text_sequences.get("1.0", tk.END).strip().split('\n') if seq.strip()]  # Get the sequences from the text widget
        if len(sequences) >= 2:  # Check if at least two sequences are provided
            alignment_results = []  # List to store alignment results
            for i in range(len(sequences)):
                for j in range(i+1, len(sequences)):
                    seq1 = sequences[i]
                    seq2 = sequences[j]
                    # alignments = pairwise2.align.localxx(seq1, seq2)
                    # alignment_result = alignments[0]
                    if alignment_option.get() == "pairwise":
                        alignments = pairwise2.align.localxx(seq1, seq2)
                        alignment_result = alignments[0]
                    # elif alignment_option.get() == "multiple":
                    #     alignments = pairwise2.align.localxx([seq1, seq2])
                    #     alignment_result = alignments[-1]  # Taking the last (best) alignment
                    align1, align2, score, begin, end = alignment_result

                    identity = sum(1 for a, b in zip(align1, align2) if a == b)
                    similarity = identity
                    gaps = align1.count('-') + align2.count('-')
                    alignment_length = len(align1)

                    identity_percentage = (identity / alignment_length) * 100
                    similarity_percentage = (similarity / alignment_length) * 100
                    gap_percentage = (gaps / alignment_length) * 100

                    result = f"Score: {score}\n" \
                             f"Alignment length: {alignment_length}\n" \
                             f"Identity: {identity}/{alignment_length} ({identity_percentage:.2f}%)\n" \
                             f"Similarity: {similarity}/{alignment_length} ({similarity_percentage:.2f}%)\n" \
                             f"Gaps: {gaps}/{alignment_length} ({gap_percentage:.2f}%)\n"
                    alignment_results.append((align1, align2, result))

            # Display alignment results
            display_alignment_results(alignment_results)

            # Perform multiple sequence alignment and generate phylogenetic tree
            # print("Sequence: " ,sequences)
            generate_phylogenetic_tree(sequences)
            # print("accuracy: ",calculate_accuracy(seq1,seq2))

        else:
            messagebox.showerror("Error", "Please enter at least two sequences.")

    tk.Button(align_window, text="Perform Alignment", command=perform_alignment).pack()

def calculate_accuracy(seq1, seq2):
    # Dummy accuracy calculation for illustration
    # print("seq1: ",seq1)
    # print("seq2: ",seq2)
    return np.random.uniform(0, 1)


# Function to generate a phylogenetic tree
def generate_phylogenetic_tree(sequences):
    try:
        # Create SeqRecord objects for each sequence
        records = [SeqRecord(Seq(seq), id=f"Sequence {i+1}") for i, seq in enumerate(sequences)]

        # Create MultipleSeqAlignment object
        alignment = MultipleSeqAlignment(records)

        # Calculate distances between sequences
        calculator = DistanceCalculator('identity')
        
        #distance matrix
        dm = calculator.get_distance(alignment)

        # print(f"DM: ",dm)
        # Construct a phylogenetic tree using UPGMA method
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)

        # Color coding based on accuracy
        clade_colors = {}
        for clade in tree.get_terminals():
            clade_name = clade.name
            clade_seq = next(record.seq for record in records if record.id == clade_name)
            max_accuracy = 0
            for other_clade in tree.get_terminals():
                if other_clade == clade:
                    continue
                other_clade_seq = next(record.seq for record in records if record.id == other_clade.name)
                accuracy = calculate_accuracy(clade_seq, other_clade_seq)
                # print(f"clade seq: %s",clade_seq)
                # print(f"other clade seq: %s",other_clade_seq)
                # print(f"accuracy: %",accuracy)
                
                max_accuracy = max(max_accuracy, accuracy)

                # print(f"maximum accuracy: %",max_accuracy)
            if max_accuracy > 0.7:  # Arbitrary threshold for high accuracy
                clade_colors[clade_name] = 'red'
            elif max_accuracy > 0.4:
                clade_colors[clade_name] = 'orange'
            else:
                clade_colors[clade_name] = 'green'

        # Display the phylogenetic tree with improvements
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax, label_colors=clade_colors)

        # Customize the tree
        ax.set_title("Phylogenetic Tree of Sequences")
        ax.set_xlabel("Branch Length")
        ax.set_ylabel("Taxa")

        plt.show()
        messagebox.showinfo("Phylogenetic Tree", "Phylogenetic tree generated successfully.")
    except Exception as e:
        messagebox.showerror("Phylogenetic Tree Error", f"Error generating phylogenetic tree: {str(e)}")


def display_alignment_results(alignment_results):
    result_window = Toplevel(root)
    result_window.title("Alignment Results")

    text = Text(result_window, height=20, width=80)
    text.pack()

    for align1, align2, result in alignment_results:
        text.insert(tk.END, "Alignment:\n")
        text.insert(tk.END, align1 + "\n")
        text.insert(tk.END, align2 + "\n")
        text.insert(tk.END, result + "\n\n")

# Function to translate DNA to protein sequence
def translate_dna_to_protein():
    dna_sequence = entry_dna.get()
    if dna_sequence:
        seq = Seq(dna_sequence)
        protein_sequence = seq.translate()
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, "w") as file:
                file.write(f"The translated protein sequence is: {protein_sequence}")
                messagebox.showinfo("Protein Translation", f"Protein sequence saved to {file_path}")
    else:
        messagebox.showerror("Error", "Please load a DNA sequence from a FASTA file.")

# Function to convert DNA to mRNA sequence
def convert_dna_to_mrna():
    dna_sequence = entry_dna.get()
    if dna_sequence:
        translation_table = str.maketrans("TACG", "AUGC")
        mrna_sequence = dna_sequence.translate(translation_table)
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, "w") as file:
                file.write(f"The mRNA sequence is: {mrna_sequence}")
                messagebox.showinfo("mRNA Sequence", f"mRNA sequence saved to {file_path}")
    else:
        messagebox.showerror("Error", "Please load a DNA sequence from a FASTA file.")

# fasta file load korbe ekhne trick hcche fasta file e duita sequence ache ekta organism eri oita select korbo
def load_fasta_file():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.fna")])
    if file_path:
        with open(file_path, "r") as file:
            records = list(SeqIO.parse(file, "fasta"))
            if len(records) > 1:
                seqs_window = Toplevel(root)
                seqs_window.title("Select a DNA Sequence")
                seqs_window.configure(bg='lightgray')

                tk.Label(seqs_window, text="Multiple sequences found. Please select one:", **style).pack(pady=10)
                
                seq_var = tk.StringVar()
                seq_var.set(str(records[0].seq))
                
                for record in records:
                    tk.Radiobutton(seqs_window, text=f"{record.id}: {record.description[:50]}...", 
                                   variable=seq_var, value=str(record.seq), bg='lightgray').pack(anchor='w', padx=10, pady=5)

                def select_sequence():
                    selected_seq = seq_var.get()
                    entry_dna.delete(0, tk.END)
                    entry_dna.insert(0, selected_seq)
                    seqs_window.destroy()
                    messagebox.showinfo("FASTA File Loaded", "Your FASTA file uploaded successfully.\nYou can check your DNA sequence if you have any doubt.")
                    save_sequence(selected_seq)

                tk.Button(seqs_window, text="Select", command=select_sequence, **button_styles).pack(pady=10)
            else:
                entry_dna.delete(0, tk.END)
                entry_dna.insert(0, str(records[0].seq))
                messagebox.showinfo("FASTA File Loaded", "Your FASTA file uploaded successfully.\nYou can check your DNA sequence if you have any doubt.")
                save_sequence(str(records[0].seq))
# operation gular resultsob file e save hye thakbe
def save_sequence(sequence):
    file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
    if file_path:
        with open(file_path, "w") as file:
            file.write(sequence)
            messagebox.showinfo("Sequence Saved", f"DNA sequence saved to {file_path}")

root = tk.Tk()
root.title("Sequence Analysis Tools")
root.geometry("950x570")

# Load background image
image_path = r'C:\Users\samiu\OneDrive\Desktop\L-3,S-1\CSM 3121 & 3122 Systems and Software Engineering\software project\AdobeStock_730210234_Preview.jpeg'

bg_image = Image.open(image_path)
bg_photo = ImageTk.PhotoImage(bg_image)

# Create a canvas
canvas = tk.Canvas(root, width=700, height=800)
canvas.pack(fill="both", expand=True)

# Set the background image
canvas.create_image(0, 0, image=bg_photo, anchor="nw")

# Style dictionary
style = {"font": ("Comic Sans MS", 14), "bg": "lightblue"}
button_styles = {"font": ("Comic Sans MS", 12), "bg": "PapayaWhip","fg": "GREEN"}

# Add trick to the canvas
canvas.create_text(450, 50, text="Sequence Analysis Tools", font=("Comic Sans MS", 24), fill="RED")
entry_dna = tk.Entry(root, font=("Comic Sans MS", 12), width=60)
canvas.create_window(450, 100, window=entry_dna)

btn_load_fasta = tk.Button(root, text="Load FASTA File", command=load_fasta_file, **button_styles)
canvas.create_window(350, 150, window=btn_load_fasta)

btn_gc_content = tk.Button(root, text="Calculate GC Content", command=calculate_gc_content, **button_styles)
canvas.create_window(550, 150, window=btn_gc_content)

btn_reverse_complement = tk.Button(root, text="Find Reverse Complement", command=find_reverse_complement, **button_styles)
canvas.create_window(340, 220, window=btn_reverse_complement)

btn_codon_usage = tk.Button(root, text="Find Codon Usage", command=find_codon_usage, **button_styles)
canvas.create_window(550, 220, window=btn_codon_usage)

btn_alignment = tk.Button(root, text="Perform Alignment", command=open_alignment_window, **button_styles)
canvas.create_window(450, 280, window=btn_alignment)

btn_translate = tk.Button(root, text="Translate DNA to Protein", command=translate_dna_to_protein, **button_styles)
canvas.create_window(330, 350, window=btn_translate)

btn_convert_mrna = tk.Button(root, text="Convert DNA to mRNA", command=convert_dna_to_mrna, **button_styles)
canvas.create_window(550, 350, window=btn_convert_mrna)

root.mainloop()
