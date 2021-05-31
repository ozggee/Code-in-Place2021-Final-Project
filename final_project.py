from tkinter import *
def main():
    root = Tk()
    root.title('DNA ANALYZER')
    canvas = Canvas(root, height=650, width=950)
    canvas.pack()
    # Frame
    frame_up = Frame(root, bg='powder blue')
    frame_up.place(relx=0.1, rely=0.08, relwidth=0.8, relheight=0.15)
    frame_down_left = Frame(root, bg='powder blue')
    frame_down_left.place(relx=0.1, rely=0.25, relwidth=0.23, relheight=0.66)
    frame_down_right = Frame(root, bg='powder blue')
    frame_down_right.place(relx=0.34, rely=0.25, relwidth=0.56, relheight=0.66)
    # label of dna sequence and dna input
    enter_label = Label(frame_up, bg='powder blue', text='DNA sequence:', font='Verdana 14')
    enter_label.pack(padx=10, pady=10, side=LEFT)
    seq = Entry(frame_up, bg='azure', font='Verdana 15', width=40)
    seq.place(x=230, y=27)
    dna = seq.get().upper()
    #functions
    len_dna(dna)
    counting(dna)
    gc(dna)
    translate_rna(dna)
    reversing_complement(dna)
    protein(dna)
    # Control buttons
    button1 = Button(frame_down_left, text='length of DNA', padx=80, pady=24, bg='medium sea green', command=len_dna(dna))
    button2 = Button(frame_down_left, text='counting bases', padx=68, pady=24, bg='medium sea green', command=counting(dna))
    button3 = Button(frame_down_left, text='gc-content', padx=77, pady=24, bg='medium sea green', command=gc(dna))
    button4 = Button(frame_down_left, text='reverse complement', padx=70, pady=24, bg='medium sea green',
                     command=reversing_complement(dna))
    button5 = Button(frame_down_left, text='RNA sequence', padx=70, pady=24, bg='medium sea green',
                     command=translate_rna(dna))
    button6 = Button(frame_down_left, text='protein sequence', padx=105, pady=24, bg='medium sea green',
                     command=protein(dna))
    button1.pack()
    button2.pack()
    button3.pack()
    button4.pack()
    button5.pack()
    button6.pack()
    root.mainloop()

def len_dna(dna):
    seq1 = len(dna)
    label_dna = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Length of DNA sequence: {seq1}')
    label_dna.pack(pady=7)

def counting(dna):
    a = 0
    t = 0
    g = 0
    c = 0
    for base in dna:
        if base == 'A':
            a = a + 1
        elif base == 'T':
            t = t + 1
        elif base == 'G':
            g = g + 1
        elif base == 'C':
            c = c + 1
    label_adenine = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Adenine: {a}')
    label_thymine = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Thymine: {t}')
    label_guanine = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Guanine: {g}')
    label_cytosine = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Cytosine: {c}')
    label_adenine.pack(pady=7)
    label_thymine.pack(pady=7)
    label_guanine.pack(pady=7)
    label_cytosine.pack(pady=7)

def gc(dna):
    g = 0
    c = 0
    for base in dna:
        if base == 'G':
            g = g + 1
        elif base == 'C':
            c = c + 1
    gc_content = (g + c) / len(dna)
    label_gc = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'GC-content: {gc_content}')
    label_gc.pack(pady=7)

def translate_rna(dna):
    rna = dna.replace('T', 'U')
    label_rna = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'RNA sequence: {rna}')
    label_rna.pack(pady=7)

def reversing_complement(dna):
    reverse_comp = dna[::-1]
    at_comp = reverse_comp.replace('A', 't').replace('T', 'a')
    gc_comp = at_comp.replace('G', 'c').replace('C', 'g')
    comp = gc_comp.upper()
    label_comp = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Reverse complement: {comp}')
    label_comp.pack(pady=7)

def protein(dna):
    code = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
            }
    last_codon_start = len(dna) - 2
    protein = ''
    for start in range(0, last_codon_start, 3):
        codon = dna[start:start + 3]
        aa = code.get(codon, '')
        protein += aa
    protein = protein
    label_protein = Label(frame_down_right, bg='alice blue', font='Verdana 13', text=f'Protein sequence: {protein}')
    label_protein.pack(pady=7)
if __name__ == "__main__":
    main()
