import pandas as pd
from Bio.Seq import MutableSeq

def mutate(sequence:str, mutation:tuple) -> str:
    """
    Takes a tuple of two indices in the following format (O, M), and the sequence of which you want to mutate
    O: An integer that denotes the index of which amino acid you wish to mutate
    M: A string that denotes which amino acid you will be mutating O into
    """

    try:
        sequence[mutation[0]]
    except:
        return sequence
    
    mutation_pieces = [sequence[0:mutation[0]], mutation[1], sequence[mutation[0]+1:]]
    return ''.join(mutation_pieces)
 
def create_gene_ndxs(ran:range) -> list:
    gene_ndxs = []
    for i in ran:
        gene_ndxs.append(i)
    return gene_ndxs

def mut_seq(sequence, mut_list):
    """
    Takes an input sequence and an array of mutatations (a list may also work),
    and changes the sequence accordingly. Returns mutated sequence.
    """
    mut_seq = list(sequence)
    for mut in list(mut_list):
        mut_ndx = int(mut[1:-1])
        mut_seq[mut_ndx-1] = mut[-1]        
    return ''.join(mut_seq)

