from itertools import chain
import numpy as np
import random

import pandas as pd
from Bio.Seq import MutableSeq

from epistasiscomponents.constants import AMINO_ACIDS, LAC_SEQ

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

#I think this can be replaced by just a simple list() call
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

def rep_sequence(sequence, no_replicates):
    new_list = []
    for acid in sequence:
        new_list.append([acid] * no_replicates)
    return list(chain(*new_list))

def generate_mutant_dataframe(sequence):
    """
    Generates a mutation notation for every possible mutation within a given sequence.
    """

    #generate all sequence indices for each amino acid. converts to string to be summable.
    seq_indices_int = list(range(1, len(sequence) + 1)) * 20
    seq_indices_int.sort()
    seq_indices_str = list(map(str, seq_indices_int))

    #creates the starting and finishing amino acids to be in the dataframe
    start_acids = rep_sequence(list(sequence), 20)
    end_acids = list(AMINO_ACIDS) * len(sequence)
    
    #creating the mutational dataframe
    mutations = pd.DataFrame({
        'start_acids' : start_acids,
        'seq_indices' : seq_indices_str,
        'end_acids' : end_acids
    })
    #sums parts along the row as to be searchable in the ddg file
    return mutations.sum(axis=1)

def get_samples(no_mutations, no_mutants, mutant_df):
    """
    This again, will have a many many things. 

    Hey, you're funny! Just thought I should let you know ;p
    """

    mut_ndx_ls = list(range(len(mutant_df)))

    mutations = []
    #change name later
    dub_dex = []

    for i in range(no_mutations * no_mutants):
        #repopulates the list after going through the correct amount of mutations
        if i % no_mutations == 0:
            mut_ndx_ls.clear()
            mut_ndx_ls += list(range(len(mutant_df)))

        pick = random.choice(mut_ndx_ls)
        upper_bound = (pick//20 + 1) * 20
        lower_bound = (pick//20) * 20

        mutations.append(pick)
        dub_dex.append(i//4)

        del mut_ndx_ls[lower_bound:upper_bound]

    arrays = [dub_dex, mutations]

    tuples = list(zip(*arrays))
    
    index = pd.MultiIndex.from_tuples(tuples, names=['mutant', 'mutation'])
   
    s = pd.Series(np.array(mutant_df.loc[mutations]), index = index)

    return s


def clever_column_rename(dataframe):
    """Im not 100% sold on the name"""
    no_col = len(dataframe.columns)
    
    names = []
    count = 0
    for i in range(no_col-3):
        mut_no = f'mut{count}'
        names.append(mut_no)
        count += 1
    names += ['hdna', 'h', 'l2e']

    return names

def screen_mutants(dataframe):
    genotype = dataframe['genotype']
    