
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
 
def create_gene_ndxs(range:range) -> list:
    gene_ndxs = []
    for i in range:
        gene_ndxs.append(i)
    return gene_ndxs