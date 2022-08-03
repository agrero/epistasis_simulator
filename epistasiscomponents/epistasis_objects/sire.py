import random

import epistasiscomponents
#make this more concise later
from epistasiscomponents.epistasis_functions.gene_functions import create_gene_ndxs, mutate

class Sire:
    """
    The parent gene that will be undergoing mutations in order to create its scion
    """
    def __init__(self, sequence:str, name:str):

        self.sequence = sequence
        self.name = name

    #in the future may want to add progeny size as a parameter
    def create_mutants(self, no_mutations):
        """
        Generates a list of indices based on the length of the given sequence.
        The function then randomly picks an amino acid to be mutated, removing the index from the list.
        The function returns a dictionary denoting its lineage, with the key being the sequence
        and the values are derived from the specific mutations that occured to create the 
        mutant.
        """
        mutable_range = range(len(self.sequence))
        gene_ndxs = create_gene_ndxs(mutable_range)

        mutants_list = []

        count = 0
        while count != no_mutations:
            mutant_ndx = random.choice(gene_ndxs)
            gene_ndxs.remove(mutant_ndx)
            mutation = (mutant_ndx, 'G')
            mutant_sequence = mutate(self.sequence, mutation)
            mutants_list.append(mutant_sequence)
            count += 1

        return mutants_list
    def __str__(self) -> str:
        return """Name: {}
Sequence: {}
        """.format(self.name, self.sequence)