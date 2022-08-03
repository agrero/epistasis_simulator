import random

import epistasiscomponents
#make this more concise later
from epistasiscomponents.epistasis_functions.gene_functions import create_gene_ndxs, mutate
from epistasiscomponents.constants import AMINO_ACIDS

class Sire:
    """
    The parent gene that will be undergoing mutations in order to create its scion
    """
    def __init__(self, sequence:str, name:str):

        self.sequence = sequence
        self.name = name
        self.deleterious_mutants = []

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
        seq_to_mutate = self.sequence

        mutants_list = []

        if len(self.sequence) < no_mutations:
            return []

        count = 0
        while count != no_mutations:
            if len(mutants_list) != 0:
                seq_to_mutate = mutants_list[-1]

            #insert energy check here for each mutant, returning to the start
            #if the conditions are not met
            #add deleterious mutations to the self.deleterious_mutants

            mutant_acid = random.choice(AMINO_ACIDS)
            mutant_ndx = random.choice(gene_ndxs)
            gene_ndxs.remove(mutant_ndx)

            mutation = (mutant_ndx, mutant_acid)
            mutant_sequence = mutate(seq_to_mutate, mutation)
            mutants_list.append(mutant_sequence)
            count += 1

        return mutants_list
    def __str__(self) -> str:
        return """Name: {}
Sequence: {}
        """.format(self.name, self.sequence)