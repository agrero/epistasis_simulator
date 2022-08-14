import random
import pandas as pd

#make this more concise later
from epistasiscomponents.epistasis_functions.gene_functions import create_gene_ndxs
from epistasiscomponents.constants import AMINO_ACIDS, DDG
from epistasiscomponents.epistasis_functions.energy_functions import dg_obs

class Sire:

    #this is also innaccurate
    """
    The parent gene that will be undergoing mutations in order to create its scion
    """
    def __init__(self, sequence:str, name:str):

        self.sequence = sequence
        self.name = name
        self.potential_progeny = []
        
    #in the future may want to add progeny size as a parameter
    def create_mutant(self , no_mutations):
        #rewrite this to how this actually works now
        """
        Generates a list of indices based on the length of the given sequence.
        The function then randomly picks an amino acid to be mutated, removing the index from the list.
        The function returns a dictionary denoting its lineage, with the key being the sequence
        and the values are derived from the specific mutations that occured to create the 
        mutant.
        """

        if len(self.sequence) < no_mutations:
            raise Exception("You cannot have more mutations than amino acids")

        mutable_range = range(len(self.sequence))
        gene_ndxs = create_gene_ndxs(mutable_range)

        mutant = []

        for mutation in range(no_mutations): #theres that enumerate thing that hsould go here

                mutant_acid = random.choice(AMINO_ACIDS)
                mutant_ndx = random.choice(gene_ndxs)
                gene_ndxs.remove(mutant_ndx)

                mutation = f"{self.sequence[mutant_ndx]}{mutant_ndx+1}{mutant_acid}" 
                mutant.append(mutation)
        
        self.potential_progeny.append(mutant)

    def mutant_energies(self):
        """
        Returns a list of mutation energies as derived from the ddg.csv. The mutations are 
        taken from the self.potential_progeny.
        """

        if len(self.potential_progeny) == 0:
            raise Exception("You cannot get a mutant energy without any mutants.\n Call .create_mutant first.")

        energies = []

        #may want to keep summary statistics seperate could make faster to parse through
        for mutant in self.potential_progeny:
            energy = pd.concat([DDG.loc[mutant], DDG.loc[mutant].sum(axis=0)], axis=1)
            energies.append(energy)


        #for energy in energies:
        #    pd.concat([energy, energy.sum(axis=0)])


        return energies

    def __str__(self) -> str:
        return """
Name: {}
Sequence: {}
        """.format(self.name, self.sequence)