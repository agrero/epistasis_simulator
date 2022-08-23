import random
import pandas as pd
import numpy as np

#make this more concise later
from epistasiscomponents.epistasis_functions.gene_functions import create_gene_ndxs, mut_seq
from epistasiscomponents.constants import AMINO_ACIDS, DDG, HIGH_IPTG, LOW_IPTG, MAX_ON, MIN_OFF
from epistasiscomponents.epistasis_functions.energy_functions import dg_obs, mutate_background, relative_populations

class Sire:

    #this is also innaccurate
    """
    The parent gene that will be undergoing mutations in order to create its scion
    """
    def __init__(self, sequence:str, name:str):

        self.sequence = sequence
        self.name = name
        self.mutation_list = []
        self.progeny = ''

        self.nrgs = []

        self.nrg_totals = []

        self.mutant_pop_distro = [] 
        
    #in the future may want to add progeny size as a parameter
    def create_mutant(self, no_mutations):
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
                if mutant_acid == self.sequence[mutant_ndx]:
                    new_list = list(AMINO_ACIDS)
                    new_list.remove(mutant_acid)
                    mutant_acid = random.choice(new_list)
                gene_ndxs.remove(mutant_ndx)

                mutation = f"{self.sequence[mutant_ndx]}{mutant_ndx+1}{mutant_acid}" 
                mutant.append(mutation)
        
        self.mutation_list.append(mutant)

    def generate_progeny(self):
        mut_input = pd.DataFrame(self.mutation_list)
        mut_input['mut seq'] = mut_input.apply(lambda x: mut_seq(self.sequence, x), axis=1)
        return mut_input

    def get_mutant_energies(self):
        """
        Returns a list of mutation energies as derived from the ddg.csv. The mutations are 
        taken from the self.mutation_list.
        """

        if len(self.mutation_list) == 0:
            raise Exception("You cannot get a mutant's energy without any mutants.\n Call .create_mutant first.")

        energy_totals = []

        #should really figure out that enumerate thing
        #this is 100% the bottleneck
        count = 1
        for mutant in self.mutation_list:
            energy = DDG.loc[mutant] 
            self.nrgs.append(energy)
            
            energy_sum = energy.sum()
            totals = {
                'hdna' : energy_sum[0],
                'L2E' : energy_sum[1],
                'h' : energy_sum[2]
            }
            total_format = pd.DataFrame(totals, index = [f'{count}'])
            energy_totals.append(total_format)

            count += 1

        self.nrg_totals.append(pd.concat(energy_totals))

    #figure out how to vectorize later (Q and D answer right here)
    #honestly this is somewhat uneccesarry
    def get_mutant_distribution(self):
        """
        Returns a population distribution of each of the 3 forms of each mutant.
        """

        totals = self.nrg_totals[0]
        totals_mixed = totals.apply(lambda row: relative_populations(
            dG_h=row['h'], dG_l2e=row['L2E'], dG_hdna=row['hdna'],
            mu_iptg=np.array([LOW_IPTG,HIGH_IPTG])
        ), axis=1)
        distro = pd.DataFrame(totals_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_mixed.index)
        self.mutant_pop_distro = distro

    def mut_background(self, low_max_on=MAX_ON, high_min_off=MIN_OFF):
        """
        From generated mutant progeny test mutants for viability and screening conditions,
        more to be added later :d.
        """

        hdna_pop_dist_split = pd.DataFrame(self.mutant_pop_distro.hdna.to_list(), 
                                           columns=['pre', 'post'], 
                                           index=self.mutant_pop_distro.index)
        #the important thing to take from here is the indices as they correlate specifically to the mutants
        return hdna_pop_dist_split.query(f'pre > {low_max_on} and post < {high_min_off}')

    def __str__(self) -> str:
        return """
Name: {}
Sequence: {}
        """.format(self.name, self.sequence)