import random
import pandas as pd

#make this more concise later
from epistasiscomponents.epistasis_functions.gene_functions import create_gene_ndxs
from epistasiscomponents.constants import AMINO_ACIDS, DDG, HIGH_IPTG, LOW_IPTG, MAX_ON, MIN_OFF
from epistasiscomponents.epistasis_functions.energy_functions import dg_obs, relative_populations

class Sire:

    #this is also innaccurate
    """
    The parent gene that will be undergoing mutations in order to create its scion
    """
    def __init__(self, sequence:str, name:str):

        self.sequence = sequence
        self.name = name
        self.progeny = []

        self.nrgs = []
        #need to edit the energy method in order to organize this in a way which total(i) is the index and 
        #the columns are hdna, h, and l2e
        #100% can make the filter a df query
        self.nrg_totals = []

        self.mutant_pop_distribution = [] 

        self.no_pass_mutants = []
        self.broad_pass_mutants = []
        self.narrow_pass_mutants = []
        
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
                gene_ndxs.remove(mutant_ndx)

                mutation = f"{self.sequence[mutant_ndx]}{mutant_ndx+1}{mutant_acid}" 
                mutant.append(mutation)
        
        self.progeny.append(mutant)

    def get_mutant_energies(self):
        """
        Returns a list of mutation energies as derived from the ddg.csv. The mutations are 
        taken from the self.progeny.
        """

        if len(self.progeny) == 0:
            raise Exception("You cannot get a mutant's energy without any mutants.\n Call .create_mutant first.")

        energy_totals = []

        #should really figure out that enumerate thing
        #this is 100% the bottleneck
        count = 1
        for mutant in self.progeny:
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
            #totals = energy.iloc[:,0:].sum()
            count += 1

        self.nrg_totals.append(pd.concat(energy_totals))

    #figure out how to vectorize later (Q and D answer right here)
    #honestly this is somewhat uneccesarry
    def get_mutant_distribution(self):
        """
        Returns a population distribution of each of the 3 forms of each mutant.
        """

        distribution = []

        #replace with apply later
        for energy in self.nrgs:
           dist = relative_populations(dG_h=energy[1],
                                       dG_hdna=energy[0],
                                       dG_l2e=energy[2])
           self.mutant_pop_distribution.append(dist)

    def mut_background(self,
                       low_iptg=LOW_IPTG, high_iptg=HIGH_IPTG,
                       low_max_on=MAX_ON, high_min_off=MIN_OFF):
        """
        From generated mutant progeny test mutants for viability and screening conditions,
        more to be added later :d
        """

        

    def __str__(self) -> str:
        return """
Name: {}
Sequence: {}
        """.format(self.name, self.sequence)