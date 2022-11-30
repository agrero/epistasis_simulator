from epistasiscomponents.constants import DDG, LAC_SEQ
from epistasiscomponents.epistasis_objects.sire import Sire
import pandas as pd
import numpy as np
import random

from epistasiscomponents.epistasis_functions.gene_functions import clever_column_rename, create_gene_ndxs, generate_mutant_dataframe, get_samples, mut_seq
from epistasiscomponents.constants import AMINO_ACIDS, DDG, G_H, G_HDNA, G_L2E, HIGH_IPTG, LAC_SEQ, LOW_IPTG, MAX_ON, MIN_OFF, NARROW_HIGH_IPTG, NARROW_LOW_IPTG
from epistasiscomponents.epistasis_functions.energy_functions import dg_obs, mutate_background, relative_populations


#inefficient but works
def get_position(mutation_index):
    """
    Takes in a mutational index (or list of indices) and returns the amino acid position starting from 0.

    mutation_index: either an integer or list of integers that will subsequently be returned as integers representing
        amino acid posiion.
    
    """
    if type(mutation_index) == set or type(mutation_index) == list:
        mut_positions = []
        for mut in mutation_index:
            if mut % 20 == 0:
                mut_positions.append(((mut - 1) // 20))
            else:
                mut_positions.append(mut // 20)
        return mut_positions

    if mutation_index % 20 == 0:
        return (mutation_index - 1) // 20
    else:
        return mutation_index // 20


def position_check(ddg, genotype, current_mutation):
    """
    Parses a given ddg file to see if the newly added mutation shares a position with any previous mutations.
    
    ddg: ddg file you are pulling your mutational energies from.
    genotype: genotype as in a set of indicies from the ddg file that makeup the current mutations to your sequence.
    current_mutation: the most recent mutation that will be checked against the genotype to see if it overlaps positions. 
    """

    mutation_position = get_position(current_mutation)
    genotype_positions = [get_position(x) for x in genotype]
    


    #checks each genotype's posiition by dividing it by the number of amino acids possible
    # removes the mutation if the new one shares the same position
    # if the current mutation isn't overlapping with any other mutations the mutation is added to the genotype
#    for i in genotype:
#        print(i)
#        if i % 20 == 0 and mutation_position == ((i - 1) % 20):
#            genotype.discard(i)
#            genotype.add(current_mutation)
#        elif mutation_position == (i % 20):
#            genotype.discard(i)
#            genotype.add(current_mutation)
#        else:
#            genotype.add(current_mutation)



def sum_from_genotype(broad_mutant_dataframe, narrow_mutant_dataframe, all_ndxs, ddg):
    """
    Used in 'evolve' in order to get energy sums from differing dataframes.

    broad_mutant_dataframe: A dataframe of your mutants that survived the broad evolutionary screen.
    narrow_mutant_dataframe: A dataframe of your mutants that survived the narrow evolutionary screen.
    all_ndxs: The total number of possible functional mutations that the sequence can possibly undergo.
    ddg: ddg file you are pulling your mutational energies from.
    """

    #generate a genotype array to be added onto with subsequent mutations later
    broad_genotypes = broad_mutant_dataframe['genotype']
    narrow_genotypes = narrow_mutant_dataframe['genotype']
    
    #This needs to change just
    #I think all you need to do is take all of your current genotypes, add on an additional column 
    # in order to show the last added mutation as to remove the not 'last added'
    if len(broad_genotypes) > len(narrow_genotypes): 
        full_new_muts = list(np.random.randint(low=0, high=len(all_ndxs), size=len(broad_genotypes)))
        length_difference = len(narrow_genotypes) - len(broad_genotypes)
        short_new_muts = full_new_muts[0:length_difference]

        for i in enumerate(broad_genotypes):
            i[1].update({full_new_muts[i[0]]})
        for i in enumerate(narrow_genotypes):
            i[1].update({short_new_muts[i[0]]})


    elif len(broad_genotypes) < len(narrow_genotypes):
        full_new_muts = list(np.random.randint(low=0, high=len(all_ndxs), size=len(narrow_genotypes)))
        length_difference = len(broad_genotypes) - len(narrow_genotypes)
        short_new_muts = full_new_muts[0:length_difference]

        for i in enumerate(broad_genotypes):
            i[1].update({short_new_muts[i[0]]})
        for i in enumerate(narrow_genotypes):
            i[1].update({full_new_muts[i[0]]})
    
    else:
        full_new_muts = list(np.random.randint(low=0, high=len(all_ndxs), size=len(narrow_genotypes)))
        for i in enumerate(broad_genotypes):
            i[1].update({full_new_muts[i[0]]})
        for i in enumerate(narrow_genotypes):
            i[1].update({full_new_muts[i[0]]})

   
        

def evolve(ddg, no_mut_bg, no_mut_per_round=1, no_evolutionary_rounds=1, broad_rep=[], narrow_rep=[], low_max_on=MAX_ON, high_min_off=MIN_OFF):
    """
    Evolutionary simulation of energetic effects of mutatations to a protein sequence.
    Mutants will be screened and placed into two different categories (broad or narrow)
    based on the distribution of their thermodynamic ensembles. Evolutionary trajectories are then
    represented as a venn diagram showing the similarities between the specific sequences located
    in each screen.

    ddg: Your ddg file specifying all energetic effects of mutations to your sequence
    no_mut_bg: The number of initial mutations for the first mutants, IE double or triple.
    no_mut_per_round: Number of mutants per round.
    """
    

    #Creating a valid background
    all_ndxs = list(range(len(ddg)))
    #chagne to list comprehension later
    ndxs = []
    for i in range(no_mut_per_round):
        ndx = []
        for i in range(no_mut_bg):
            ndx.append(all_ndxs[random.randint(0, len(all_ndxs)-1)])
        ndxs.append(ndx)

    init_trajectories = []
    for i in ndxs:
        if len(i) != no_mut_bg:
            i.append(all_ndxs[random.randint(0, len(all_ndxs)-1)])
        mut = ddg.iloc[i, :]
        genotype = set((i))
        mut_nrgs = mut.iloc[:,1:]
        mut_sum = mut_nrgs.sum()

        mut_sum['genotype'] = genotype

        init_trajectories.append(mut_sum)

    y = pd.concat(init_trajectories, axis=1)
    z = pd.DataFrame(y)
    g = z.swapaxes('columns', 'index')
    #background created

    #get mutant distributions
    #broad trajectories
    totals_broad_mixed = g.apply(lambda row: relative_populations(
        dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
        mu_iptg=np.array([LOW_IPTG, HIGH_IPTG])
        ), axis=1)
    temp_broad_distro = pd.DataFrame(totals_broad_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_broad_mixed.index)
    l2e_broad_distro = pd.DataFrame(temp_broad_distro.l2e.to_list(), index = temp_broad_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
    h_broad_distro = pd.DataFrame(temp_broad_distro.h.to_list(), index = temp_broad_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
    hdna_broad_distro = pd.DataFrame(temp_broad_distro.hdna.to_list(), index= temp_broad_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])
    init_broad_distro = pd.concat([l2e_broad_distro, h_broad_distro, hdna_broad_distro], axis=1)

    broad_rep.append(init_broad_distro)
    #narrow trajectories
    totals_narrow_mixed = g.apply(lambda row: relative_populations(
        dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
        mu_iptg=np.array([NARROW_LOW_IPTG, NARROW_HIGH_IPTG])
        ), axis=1)
    temp_narrow_distro = pd.DataFrame(totals_narrow_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_narrow_mixed.index)

    l2e_narrow_distro = pd.DataFrame(temp_narrow_distro.l2e.to_list(), index = temp_narrow_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
    h_narrow_distro = pd.DataFrame(temp_narrow_distro.h.to_list(), index = temp_narrow_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
    hdna_narrow_distro = pd.DataFrame(temp_narrow_distro.hdna.to_list(), index= temp_narrow_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])
            
    init_narrow_distro = pd.concat([l2e_narrow_distro, h_narrow_distro, hdna_narrow_distro], axis=1)
    
    narrow_rep.append(init_narrow_distro)

    n_and_b_distros = [init_broad_distro, init_narrow_distro]
    #screen mutants
    # the 0 index is for broad distro passes and the 1 index is for narrow distro passes    
    init_screen_pass_ndxs = []

    for mutant_frame in n_and_b_distros:
        filter1 = mutant_frame["HDNA LOW IPTG"] < low_max_on
        filter2 = mutant_frame["HDNA HIGH IPTG"] > high_min_off
        mutant_frame.where(filter1, inplace=True)
        mutant_frame.where(filter2, inplace=True)
        init_screen_pass_ndxs.append(mutant_frame.dropna().index)

        #energies = self.trajectories[mutant_frame].loc[mutant_frame.dropna().index]
        #self.broad_pass_energies.append(energies)


    initial_broad_mutants = g.loc[init_screen_pass_ndxs[0]]
    initial_narrow_mutants = g.loc[init_screen_pass_ndxs[1]]
    
    print(initial_broad_mutants)
    print(initial_narrow_mutants)
 
    #iterate and make more mutant progeny!
    broad_muts = []
    narrow_muts = []

    for iteration in range(no_evolutionary_rounds):
        current_broad_mutants = initial_broad_mutants
        current_narrow_mutants = initial_narrow_mutants

        sum_from_genotype(current_broad_mutants, current_narrow_mutants, all_ndxs)



  

broad_r = []
narrow_r = []

#evolve(DDG, no_mut_bg=2, no_mut_per_round=100, no_evolutionary_rounds=1, broad_rep=broad_r, narrow_rep=narrow_r)
index = [7, 17, 34]
values = [[set({3520, 717})], [set({905, 467})], [set({5364, 20})]]
position_check_frame = pd.DataFrame(data=values, index=index)

position_check_frame.apply(lambda row: position_check(
    ddg=DDG, genotype=row, current_mutation=20
))


