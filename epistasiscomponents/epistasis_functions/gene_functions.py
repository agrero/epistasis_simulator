from itertools import chain
import numpy as np
import pandas as pd

from epistasiscomponents.constants import HIGH_IPTG, LOW_IPTG, MAX_ON, MIN_OFF, NARROW_HIGH_IPTG, NARROW_LOW_IPTG
from epistasiscomponents.epistasis_functions.energy_functions import relative_populations
import epistasiscomponents.epistasis_functions.gene_functions as gf

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


def mutate(sequence:str, mutation:tuple) -> str:
    """
    Takes in a sequence and mutates the specific position
    Sequence: Sequence to mutate.
    Mutation: Amino acid change.
    """

    try:
        sequence[mutation[0]]
    except:
        return sequence
    
    mutation_pieces = [sequence[0:mutation[0]], mutation[1], sequence[mutation[0]+1:]]
    return ''.join(mutation_pieces)

def create_gene_ndxs(ran:range) -> list:
    """
    Creates a list of indexes for a given range.
    """
    gene_ndxs = []
    for i in ran:
        gene_ndxs.append(i)
    return gene_ndxs

def mut_seq(sequence, mut_list):
    """
    Takes an input sequence and an array of mutatations (a list may also work),
    and changes the sequence accordingly. Returns mutated sequence.
    Different from mutate as it takes in a list of mutations instead of single amino acid.
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

def degenerate_mutation_check(mutation_list):
    """
    Takes in a list of mutational indices from a ddg file. These are then tested to see if there is any repeat amino acid positions.
    If there is a degenerate position, the older mutation is removed.
    """
    #generate a list of aa positions from the list of mutations
    temp_mutation_list = list(mutation_list)[0]
    position_list = [get_position(x) for x in temp_mutation_list]

    #check the position list for degenerate mutations to be removed
    repeats = []
    for i,item in enumerate(position_list):
        test_dex = position_list[:i] + position_list[i+1:]
        if item in test_dex:
            repeats.append(i)
    if len(repeats) > 0:
        del temp_mutation_list[repeats[0]]

    #if there is degeneracy it'll return the truncated list, otherwise nothing happens
    return temp_mutation_list



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


def position_check(genotype, current_mutation, first_pass, ndxs):
    """
    Parses a given ddg file to see if the newly added mutation shares a position with any previous mutations.
    
    ddg: ddg file you are pulling your mutational energies from.
    genotype: genotype as in a set of indicies from the ddg file that makeup the current mutations to your sequence.
    current_mutation: the most recent mutation that will be checked against the genotype to see if it overlaps positions.
    first_pass: on the first pass of the evolutionary screen make true
    ndxs: indices to work with if the first pass has degenerate mutations
    """
    if first_pass:
        if len(genotype) == 0:
            return random.choice(ndxs)

        else: 
            mutation_position = get_position(current_mutation)
            genotype_positions = [get_position(x) for x in genotype]

            if mutation_position in genotype_positions:

                lower_bound = (mutation_position - 1) * 20
                upper_bound = mutation_position * 20
                used_ndxs = ndxs[lower_bound:upper_bound]
                used_and_all_ndxs = used_ndxs + ndxs

                return random.choice(list(set((used_and_all_ndxs))))
        
            else:
                return random.choice(ndxs)

    mutation_position = get_position(current_mutation)
    genotype_positions = [get_position(x) for x in genotype]

    new_genotype = []

    #do this in a more elegant way later
    for position in enumerate(genotype_positions):
        mutant = []
        if mutation_position == position[1]:
            mutant.append(current_mutation)
        else:
            new_genotype.append()

def get_energies_from_genotype(genotype, ddg):
    """
    Takes in a list of dgg indices and returns their sum.
    
    genotype: a list of indices from the ddg file
    ddg: ddg file of a specific protein
    """

    try:
        mut = ddg.iloc[genotype]
    except:
        broken_array = np.array(genotype)
        temp_genotype = broken_array[np.logical_not(np.isnan(broken_array))]
        mut = ddg.iloc[temp_genotype]
    mut_sum = mut.sum()

    return mut_sum[1:]    

def evolve(ddg, 
            no_mutations_background, no_mutations_per_round, no_evolutionary_rounds,
            low_max_on=MAX_ON, high_min_off=MIN_OFF,
            denaturation_check=False):
            """
            Evolutionary simulation of energetic effects of mutatations to a protein sequence.
            Mutants will be screened and placed into two different categories (broad or narrow)
            based on the distribution of their thermodynamic ensembles. Each sequence that passes will be returned as a dictionary that contains 
            a list of list for all of broad and narrow screen passes for each round.

            ddg: Your ddg file specifying all energetic effects of mutations to your sequence
            no_mutations_background: The number of initial mutations for the first mutants, IE double or triple.
            no_mutations_per_round: Number of mutants per round.
            no_evolutionary_rounds: The number of evolutionary rounds you want to run. Each round adds another mutation on to your base sequences.
            denaturation_check: Unstable, if True it will eliminate any proteins of which should be 'destabilized.'
            """


            #creating a valid background
            all_ndxs = list(range(len(ddg)))
            
            bg_ndxs = []
            for i in range(no_mutations_per_round):
                ndx=[]
                for i in range(no_mutations_background):
                    ndx.append(position_check(genotype=ndx, current_mutation=all_ndxs[random.randint(0, len(all_ndxs)-1)],
                    first_pass=True, ndxs=all_ndxs))
                bg_ndxs.append(ndx)
            
            initial_trajectories = []
            for i in bg_ndxs:
                if len(i) != no_mutations_background:
                    i.append(all_ndxs[random.randint(0, len(all_ndxs)-1)])

                mut = ddg.iloc[i, :]

                genotype = i
                mut_nrgs = mut.iloc[:,1:]
                mut_sum = mut_nrgs.sum()
                mut_sum['genotype'] = genotype

                initial_trajectories.append(mut_sum)
            
            joined_initial_trajectories = pd.concat(initial_trajectories, axis=1)
            joined_dataframe = pd.DataFrame(joined_initial_trajectories)
            mutational_background = joined_dataframe.swapaxes('columns', 'index')
            #background created
                
            #getting mutant distributions
            # broad distributions
            totals_broad_mixed = mutational_background.apply(lambda row: relative_populations(
                dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
                mu_iptg=np.array([LOW_IPTG, HIGH_IPTG])
                ), axis=1)
            temp_broad_distro = pd.DataFrame(totals_broad_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_broad_mixed.index)

            l2e_broad_distro = pd.DataFrame(temp_broad_distro.l2e.to_list(), index = temp_broad_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
            h_broad_distro = pd.DataFrame(temp_broad_distro.h.to_list(), index = temp_broad_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
            hdna_broad_distro = pd.DataFrame(temp_broad_distro.hdna.to_list(), index= temp_broad_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])

            init_broad_distro = pd.concat([l2e_broad_distro, h_broad_distro, hdna_broad_distro], axis=1)
            
            # narrow distributions
            totals_narrow_mixed = mutational_background.apply(lambda row: relative_populations(
                dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
                mu_iptg=np.array([NARROW_LOW_IPTG, NARROW_HIGH_IPTG])
                ), axis=1)
            temp_narrow_distro = pd.DataFrame(totals_narrow_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_narrow_mixed.index)

            l2e_narrow_distro = pd.DataFrame(temp_narrow_distro.l2e.to_list(), index = temp_narrow_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
            h_narrow_distro = pd.DataFrame(temp_narrow_distro.h.to_list(), index = temp_narrow_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
            hdna_narrow_distro = pd.DataFrame(temp_narrow_distro.hdna.to_list(), index= temp_narrow_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])
            
            init_narrow_distro = pd.concat([l2e_narrow_distro, h_narrow_distro, hdna_narrow_distro], axis=1)
            #initial distributions created

            #screen mutants
            narr_and_broad_distros = [init_broad_distro, init_narrow_distro]
            
            # the 0 index is for the broad passes and the 1 index is for the narrow passes
            init_screen_pass_ndxs = []

            for mutant_frame in narr_and_broad_distros:
                filter1 = mutant_frame["HDNA LOW IPTG"] < low_max_on
                filter2 = mutant_frame["HDNA HIGH IPTG"] > high_min_off
                mutant_frame.where(filter1, inplace=True)
                mutant_frame.where(filter2, inplace=True)
                init_screen_pass_ndxs.append(mutant_frame.dropna().index)

            initial_broad_mutants = mutational_background.loc[init_screen_pass_ndxs[0]]
            initial_narrow_mutants = mutational_background.loc[init_screen_pass_ndxs[1]]
            #initial screened mutant background created

            #begin iterating 
            # repository for mutants
            broad_muts = []
            narrow_muts = []

            broad_muts.append(initial_broad_mutants)
            narrow_muts.append(initial_narrow_mutants)

            print('finished creating background')

            for iteration in range(no_evolutionary_rounds):

                it_no = 0
                print(f'starting iteration {it_no}')

                current_broad_mutants = broad_muts[-1]
                current_narrow_mutants = narrow_muts[-1]
                
                denaturation_pass_ndxs = []
                #checkign for denaturation
                #not working, needs some fixing, but kinda helped us see that it's actually working
                
                if denaturation_check:
                    current_muts = [current_broad_mutants, current_narrow_mutants]
                    for mutant_frame in current_muts:
                        filter_1 = mutant_frame['hdna'] < 20
                        filter_2 = mutant_frame['h'] < 20
                        filter_3 = mutant_frame['l2e'] < 20
                        mutant_frame.where(filter_1, inplace=True)
                        mutant_frame.where(filter_2, inplace=True)
                        mutant_frame.where(filter_3, inplace=True)
                        denaturation_pass_ndxs.append(mutant_frame.dropna().index)

                    current_broad_mutants = broad_muts[-1].loc[denaturation_pass_ndxs[0]]
                    current_narrow_mutants = narrow_muts[-1].loc[denaturation_pass_ndxs[0]]


                ndxs = []
                for i in range(no_mutations_per_round):
                    ndxs.append(random.choice(all_ndxs))
                
                # expanding mutants for merging                
                broad_repeat = pd.DataFrame(np.repeat(current_broad_mutants.values, len(ndxs), axis=0), columns=current_broad_mutants.columns)
                narrow_repeat = pd.DataFrame(np.repeat(current_narrow_mutants.values, len(ndxs), axis=0), columns=current_narrow_mutants.columns)
           
                # expanding ndxs for merging
                broad_ndxs_repeat = ndxs * len(current_broad_mutants)

                #broad_ndxs_repeat.sort()
                narrow_ndxs_repeat = ndxs * len(current_narrow_mutants)
                #narrow_ndxs_repeat.sort()
                ndxs.clear()

                # merging new and current mutations
                broad_mixed = pd.concat([broad_repeat['genotype'], pd.Series(broad_ndxs_repeat)], axis=1)
                narrow_mixed = pd.concat([narrow_repeat['genotype'], pd.Series(narrow_ndxs_repeat)], axis=1)

                broad_merged = pd.DataFrame(broad_mixed.loc[:,'genotype'].tolist())
                broad_merged['end'] = broad_ndxs_repeat
                narrow_merged = pd.DataFrame(narrow_mixed.loc[:,'genotype'].tolist())
                narrow_merged['end'] = narrow_ndxs_repeat

                #  unpacking the merged dataframe to then input it as a list
                broad_merged_listed = broad_merged.values.tolist()
                narrow_merged_listed = narrow_merged.values.tolist()

                broad_remerged_wrong_ax = pd.DataFrame([broad_merged_listed])
                narrow_remerged_wrong_ax = pd.DataFrame([narrow_merged_listed])
                
                #it works up til here
                broad_remerged = broad_remerged_wrong_ax.swapaxes('columns', 'index')
                narrow_remerged = narrow_remerged_wrong_ax.swapaxes('columns', 'index')

                broad_position_checked = broad_remerged.apply(lambda row: gf.degenerate_mutation_check(row), axis=1)
                narrow_position_checked = narrow_remerged.apply(lambda row: gf.degenerate_mutation_check(row), axis=1)
                
                broad_energies = broad_position_checked.apply(lambda row: get_energies_from_genotype(row, ddg))
                broad_energy_genotype = pd.concat([broad_energies, broad_position_checked], axis=1)

                narrow_energies = narrow_position_checked.apply(lambda row: get_energies_from_genotype(row, ddg))
                narrow_energy_genotype = pd.concat([narrow_energies, narrow_position_checked], axis=1)
   
                broad_energy_genotype.rename(columns={0:'genotype'}, inplace=True)
                narrow_energy_genotype.rename(columns={0:'genotype'}, inplace=True)
                
                #getting mutant distributions
                # broad distributions
                totals_broad_mixed = mutational_background.apply(lambda row: relative_populations(
                    dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
                    mu_iptg=np.array([LOW_IPTG, HIGH_IPTG])
                    ), axis=1)
                temp_broad_distro = pd.DataFrame(totals_broad_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_broad_mixed.index)

                l2e_broad_distro = pd.DataFrame(temp_broad_distro.l2e.to_list(), index = temp_broad_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
                h_broad_distro = pd.DataFrame(temp_broad_distro.h.to_list(), index = temp_broad_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
                hdna_broad_distro = pd.DataFrame(temp_broad_distro.hdna.to_list(), index= temp_broad_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])

                current_broad_distro = pd.concat([l2e_broad_distro, h_broad_distro, hdna_broad_distro], axis=1)
                
                
                # narrow distributions
                totals_narrow_mixed = mutational_background.apply(lambda row: relative_populations(
                    dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
                    mu_iptg=np.array([NARROW_LOW_IPTG, NARROW_HIGH_IPTG])
                    ), axis=1)
                temp_narrow_distro = pd.DataFrame(totals_narrow_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_narrow_mixed.index)

                l2e_narrow_distro = pd.DataFrame(temp_narrow_distro.l2e.to_list(), index = temp_narrow_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
                h_narrow_distro = pd.DataFrame(temp_narrow_distro.h.to_list(), index = temp_narrow_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
                hdna_narrow_distro = pd.DataFrame(temp_narrow_distro.hdna.to_list(), index= temp_narrow_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])
                
                current_narrow_distro = pd.concat([l2e_narrow_distro, h_narrow_distro, hdna_narrow_distro], axis=1)
                #initial distributions created

                current_distros = [current_broad_distro, current_narrow_distro]
                curr_screen_pass_ndxs = []
                
                #screen mutants
               
                for mutant_frame in current_distros:
           
                    filter1 = mutant_frame["HDNA LOW IPTG"] < low_max_on
                    filter2 = mutant_frame["HDNA HIGH IPTG"] > high_min_off
                    mutant_frame.where(filter1, inplace=True)
                    mutant_frame.where(filter2, inplace=True)
                    
                    curr_screen_pass_ndxs.append(mutant_frame.dropna().index)
                
                current_broad_mutants = broad_energy_genotype.loc[curr_screen_pass_ndxs[0]]
                current_narrow_mutants = narrow_energy_genotype.loc[curr_screen_pass_ndxs[1]]
                
                broad_muts.append(current_broad_mutants)
                narrow_muts.append(current_narrow_mutants)

                print(f'finished iteration {it_no}')
                it_no += 1
                
            
            dict_return = {
                'broad' : broad_muts,
                'narrow' : narrow_muts
            }
            return dict_return
