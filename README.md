# epistasis_simulator

Example usage:

mutant_dictionary = evolve(ddg=con.DDG, no_mutations_background=2, no_mutations_per_round=20, no_evolutionary_rounds=4, 
        low_max_on=con.MAX_ON, high_min_off=con.MIN_OFF,
        denaturation_check=False)

This would output a dictionary with two different items with the keys 'broad' and narrow' of which correspond to the broad and narrow screened mutants.

Each dictionary item is a list (with this specific one being of length 4 due to the number of evolutionary rounds) with the 0th index corresponding
to the background, and each increasing index is the next round of mutations added.

An example of how to check for differences in sequences through each round can be found below: 

for i, item in enumerate(panda['broad']):
    comb = pd.concat([item['genotype'], panda['narrow'][i]['genotype']])
    
    no_duplicates = comb.duplicated().sum()
    no_original_broad = len(comb.drop_duplicates(keep=False))
    no_original_narrow = len(panda['narrow'][i]['genotype']) - no_duplicates

You can track the number of shared sequences (no_duplicates), broad mutants (no_original_broad), and narrow (no_original_narrow). 
Preliminary screening showed a pretty homogenous sequence space throughout each round of evolution, but that is likely due to most mutations are somewhat nonsense.

For the Future:
I would like to filter out most of the nonsense mutations for each round of evolution as it would likely give us more fruitful results.

Taking a more object oriented approach may help with some of the organization. 

Some features are also still specified for the 3 states of laci only, however I would like to make the program more generic where it could analyze any number of 
states that may be present in the ddg file.
