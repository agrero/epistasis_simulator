import epistasiscomponents.constants as con

import pandas as pd


from clean import evolve_v2

import matplotlib.pyplot as plt



#matching = []
#only_narrow = []
#only_broad = []


#put on to do list, change how you're looking at which sequences stay per round, the indexes change with each subsequent mutation added,
#even though the original base sequences are the same


#absolutely hideous
#fix later

mutant_dictionary = evolve_v2(ddg=con.DDG, no_mutations_background=2, no_mutations_per_round=20, no_evolutionary_rounds=4, 
        low_max_on=con.MAX_ON, high_min_off=con.MIN_OFF,
        denaturation_check=False)

for i, item in enumerate(panda['broad']):
    comb = pd.concat([item['genotype'], panda['narrow'][i]['genotype']])
    
    no_duplicates = comb.duplicated().sum()
    no_original_broad = len(comb.drop_duplicates(keep=False))
    no_original_narrow = len(panda['narrow'][i]['genotype']) - no_duplicates
    
    print(no_duplicates)
    print(no_original_broad)
    print(no_original_narrow)

    


