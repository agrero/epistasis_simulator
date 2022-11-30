from epistasiscomponents.constants import DDG, LAC_SEQ
from epistasiscomponents.epistasis_objects.sire import Sire
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

sire = Sire(LAC_SEQ, 'lac', DDG)

no_mutants = 100

sire.create_valid_background(no_mutants)

sire.new_iter_trajectory(10, 10, 3)
#sire.get_mutant_distribution_v2()
#sire.screen_mutants()




#matching = []
#only_narrow = []
#only_broad = []


#put on to do list, change how you're looking at which sequences stay per round, the indexes change with each subsequent mutation added,
#even though the original base sequences are the same


#absolutely hideous
#fix later

"""for i in enumerate(sire.trajectories):
    print(sire.broad_pass_energies[i[0]])
    series1 = list(sire.broad_pass_energies[i[0]].index)
    series2 = list(sire.narrow_pass_energies[i[0]].index)

    for i in series1:
        if i in series2:
            matching.append(i)
        else:
            only_narrow.append(i)

    for i in series2:
        if not i in series1:
            only_broad.append(i)

    print(matching)
    print(only_narrow)
    print(only_broad)

    venn2(subsets = (len(only_broad), len(only_narrow), len(matching)), set_labels = ['Only Broad', 'Only Narrow'])
    plt.show()

    matching.clear()
    only_narrow.clear()
    only_broad.clear()
    series1.clear()
    series2.clear()"""
    
