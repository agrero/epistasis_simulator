import epistasiscomponents.constants as con

import pandas as pd
import numpy as np

from epistasiscomponents.epistasis_functions.gene_functions import position_check


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

init_stuff = pd.Series([[5183, 2336], [3471, 3339], [1616, 743]])
new_stuff = pd.Series([5182, 428, 1616])
new_frame = pd.concat([init_stuff, new_stuff], axis=1)

new_f = pd.DataFrame(new_frame.loc[:,0].tolist())
new_f['end'] = new_stuff

new_data = new_f.values.tolist()

new_new_f = pd.DataFrame([new_data])
swap = new_new_f.swapaxes('columns', 'index')

position_checked = swap.apply(lambda row: position_check(row), axis=1)


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
    
