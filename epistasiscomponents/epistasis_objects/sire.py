from array import array
import random
import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go


#make this more concise later
from epistasiscomponents.epistasis_functions.gene_functions import clever_column_rename, create_gene_ndxs, generate_mutant_dataframe, get_samples, mut_seq
from epistasiscomponents.constants import AMINO_ACIDS, DDG, G_H, G_HDNA, G_L2E, HIGH_IPTG, LAC_SEQ, LOW_IPTG, MAX_ON, MIN_OFF, NARROW_HIGH_IPTG, NARROW_LOW_IPTG
from epistasiscomponents.epistasis_functions.energy_functions import dg_obs, mutate_background, relative_populations


class Sire:
    """
    The Sire class takes in a sequence and known mutative ddG values to create a sampling of mutant scion.
    The scion are then analyzed based on the presence of each component of their thermodynamic ensemble, 
    which can then be screened based on effector concentrations. 
    """
    def __init__(self, sequence:str, name:str, ddG_file:str):
        
        self.sequence = sequence
        self.name = name
        self.ddg = ddG_file
        self.mutation_list = ''
        self.progeny = ''

        self.nrgs = []
        self.nrg_totals = []
        self.screened_nrg_totals = ''

        self.mutant_pop_distro = [] 
        self.valid_background = ''
        self.trajectories = []

        self.broad = []
        self.narrow = []

    def create_valid_background(self, no_init_bgs):
        valid_m1_bg = self.ddg.query(f'hdna < {G_HDNA} and h < {G_H} and l2e < {G_L2E}')
        
        all_ndxs = list(range(len(valid_m1_bg)))
        ndx = []
        no_iter = no_init_bgs
        if no_init_bgs > len(all_ndxs):
            no_iter = len(all_ndxs)
            print(f'The number of backgrounds chosen is more than there are available, defaulting to: {len(all_ndxs)}')

        for i in range(no_iter):
            ndx.append(all_ndxs.pop(random.randint(0, len(all_ndxs)-1)))
        self.valid_background = valid_m1_bg.iloc[ndx,:]
   
    def new_iter_trajectory(self, no_init_bgs, no_iter_bgs, no_generations=1):
        valid_bg = pd.DataFrame(np.repeat(self.valid_background.values, len(self.ddg), axis=0), columns=self.valid_background.columns)
        y = valid_bg.groupby(by='mut', group_keys=True).apply(lambda x: x)
        count = no_generations - 1
        no_init_bgs = no_init_bgs
        #switch to len(valid_bg) once two iterations work

        ddg_replicate_2 = pd.concat([self.ddg]*no_init_bgs, ignore_index=True).rename({'mut':'mut1','hdna':'hdna1','h':'h1','l2e':'l2e1'}, axis=1)
        new_trajectory = pd.concat([y, ddg_replicate_2], axis=1, ignore_index=True).rename({0:'bg',1:'hdna',2:'h',3:'l2e',
                                                                                            4:'mut1',5:'hdna1',6:'h1',7:'l2e1'}, axis=1)
        new_trajectory = new_trajectory.assign(hdna_sum=lambda x: x.hdna + x.hdna1,
                                               h_sum=lambda x: x.h + x.h1,
                                               l2e_sum=lambda x: x.l2e + x.l2e1).query(f'hdna_sum < {G_HDNA} and h_sum < {G_H} and l2e_sum < {G_L2E}')
        while count != 0:

            #dub_dex = [new_trajectory.loc[:,'bg'], list(new_trajectory.index)]
            #tuples = list(zip(*dub_dex))
            #index = pd.MultiIndex.from_tuples(tuples, names=['bg', 'ndx'])

            to_drop = ['hdna', 'h', 'l2e', 'hdna1', 'h1', 'l2e1']
            new_values = new_trajectory.drop(to_drop, axis=1)

            col_names = clever_column_rename(new_values)
            s = pd.DataFrame(np.array(new_values), columns=col_names)
            s = s.drop_duplicates()

            no_init_bgs = no_iter_bgs
            
            ndx_choices = []
            for i in range(no_iter_bgs):
                ndx_choices.append(random.choice(list(range(len(s))))) 

            pre_new_trajectory = pd.DataFrame(s.iloc[ndx_choices])
            self.trajectories.append(pre_new_trajectory)
            #change the input here 

            no_init_bgs = no_iter_bgs

            mut_no = f'mut{pre_new_trajectory.shape[1]-3}'
            ddg_replicate = pd.concat([self.ddg]*len(pre_new_trajectory.index), ignore_index=True).rename({'mut':mut_no,'hdna':'hdna1','h':'h1',
                                                                                                           'l2e': 'l2e1',}, axis=1)
            pre_new_trajectory = pd.DataFrame(np.repeat(pre_new_trajectory.values, len(self.ddg), axis=0), columns=pre_new_trajectory.columns)
            new_trajectory = pd.concat([pre_new_trajectory, ddg_replicate], axis=1, ignore_index=False).reset_index(drop=True)
            new_trajectory = new_trajectory.assign(hdna_sum=lambda x: x.hdna + x.hdna1,
                                               h_sum=lambda x: x.h + x.h1,
                                               l2e_sum=lambda x: x.l2e + x.l2e1).query(f'hdna_sum < {G_HDNA} and h_sum < {G_H} and l2e_sum < {G_L2E}')

            count -= 1
        to_drop = ['hdna', 'h', 'l2e', 'hdna1', 'h1', 'l2e1']
        new_values = new_trajectory.drop(to_drop, axis=1)

        col_names = clever_column_rename(new_values)
        s = pd.DataFrame(np.array(new_values), columns=col_names)
        s = s.drop_duplicates()

        ndx_choices = []

        for i in range(no_iter_bgs):
            ndx_choices.append(random.choice(list(range(len(s)))))

        pre_new_trajectory = pd.DataFrame(s.iloc[ndx_choices]) 

        self.trajectories.append(pre_new_trajectory)

        
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

        for mutation in range(no_mutations): 
                mutant_acid = random.choice(AMINO_ACIDS)
                mutant_ndx = random.choice(gene_ndxs)
                if mutant_acid == self.sequence[mutant_ndx]:
                    new_list = list(AMINO_ACIDS)
                    new_list.remove(mutant_acid)
                    mutant_acid = random.choice(new_list)
                gene_ndxs.remove(mutant_ndx)

                mutation = f"{self.sequence[mutant_ndx]}{mutant_ndx+1}{mutant_acid}" 
                mutant.append(mutation)
        
        self.mutation_list = mutant

    def create_mutant_v2(self, no_mutations:int, progeny_size:int):
        """
        Stuff goes here! Much stuff. Informative stuff! More than what is currently here!
        """
        mut_df = generate_mutant_dataframe(self.sequence)
        samples = get_samples(no_mutations, progeny_size, mut_df)

        self.mutation_list = samples.drop_duplicates()
        

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
            energy = self.ddg.loc[mutant] 
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

    def get_mutant_energies_v2(self):

        if len(self.mutation_list) == 0:
            raise Exception("You cannot get a mutant's energy without any mutants.\n Call .create_mutant first.")

        df3 = self.mutation_list.to_frame('mut').join(self.ddg, on='mut')
        self.nrgs = df3
        self.nrg_totals = df3.groupby(by='mutant').sum()

    def get_mutant_distribution(self):
        """
        Returns a population distribution of each of the 3 forms of each mutant.
        """

        totals = self.nrg_totals
        totals_mixed = totals.apply(lambda row: relative_populations(
            dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
            mu_iptg=np.array([LOW_IPTG,HIGH_IPTG])
        ), axis=1)
        distro = pd.DataFrame(totals_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_mixed.index)
        self.mutant_pop_distro = distro

    def get_mutant_distribution_v2(self):
        totals = self.trajectories
        for i in totals:
            #broad trajectories
            totals_broad_mixed = i.apply(lambda row: relative_populations(
                dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
                mu_iptg=np.array([LOW_IPTG, HIGH_IPTG])
            ), axis=1)
            temp_broad_distro = pd.DataFrame(totals_broad_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_broad_mixed.index)
            l2e_broad_distro = pd.DataFrame(temp_broad_distro.l2e.to_list(), index = temp_broad_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
            h_broad_distro = pd.DataFrame(temp_broad_distro.h.to_list(), index = temp_broad_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
            hdna_broad_distro = pd.DataFrame(temp_broad_distro.hdna.to_list(), index= temp_broad_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])
            broad_distro = pd.concat([l2e_broad_distro, h_broad_distro, hdna_broad_distro], axis=1)

            #wuery can't search for low
            #broad_distro_2 = broad_distro.query(f'HDNA LOW IPTG > {MIN_OFF} and HDNA HIGH IPTG < {MAX_ON}')
            self.broad.append(broad_distro)


            totals_narrow_mixed = i.apply(lambda row: relative_populations(
                dG_h=row['h'], dG_l2e=row['l2e'], dG_hdna=row['hdna'],
                mu_iptg=np.array([NARROW_LOW_IPTG, NARROW_HIGH_IPTG])
            ), axis=1)
            temp_narrow_distro = pd.DataFrame(totals_narrow_mixed.to_list(), columns=['h', 'l2e', 'hdna'], index=totals_narrow_mixed.index)

            l2e_narrow_distro = pd.DataFrame(temp_narrow_distro.l2e.to_list(), index = temp_narrow_distro.index, columns=['L2E LOW IPTG', 'L2E HIGH IPTG'])
            h_narrow_distro = pd.DataFrame(temp_narrow_distro.h.to_list(), index = temp_narrow_distro.index, columns=['H LOW IPTG', 'H HIGH IPTG'])
            hdna_narrow_distro = pd.DataFrame(temp_narrow_distro.hdna.to_list(), index= temp_narrow_distro.index, columns=['HDNA LOW IPTG', 'HDNA HIGH IPTG'])
            
            narrow_distro = pd.concat([l2e_narrow_distro, h_narrow_distro, hdna_narrow_distro], axis=1)
            self.narrow.append(narrow_distro)

            

    def mut_background(self, low_max_on=MAX_ON, high_min_off=MIN_OFF):
        """
        From generated mutant progeny test mutants for viability and screening conditions,
        more to be added later :d.
        """

        hdna_pop_dist_split = pd.DataFrame(self.mutant_pop_distro.hdna.to_list(), 
                                           columns=['pre', 'post'], 
                                           index=self.mutant_pop_distro.index)
        #the important thing to take from here is the indices as they correlate specifically to the mutants
        screen_pass = hdna_pop_dist_split.query(f'pre > {low_max_on} and post < {high_min_off}')
        self.screened_nrg_totals = self.nrg_totals.loc[screen_pass.index]

    def plot_3d(self, rotate=False, lin_max=100):
        fig = plt.figure(figsize= plt.figaspect(0.5))
        ax = fig.add_subplot(1, 2, 1, projection = '3d')

        xs = self.screened_nrg_totals.loc[:,'h']
        ys = self.screened_nrg_totals.loc[:,'hdna']
        zs = self.screened_nrg_totals.loc[:,'l2e']
        ax.scatter(xs, ys, zs)

        xf = np.linspace(0, lin_max)
        yf = np.linspace(0, lin_max)
        zf = np.linspace(0, lin_max)
        ax.plot(xf,yf,zf, linewidth=2, color='r')

        ax.set_xlabel('h')
        ax.set_ylabel('hdna')
        ax.set_zlabel('L2E')

        ax = fig.add_subplot(1,2,2, projection='3d')
        ax.scatter(xs,ys,zs)

        ax.set_xlabel('h')
        ax.set_ylabel('hdna')
        ax.set_zlabel('L2E')

        xf = np.linspace(0, lin_max)
        yf = np.linspace(0, lin_max)
        zf = np.linspace(0, lin_max)
        ax.plot(xf,yf,zf, linewidth=2, color='r')

        if rotate:
            for angle in range(0, 360):
                ax.view_init(10,angle)
                plt.draw()
                plt.pause(.001)
        else:
            plt.show()

        plt.close()

    def violin_plot(self):
        fig = go.Figure()

        narrow = self.narrow[0]
        broad = self.broad[0]

        fig.add_trace(go.Violin(
            x=list(narrow.columns[2]),
            y=narrow['HDNA LOW IPTG'],
            legendgroup='Yes', scalegroup='Yes', name='Yes',
            side='negative',
            line_color='blue'))

        fig.add_trace(go.Violin(x=list(broad.columns[2]),
                                y=broad['HDNA LOW IPTG'],
                                legendgroup='No', scalegroup='No', name='No',
                                side='positive',
                                line_color='orange')
        )
        fig.update_traces(meanline_visible=True)
        fig.update_layout(violingap=0, violinmode='overlay')
        fig.show()

    def __str__(self) -> str:
        return """
Name: {}
Sequence: {}
        """.format(self.name, self.sequence)