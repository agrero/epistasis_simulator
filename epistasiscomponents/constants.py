import numpy as np
import pandas as pd

def format_ddg(ddg, column_headers):
    """
    Reformat the indexing for ddg files to include an 'amino acid position' index.
    """
    individual_ndx = list(range(0,len(ddg)))
    aa_ndx = list(range(0, len(ddg), 20)) * 20
    aa_ndx.sort()
    arrays = [aa_ndx, individual_ndx]
    tuples = list(zip(*arrays))
    
    multi_dex = pd.MultiIndex.from_tuples(tuples, names=['amino acid', 'mut ndx'])

    return pd.DataFrame(ddg.values, index=multi_dex, columns=column_headers)

BETA =  (1/(298*0.001987)) # energy units in kcal/mol
MU_OPERATOR =  0 # operator chemical potential is 10 kcal/mol
G_H =          0 # h form is favored 
G_L2E  =       10 # 12e form is destabilized by 10 kcal/mol (@ IPTG Kd)
G_HDNA =       0 # hdna form destabilized by 20 kcal/mol (@ IPTG Kd and peptide Kd)
MU_IPTG_RANGE = np.arange(-2,12.1,.1) # IPTG chemical potential range for calculations in kcal/mol
conv = 1 # no conversion to other energy units (assumes Rosetta energy units are ~ kcal/mol. Kellogg et al, 2011).

#IPTG chemical potentials where cutoffs will be calculated for broad screen

LOW_IPTG = 0
HIGH_IPTG = 5


#IPTG chemical potentials where cutoffs will be calculated for narrow screen

NARROW_LOW_IPTG = .4
NARROW_HIGH_IPTG = .6

#MAX_ON: fx_hdna must be at least this at low screen point to pass
#MIN_OFF: fx_hdna must be below this at high screen point to pass

MAX_ON = 0.49
MIN_OFF = 0.01

#ddg for all possible mutations in the lac repressor
ddg_read = pd.read_csv("epistasiscomponents//ddg.csv")
DDG = format_ddg(ddg_read, column_headers=['mut','hdna', 'h', 'l2e'])

#The lac repressor amino sequence is here becuase it's ugly to look at
LAC_SEQ = "LLIGVATSSLALHAPSQIVAAIKSRADQLGASVVVSMVERSGVEACKTAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIALLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDDTEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLA"
