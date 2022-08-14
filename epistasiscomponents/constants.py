import numpy as np
import pandas as pd

BETA =  (1/(298*0.001987)) # energy units in kcal/mol
MU_OPERATOR =  0 # operator chemical potential is 10 kcal/mol
G_H =          0 # h form is favored 
G_L2E  =       10 # 12e form is destabilized by 10 kcal/mol (@ IPTG Kd)
G_HDNA =       0 # hdna form destabilized by 20 kcal/mol (@ IPTG Kd and peptide Kd)
MU_IPTG_RANGE = np.arange(-2,12.1,.1) # IPTG chemical potential range for calculations in kcal/mol
conv = 1 # no conversion to other energy units (assumes Rosetta energy units are ~ kcal/mol. Kellogg et al, 2011).

#IPTG chemical potentials where cutoffs will be calculated for broad screen

LOW_IPTG = 0
HIGH_IPTG = 10


#IPTG chemical potentials where cutoffs will be calculated for narrow screen

NARROW_LOW_IPTG = 3
NARROW_HIGH_IPTG = 5

#MAX_ON: fx_hdna must be at least this at low screen point to pass
#MIN_OFF: fx_hdna must be below this at high screen point to pass

MAX_ON = 0.49
MIN_OFF = 0.01

#Tuple of all amino acid single letter codes
AMINO_ACIDS = ('G', 'A', 'L', 'M', 'F',
               'W', 'K', 'Q', 'E', 'S',
               'P', 'V', 'I', 'C', 'Y',
               'H', 'R', 'N', 'D', 'T')

#ddg for all possible mutations in the lac repressor
DDG = pd.read_csv("epistasiscomponents//ddg.csv", index_col=0)

#The lac repressor amino sequence is here becuase it's ugly to look at
LAC_SEQ = "LLIGVATSSLALHAPSQIVAAIKSRADQLGASVVVSMVERSGVEACKTAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIALLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDDTEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLA"