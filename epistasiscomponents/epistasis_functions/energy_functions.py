import numpy as np
import epistasiscomponents
from epistasiscomponents.constants import G_H, G_L2E, G_HDNA, MU_IPTG_RANGE, BETA, MU_OPERATOR, LOW_IPTG, HIGH_IPTG, MAX_ON, MIN_OFF
import pandas as pd

def relative_populations(dG_h=0,
         dG_l2e=0,
         dG_hdna=0,   
         G_h=G_H,
         G_l2e=G_L2E,
         G_hdna=G_HDNA,          
         mu_operator=3.5,
         mu_iptg=MU_IPTG_RANGE,                      
         beta=BETA):
    """
    Calculate the relative fraction of h, l2e, hdna
    
    dG_h: effect of mutations on H-state stability
    dG_l2e: effect of mutations on L-state + IPTG stability
    dG_hdna: effect of mutations on H-state + operator stability
    G_h: free energy of H-state at reference conditions
    G_l2e: free energy of L-state + IPTG state at reference conditions
    G_hdna: free energy of H-state + operator state at reference conditions
    mu_operator: chemical potential of operator
    mu_iptg: chemical potential of IPTG
    beta: 1/RT 
    
    Returns the fractional populations of h, l2e, and hdna
    """
    
    w_h =    np.exp(-(beta)*(G_h + dG_h))
    w_l2e =  np.exp(-(beta)*(G_l2e  + dG_l2e - 2*mu_iptg))
    w_hdna = np.exp(-(beta)*(G_hdna + dG_hdna - mu_operator))
    
    Z = w_h + w_l2e + w_hdna

    return w_h/Z, w_l2e/Z, w_hdna/Z

def dg_obs(dG_h=0,
           dG_l2e=0,
           dG_hdna=0,             
           G_h=G_H,
           G_l2e=G_L2E,
           G_hdna=G_HDNA,   
           mu_operator=3.5,
           mu_iptg=MU_IPTG_RANGE, 
           beta=BETA):      
    """
    Calculate the free energy of the operator bound state relative to 
    others. 
    
    dG_h: effect of mutations on H-state stability
    dG_l2e: effect of mutations on L-state + IPTG stability
    dG_hdna: effect of mutations on H-state + DNA stability
    G_h: free energy of H-state at reference conditions
    G_l2e: free energy of L-state + IPTG at reference conditions
    G_hdna: free energy of H-state + operator at reference conditions
    mu_operator: chemical potential of operator (single value)
    mu_iptg: chemical potential of iptg (range of values)
    beta: 1/RT
    
    Returns dG (H + L2E --> HDNA)
    """
    
    w_h, w_l2e, w_hdna = relative_populations(dG_h,
                              dG_l2e,
                              dG_hdna,
                              G_h,
                              G_l2e,
                              G_hdna,
                              mu_operator,
                              mu_iptg,
                              beta)
    
    return -(1/beta)*np.log(w_hdna/(w_h + w_l2e))



#I think to mimimize hte requ
def epistasis_vs_iptg(dG_h_m1,
                      dG_l2e_m1,
                      dG_hdna_m1,
                      dG_h_m2,
                      dG_l2e_m2,
                      dG_hdna_m2,
                      dG_h_m1m2=None,
                      dG_l2e_m1m2=None,
                      dG_hdna_m1m2=None,
                      G_h=G_H,
                      G_l2e=G_L2E,
                      G_hdna=G_HDNA,
                      mu_operator=MU_OPERATOR,
                      mu_iptg_range=MU_IPTG_RANGE,
                      beta=BETA):
    """
    Calculate epistasis for the system as a function of iptg concentration.
    
    dG_h_m1: effect of mutation 1 on h stability
    dG_l2e_m1: effect of mutation 1 on l2e stability
    dG_hdna_m1: effect of mutation 1 on hdna stability
    dG_h_m2: effect of mutation 2 on h stability
    dG_l2e_m2: effect of mutation 2 on l2e stability
    dG_hdna_m2: effect of mutation 2 on hdna stability
    dG_h_m1m2: effect of mutation 1 and 2 on h stability (if none, add effects)
    dG_l2e_m1m2: effect of mutation 1 and 2 on l2e stability (if none, add effects)
    dG_hdna_m1m2: effect of mutation 1 and 2 on hdna stability (if none, add effects)
    G_h: free energy of h state
    G_l2e: free energy of l2e state
    G_hdna: free energy of hdna state
    mu_operator: chemical potential of operator
    mu_iptg_range: chemical potential of iptg (range)
    beta: 1/RT
                        
    returns magnitue, sign1, and sign2 as a function of iptg
    """
    
    # If we are assuming mutations are additive within a state, 
    # throw out the (possibly) calculated dG12 value. 
    if dG_h_m1m2 is None:
        dG_h_m1m2 = dG_h_m1 + dG_h_m2
        dG_l2e_m1m2 = dG_l2e_m1 + dG_l2e_m2
        dG_hdna_m1m2 = dG_hdna_m1 + dG_hdna_m2
        
    # free energy of hdna for each genotype
    wt =  dg_obs(dG_h=0,
                 dG_l2e=0,
                 dG_hdna=0,                          
                 G_h=G_h,
                 G_l2e=G_l2e,
                 G_hdna=G_hdna,
                 mu_operator=mu_operator,
                 mu_iptg=mu_iptg_range,
                 beta=beta)

    m1 =  dg_obs(dG_h=dG_h_m1,
                 dG_l2e=dG_l2e_m1,
                 dG_hdna=dG_hdna_m1,  
                 G_h=G_h,
                 G_l2e=G_l2e,
                 G_hdna=G_hdna,
                 mu_operator=mu_operator,
                 mu_iptg=mu_iptg_range,
                 beta=beta)
    
    m2 =  dg_obs(dG_h=dG_h_m2,
                 dG_l2e=dG_l2e_m2,
                 dG_hdna=dG_hdna_m2, 
                 G_h=G_h,
                 G_l2e=G_l2e,
                 G_hdna=G_hdna,
                 mu_operator=mu_operator,
                 mu_iptg=mu_iptg_range,
                 beta=beta)
    
    m12 = dg_obs(dG_h=dG_h_m1m2,
                 dG_l2e=dG_l2e_m1m2,
                 dG_hdna=dG_hdna_m1m2,  
                 G_h=G_h,
                 G_l2e=G_l2e,
                 G_hdna=G_hdna,
                 mu_operator=mu_operator,
                 mu_iptg=mu_iptg_range,
                 beta=beta)  

    
    # magnitude of epistasis (signed)
    mag = (m12 - m2) - (m1 - wt)
    
    # Sign of mutation 1 (will be False if mutation effect has same
    # sign in both backgrounds; True if opposite signs)
    sign1 = (m12 - m2)/(m1 - wt)
    sign1 = sign1 < 0

    # Sign of mutation 2 (will be False if mutation effect has same
    # sign in both backgrounds; True if opposite signs)
    sign2 = (m12 - m1)/(m2 - wt)
    sign2 = sign2/np.abs(sign2) < 0
    
    return mag, sign1, sign2

"""
# Create dictionary keying mutation to their energetic effects on each state
unique_muts, unique_indexes = np.unique(df.m1,return_index=True)#np.unique(muts[:,1],return_index=True)
sorter = []
for i, m in enumerate(unique_muts):
    index = unique_indexes[i]
    sorter.append((int(m[1:-1]),
                   m,
                   df.loc[index,"h_mean_dG_m1"],
                   df.loc[index,"l2e_mean_dG_m1"],
                   df.loc[index,"hdna_mean_dG_m1"]))
sorter.sort()

mutant_energies = {}
for s in sorter:
    mutant_energies[s[1]] = [s[2],s[3],s[4]]
"""

def mutate_background(mutant_energies,
                      current_mutant=None,
                      low_iptg=LOW_IPTG,high_iptg=HIGH_IPTG,
                      low_max_on=MAX_ON,high_min_off=MIN_OFF,
                      G_h=G_H,G_l2e=G_L2E,G_hdna=G_HDNA,
                      mu_operator=MU_OPERATOR,
                      beta=BETA):
    """
    Introduce all possible mutations into a given mutant background and 
    return a list of mutations that did or did not pass 
    the screening conditions.
    
    mutant_energies: energetic effects of each mutation
    current_mutation: current mutation(s) (tuple of keys) in background
    G_h: free energy of H-state at reference conditions
    G_l2e: free energy of L-state + IPTG state at reference conditions
    G_hdna: free energy of H-state + operator state at reference conditions
    mu_operator: operator chemical potential
    beta: 1/RT
    low_iptg: IPTG chemical potential of "off" (low) screen
    high_iptg: IPTG chemical potential for "on" (high) screen
    low_max_on: fx_hdna must be at least at low IPTG this to pass
    high_min_off: fx_hdna must be below this at high IPTG to pass
    """
    
    no_pass = []
    screen_pass = []
    
    # Populate mutations
    if current_mutant is None:
        current_mutant = []
    else:
        current_mutant = list(current_mutant)
        
    current_sites = [m[1:-1] for m in current_mutant]
        
    # Get current mutant perturbation to dG for each state
    dG_h = 0
    dG_l2e = 0
    dG_hdna = 0
    for m in current_mutant:
        dG_h += mutant_energies[m][0]
        dG_l2e += mutant_energies[m][1]
        dG_hdna += mutant_energies[m][2]
        
    # Go through every possible mutant
    for m in mutant_energies:
        
        # Skip mutation if mutation already seen at this site
        if m[1:-1] in current_sites:
            continue
            

        # Update mutant 
        this_mutant = current_mutant[:]
        this_mutant.append(m)
        this_mutant = tuple(this_mutant)
            
        # Add mutant energies to the existing energies
        me =  mutant_energies[m]
        mutant_dG_h = dG_h + me[0] 
        mutant_dG_l2e = dG_l2e + me[1]
        mutant_dG_hdna = dG_hdna + me[2]
        
        # Calculate the fractional population of each conformation
        fx_h, fx_l2e, fx_hdna = relative_populations(dG_h=mutant_dG_h,
                                     dG_l2e=mutant_dG_l2e,
                                     dG_hdna=mutant_dG_hdna,
                                     G_h=G_h,
                                     G_l2e=G_l2e,
                                     G_hdna=G_hdna,
                                     mu_operator=mu_operator,
                                     mu_iptg=np.array([low_iptg,high_iptg]),
                                     beta=beta)
            
        # Record whether this would pass the screen conditions
        screen_state = (fx_hdna[0] > low_max_on,
                        fx_hdna[1] < high_min_off)
                
        # Update lists of mutants that do not pass screen, pass mid screen, or 
        # pass high screen
        total_fail = True
        if screen_state[0]:
            
            if screen_state[1]:
                screen_pass.append(this_mutant)
                total_fail = False
        
        if total_fail:
            no_pass.append(this_mutant)
            
    return no_pass, screen_pass


#This is what takes up all the time and will probably jsut end up geting gotten rid of
""" # Go through single and double point mutants. Start with "None" mutations
fails = [[None]]
broad = [[None]]
narrow = [[None]]

for mutant_round in range(2):
    
    # Screen all mutant pairs that passed the mid condition and put them 
    # through the mid condition again
    broad.append([])
    for b in tqdm(broad[-2]):
        no_pass, broad_pass = mutate_background(mutant_energies,
                                                current_mutant=b,
                                                low_iptg=LOW_IPTG,high_iptg=HIGH_IPTG)
        broad[-1].extend(broad_pass)
        
for mutant_round in range(3):

    # Screen all mutant pairs that passed the high condition and put them 
    # through the high condition again
    narrow.append([])
    for n in tqdm(narrow[-2]):
        no_pass, narrow_pass = mutate_background(mutant_energies,
                                                 current_mutant=n,
                                           low_iptg=NARROW_LOW_IPTG,high_iptg=NARROW_HIGH_IPTG)
        narrow[-1].extend(narrow_pass) """

def create_totals(h_total, l2e_total, hdna_total):
    nrg_sums = {
        'h' : h_total,
        'l2e' : l2e_total,
        'hdna' : hdna_total
    }
    return pd.DataFrame(nrg_sums)