import numpy as np
from . import size

def load_constants(lambda_bounds = [0, 2]):
    """
    Loads and returns a dictionary with a variety of numeric constants

    Parameters
    ----------
    lambda_bounds : list
        List of length=2 entries that define the upper and lower bounds of the
        grwoth rate range to consider. Should be in units of hr^-1

    Returns
    -------
    constants : dict
        A dictonary of commonly used constants. Each entry is a dictionary 
        with three keys -- value, units, and  source. Source will include either the BNID 
        associated with that entry or the literature source
        
        density : Cell density in units of pg / cubic micron
        volume : Cell volume in units of cubic micron. Given as a function of
                growth rate 
        growth_rate : User provided range of growth rates
        dry_mass_frac : Fraction of dry mass 
        theta_C : Fraction of dry mass that is carbon
        theta_P : Fraction of dry mass that is phosphorus
        theta_S : Fraction of dry mass that is sulfur
        theta_prot: Fraction of dry mass that is protein        
        theta_DNA: Fraction of dry mass that is DNA
        theta_RNA : Fraction of dry mass that is RNA
        theta_lipid: Fraction of dry mass that is lipid
        N_ori : Number of origins as calculated in Si et al. 2017
        rate_DNAP: Polymerization rate of DNAP
        rate_RNAP: Polymerization rate of RNAP
        L_genome: Length of single E. coli chromosome
        t_div: Time required to copy a single genome
    """
    growth_rate = np.linspace(lambda_bounds[0], lambda_bounds[1], 200)

    # Define basal sources
    constants = {
        'density': {'value':1.1, 'units': 'pg/um^3', 'source': 'BNID:103875'},
        'volume': {'value':size.lambda2size(growth_rate), 
                    'units':'um^3', 'source':'defined quantity'},
        'growth_rate': {'value': growth_rate, 
                    'units':'hr^-1', 'source':'defined quantity'},
        'dry_mass_frac': {'value': 0.3, 'units':'fractional', 'source': 'BNID:100214'},
        'theta_C':  {'value': 0.5, 'units':'fractional', 'source': 'BNID:100649'},
        'theta_P': {'value': 0.03, 'units': 'fractional', 'source': 'BNID:100653'},
        'theta_S': {'value': 0.01, 'units': 'fractional', 'source': 'BNID:100655'},
        'theta_prot': {'value': 0.5, 'units': 'fractional', 'source': 'BNID:101436'},
        'theta_DNA': {'value': 0.03, 'units': 'fractional', 'source': 'BNID:101436'},
        'theta_lipid': {'value':0.1, 'units': 'fractional', 'source': 'BNID:101436'},
        'theta_RNA': {'value': 0.2, 'units': 'fractional', 'source': 'BNID:101436'},
        'rate_DNAP': {'value': 600, 'units':'nt/s', 'source': 'BNID:104928'},
        'rate_RNAP': {'value': 40, 'units': 'nt/s', 'source': 'BNID:101904'},
        'L_genome': {'value': 4.6E6, 'units': 'bp', 'source': 'BNID:100269'},
        }

    # Define calculated quanities
    constants['t_double'] = {'value':3600 * np.log(2) / constants['growth_rate']['value'],
                              'units': 's', 'source':'calculated quantity'}
    constants['cell_mass'] = {'value': constants['density']['value'] * constants['volume']['value'],
                              'units': 'pg', 'source':'calculated quantity'}
    constants['t_div'] = {'value': constants['L_genome']['value'] / constants['rate_DNAP']['value'] / 2,
                          'units': 's', 'source': 'calculated quantity'}
    constants['N_ori'] = {'value': 2**(3600 * np.log(2)/growth_rate / constants['t_div']['value']),
                          'units': 'number', 'source':'calculated quantity (see Si et al. 2017)'}
    return constants

     
