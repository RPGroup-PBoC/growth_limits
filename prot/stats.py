import numpy as np
import pandas as pd

# Compute the total fractional occupancy. 
def compute_fraction(df, groupby, mass_key='fg_per_cell', count_key='tot_per_cell'):
    """
    Computes the fraction of the proteome occupied by each constitutent member
    both in terms of mass fraction and total copy number.

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe with proteomic information
    groupby: str
        Key by which to group the data to compute the fraction.
    mass_key: str
        Key by which to compute the mass fraction.
    count_key: str
        Key by which to compute the count fraction.

    Returns
    -------
    frac_df: pandas DataFrame
        DataFrame with columns corresponding to the group, fraction of proteome 
        by mass, and fraction of proteome by count. 
    """
    if type(groupby) != str:
       raise TypeError(f'Groupkey must be a string. A {type(groupby)} was passed') 
    if groupby.lower() not in df.keys():
        raise ValueError(f'Group not found in dataframe keys')

    # Compute the total mass
    tot_mass = df[mass_key].sum()
    tot_number = df[count_key].sum()

    # Instantiate a storage list of dataframes
    frac_df = pd.DataFrame([])
    for g, d in df.groupby(groupby):
        # Compute the fractional mass
        frac_mass = d[mass_key].sum() / tot_mass
        frac_count = d[count_key].sum() / tot_number
        frac_df = frac_df.append({'frac_mass':frac_mass,
                        'frac_count':frac_count,
                        'group':g}, ignore_index=True)
    return frac_df
                    