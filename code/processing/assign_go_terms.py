"""
This script takes a dataframe of absolute protein measurements and breaks
each gene down by its annotated go term. This results in a large dataframe that 
is highly redundant, but will be useful for generating plots of molecular 
complex abundance as a function of growth rate. 

TODO: Figure out how to properly deal with quaternary structure and subunit
abundance
"""
#%%
import numpy as np
import pandas as pd 
import tqdm as tqdm.