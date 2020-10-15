# `complex_explorer`
This directory contains the code used to generate the [data explorer](https://rpgroup.caltech.edu/data_explorer) interactive figure. The contents are as follows:


*  `complex_explorer.py` generates the actual figure. 
*  `COG_selector.js` provides interactivity for filtering the complexes by their secondary COG category. 
*  `complex_selector.js` and `protein_selector.js` provide interactivity for selecting complexes or individual proteins, respectively. 
*  `complex_annotations.csv` is a data file which acts as a look up table of complex to the secondary COG category.  
*  `complexes_compressed.csv` and `proteins_compressed.csv` are the data sources for the complex and protein explorers, respectively. They are processed in a way that makes plotting easy. 
*  `explorer_munging.py` loads in the raw data and generates the two `csv` files described above. 