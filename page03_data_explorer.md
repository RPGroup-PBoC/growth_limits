---
layout: page
title: Dataset Explorer
permalink: data_explorer
sidebar: true
interactive: data_explorer.html
---
---

## Data Explorer
This interactive figure presents the data assembled from the corresponding data
sets<sup>1-4</sup> after association of each protein with the appropriate
macromolecular complex, functional assignment, and standardization (i.e. scaling
to cell volume). This figure has two tabs -- 'Explore By Complex Function' and
'Explore By Individual Protein' -- which are clickable at the top of the
interactive and are described below:

### Explore By Complex Function

In this panel, you can explore each macromolecular complex detected across
the data sets. These complexes are broken down by their Clusters of
Orthologous Groups (COG) functional assignments. Selecting a COG category
and a macromolecular complex triggers visualization of several features of
the complex and its association with the data. 

**Plot of complex abundance vs growth rate**. This plot shows the number of
complexes present in each data set as a function of the growth rate. The number
of complexes is computed by default as the mean number functional complexes that
could be assembled given the observed abundances of each individual protein and
its stoichiometry in the macromolecular complex. The method of aggregation can
be selected by the aggregation button group to the right of the plot (minimum,
maximum, median, and mean are provided aggregation functions). Note that for
complexes consisting of a single protein, all aggregation functions produce the
same result

**EcoCyc complex record.** The complex identification number is given on the
right-hand side below the aggregation method buttons. This is presented as a
clickable link to the  EcoCyc web page associated with that complex

**Stoichiometric breakdown.** For each complex, the expected stoichiometric
ratios of each complex is given for each subunit in the table on the right-hand
side. Hovering your mouse over a point in the plot on the left-hand side will
populate this table with the observed stoichiometry for that data set and growth
condition. Note that the expected stoichiometry is computed relative to the
minimum subunit copy number, identified as the first subunit with stoichiometry
equal to 1 from the top of the table.


### Explore by Individual Protein

Like the complex visualization, all individual proteins present in the data sets
are assigned to a particular COG classification. Selecting a COG subclass
populates another drop-down menu where individual proteins can be selected,
populating the plot on the left-hand side. Below this dropdown is information
about the annotated function of the protein of interest.
<br/>
<center>

{% include_relative interactives/{{page.interactive}} %}

</center>
This interactive figure was made using the [Bokeh Python plotting library]().
The code and data used to generate this figure can be found here:

### Code
* [**`data_explorer.py`**](code/data_explorer.py) \| Python file for defining layout and plotting canvas
* [**`explorer_complex.js`**](code/explorer_complex.js) \| Javascript code for
  defining interaction with complex annotations.
* [**`explorer_subunits.js`**](code/explorer_subunits.js) \| Javascript code for
  defining hovering feature and population of stoichiometry table.
* [**`explorer_prot.js`**](code/explorer_prot.js) \| Javascript code for
  defining interaction with individual protein annotations.
* [**class_selection.js`**](code/class_selection.js) \| Javascript code for
  filtering data by COG class associations.
* [**`explorer_data_munging.py`**](code/explorer_data_munging.py) \| Python file
  used to generate data files for visualization. This script requires access to
  `compiled_annotated_complexes.csv` and `compiled_absolute_measurements.csv`,
  accessible on the [data]({{site.url}}/data) page of this website. 

### Data 
* [**`cplx_desc.csv`**](data/cplx_desc.csv) \| Data file containing list of
  complexes, their COG classifications, and their EcoCyc complex identification
  number.
* [**`cplx_numeric.csv`**](data/cplx_numeric.csv) \| Data file containing all
  measurements and aggregations of complex abundances.
* [**`prot_desc.csv`**](data/prot_desc.csv) \| Data file containing list of
  individual proteins, their COG classification, and their EcoCyc gene product
  annotation.
* [**`prot_numeric.csv`**](data/prot_numeric.csv) \| Data file containing all
  measurements of protein abundances.

## References 

1. Li, G.-W., Burkhardt, D., Gross, C., and Weissman, J. S. (2014). Quantifying absolute protein synthesis rates reveals principles underlying allocation of cellular resources. Cell, 157(3):624–635.
2. Peebo, K., Valgepea, K., Maser, A., Nahku, R., Adamberg, K., and Vilu, R. (2015). Proteome reallocation in *Escherichia coli* with increasing specific growth rate. Molecular BioSystems, 11(4):1184–1193.
3. Schmidt, A., Kochanowski, K., Vedelaar, S., Ahrné, E., Volkmer, B., Callipo, L., Knoops, K., Bauer, M., Aebersold, R., and Heinemann, M. (2016). The quantitative and condition-dependent *Escherichia coli* proteome. Nature Biotechnology, 34(1):104–110.
4. Valgepea, K., Adamberg, K., Seiman, A., and Vilu, R. (2013). *Escherichia coli* achieves faster growth by increasing catalytic and translation rates of proteins. Molecular BioSystems, 9(9):2344.
