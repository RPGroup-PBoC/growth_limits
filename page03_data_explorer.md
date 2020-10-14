---
layout: page
title: Dataset Explorer
permalink: data_explorer
sidebar: true
interactive: complex_explorer.html
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
This interactive figure was made using the [Bokeh Python plotting library](https://docs.bokeh.org/en/latest/index.html).
The code and data used to generate this figure can be found on the [GitHub Repository](https://github.com/RPGroup-PBoC/growth_limits/tree/publication/code/figures/complex_explorer).


## References 

1. Li, G.-W., Burkhardt, D., Gross, C., and Weissman, J. S. (2014). Quantifying absolute protein synthesis rates reveals principles underlying allocation of cellular resources. Cell, 157(3):624–635.
2. Peebo, K., Valgepea, K., Maser, A., Nahku, R., Adamberg, K., and Vilu, R. (2015). Proteome reallocation in *Escherichia coli* with increasing specific growth rate. Molecular BioSystems, 11(4):1184–1193.
3. Schmidt, A., Kochanowski, K., Vedelaar, S., Ahrné, E., Volkmer, B., Callipo, L., Knoops, K., Bauer, M., Aebersold, R., and Heinemann, M. (2016). The quantitative and condition-dependent *Escherichia coli* proteome. Nature Biotechnology, 34(1):104–110.
4. Valgepea, K., Adamberg, K., Seiman, A., and Vilu, R. (2013). *Escherichia coli* achieves faster growth by increasing catalytic and translation rates of proteins. Molecular BioSystems, 9(9):2344.
