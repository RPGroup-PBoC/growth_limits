---
layout: page
title: Model Explorer
permalink: model_explorer
sidebar: true
interactive: model_explorer.html
---
---

## Model Explorer

In the manuscript, we present a minimal model describing the relationship
between cell volume, growth rate, and ribosome abundance. We consider that the
rate of translation, termed the elongation rate $$r_t$$ is dependent on the
effective concentration of free amino acids $$[AA]_{eff}$$ and their affinity to
the elongating ribosomes, defined by a generalized dissociation constant
$$K_D$$. This magnitude of $$[AA]_{eff}$$ relative to $K_D$ tunes the elongation
rate between $$\approx 0$$ (complete cessation of translation) and some maximum
value $$r_t^{(max)}$$. This is defined mathematically as 

$$r_t = \frac{r_t^{(max)}}{1 + \frac{K_D}{[AA]_{eff}}}.$$ (1)


The effective amino acid concentration $$[AA]_{eff}$$ is defined as the
difference between the rate at which amino acids are produced $$r_{AA}$$ and the
rate at which they are consumed,

$$[AA]_{eff} = {1 \over V N_A} (r_{AA} - r_t \times R \times f_a),$$ (2)

where $$R$$ is the number of ribosomes, $$f_a$$ is the fraction of the ribosome
pool that is actively translating, and $$V N_A$$ is the product of the volume and
Avogadro's number to maintain units of concentration.

The cell growth rate $$\lambda$$ is related to the elongation rate $$r_t$$ by the number of
peptide bonds that need to be formed for a cell to divide, $$N_{pep}$$, via 

$$\lambda = \frac{r_t \times R \times f_a}{N_{pep}}.$$ (3)

In the interactive figure below, we provide a means to get a sense for the
quantitative predictions of this model. Use the sliders at the top of the figure
to adjust the corresponding parameter value and see how that changes the
translation rate and the cellular growth rate. 

<center>

{% include_relative interactives/{{page.interactive}} %}

</center>

This figure was made using the [Bokeh Python plotting framework](). The code
used to generate this figure can be downloaded below:

1. [**`model_explorer.py`**](code/model_explorer.py) \| This is a Python file
   which defines the layout and interactive widgets of the figure. 
2. [**`model_explorer.js`**](code/model_explorer.js) \| This is the Javascript
   code which defines the mathematics of the model and accommodates all
   interactions with the widgets.