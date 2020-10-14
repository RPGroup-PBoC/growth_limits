---
layout: page
title: Data
permalink: code
img: seg.png
sidebar: true
---

---

The data used in the main text of this work primarily come from four unique recent 
proteomic measurements. These publications are 

* Schmidt et al. 2016. *The quantitative condition-dependent Escherichia coli proteome*. Nature Biotechnology **34**(1). DOI: [10.1038/nbt.3418](https:/doi.org/10.1038/nbt.3418)
* Peebo et al. 2015. *Proteome reallocation in Escherichia coli with increasing specific growth rate.* Molecular BioSystems **11**. DOI: [10.1039/C4MB00721B](https://doi.org/10.1039/C4MB00721B)
* Li et al. 2014. *Quantifying absolute protein synthesis rates reveals principles underlying allocation of cellular resources.* Cell **157**(3). DOI:[10.1016/j.cell.2014.02.033](https://doi.org/10.1016/j.cell.2014.02.033)
* Valgepea et al. 2013. *Escherichia coli achieves faster growth by increasing catalytic and translation rates of proteins*. Molecular BioSystems **9**. DOI: [10.1039/C3MB70119K](https://doi.org/10.1039/C3MB70119K)

All of the code used in the processing of the raw data sets to make them directly comparable
and comparable between growth rates can be found on the corresponding [GitHub Repository](https://github.com/rpgroup-pboc/growth_limits). The result of this processing are three different files with measurements of the abundances of individual proteins, complexes, or classes of complexes. These (tidy) data sets, with headers describing the annotation, can be downloaded below:

{% if site.data.datasets %}
{% for ds in site.data.datasets %}
* [{{ds.name}}]({%if ds.storage !=
  'remote'%}{{site.url}}/{{site.baseurl}}/datasets/{{ds.link}}{%
  else%}{{site.link}}{% endif %}) \| {% if ds.filetype %}(filetype:
  {{ds.filetype}}){%endif%}{% if ds.filesize %}({{ds.filesize}}){%endif%}{%
  if ds.storage == remote %} DOI: {{ds.DOI}}{%endif%}
{% endfor %}
{% endif %}
