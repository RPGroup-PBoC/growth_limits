# Fundamental Limits on the Rate of Bacterial Growth


## Branches
This repository contains three branches -- `master`, `gh-pages`, and
`publication`. The `master` branch contains the entire project history for
start to finish and includes code and data that were exploratory or may not
have made it into the final publication as well as the LaTeX document files.
The `gh-pages` branch hosts all files associated with the [paper
website](https://rpgroup.caltech.edu/growth_limits) and contains the source
files for the final compiled dataset as well as all scripts used for the
figure generation. The final branch, `publication` (where you are now)
strikes a middle ground between `master` and `gh-pages`. This branch contains
all data and code used to generate figures, clean and process data, and
compute estimations as necessary. This branch should be sufficient to
reproduce everything presented in this work.

## Installation
To run the code used in this work, you will need to install the `prot` module. `prot` 
is a custom Python package written specifically for this project and contains an 
array of functions used for everything from stylizing the plots to computing cell 
volumes as a function of growth rate. The requirements can be installed by executing the folliwng command using `pip` in the command line:

> pip install -r requirements.txt

The software module itself can be installed locally by executing the command in the root directory,

> pip install -e ./

When installed, a new folder `prot.egg-info` will be made in the root directory and is necessary 
to run any of the code. 

## `PathwayTools`
A key component of our work is scraping the [EcoCyc](https://ecocyc.org/) database to associated individual  proteins present in the datasets with the associated complexes. In order to do so,
we used the [PathwayTools](http://bioinformatics.ai.sri.com/ptools/) software API. 
In order to run this part of the analysis (`code/processing/ecocyc_subunit_munging/munge.py`), 
you must have an active PathwayTools license. 

## Repository Architecture
This repository is broken up into several directories and subdirectories. Please
see each directory for information regarding each file.

### `code`
All Python code used in this work is located here and is written to be
executed from its local directory. This directory contains two subdirectories
in which the code is organized.

1. **`processing`** | All code that was used in the cleaning of the raw datasets,
   interaction with databases, and collation of the complete datasets. This
   directory is broken up further into directories associated with each individual 
   dataset or process.
2. **`figures`** | All code used to generate the figures, both main text and supplemental.
    The scripts are written to be executed from this directory.

### `data` 
This self-explanatory folder contains all of the data, both raw and processed,
used in this research. 

## License
![](https://licensebuttons.net/l/by/3.0/88x31.png)


All creative work (text, plots, images, etc.) are licensed under the
[Creative Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/)
license. All software is distributed under the standard MIT license as
follows:

```
Copyright (c) 2020 The Authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```