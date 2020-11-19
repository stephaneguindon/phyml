
[![Build Status](https://travis-ci.org/stephaneguindon/phyml.svg?branch=master)](https://travis-ci.org/stephaneguindon/phyml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/phyml/README.html) 


### Overview

PhyML is a software package that uses modern statistical approaches to analyse alignments of nucleotide or amino acid sequences in a phylogenetic
framework. The main tool in this package builds phylogenies under the maximum likelihood criterion. It implements a large number of substitution
models coupled to efficient options to search the space of phylogenetic tree topologies. PhyTime is another tool in the PhyML package that focuses
on divergence date estimation in a Bayesian setting. The main strengths of PhyTime lies in its ability to accommodate for uncertrainty in the placement of fossil calibration
and the use of realistic models of rate variation along the tree. Finally, PhyREX fits the spatial-Lambda-Fleming-Viot
model to geo-referenced genetic data. This model is similar to the structured coalescent but assumes that individuals are distributed along a spatial continuum rather
than discrete demes. PhyREX can be used to estimate population densities and rates of dispersal. Its output can be processed by treeannotator (from the BEAST package) as well as
SPREAD.

### Citations

- New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0 S Guindon, JF Dufayard, V Lefort, M Anisimova, W Hordijk, O Gascuel *Systematic Biology* 59 (3), 307-321, 2010.
- Accounting for calibration uncertainty: Bayesian molecular dating as a “doubly intractable” problem S Guindon *Systematic Biology* 67 (4), 651–661, 2018.
- Demographic inference under the coalescent in a spatial continuum S Guindon, H Guo, D Welch *Theoretical Population Biology* 111, 43-50, 2016.


### Installation

To install any program that is part of the PhyML package, type the following command:

```bash
sh ./autogen.sh;
```

If you are using a Mac computer or running a Unix-like operating system, you will need 
to install the packages autoconf automake and pkg-config. On a Mac, the following command 
should set you up (provided Homebrew is installed on your Mac...): brew install pkg-config
autoconf automake;

Next, to install any program that is part of the PhyML package, type the following commands:

```bash
./configure --enable-XXXX;
make;
```
where XXXX is phyml or phyrex or phytime.

To compile a Windows executable, install MinGW and run:

```bash
./configure --enable-win --enable-XXXX;
make;
```

To install the MPI version of PhyML, type the following commands:

```bash
autoreconf -i;
./configure --enable-phyml-mpi;
make;
```

