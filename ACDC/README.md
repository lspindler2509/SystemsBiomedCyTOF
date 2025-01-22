# **A**utomated **C**ell type **D**iscovery and **C**lassification through knowledge transfer #
### ACDC for CyTOF data processing ###

ACDC is a Python library designed for single-cell data annotation, utilizing a combination of clustering, Gaussian mixture models, and random walk classifiers.
ACDC leverages a two-step process:
* Landmark Generation: Prior knowledge of cell types is encoded using a marker-cell table and transformed into high-dimensional landmarks.
* Semi-Supervised Classification: A random walk algorithm propagates information from these landmarks to classify all single-cell events in the data.

This is a Python 3 implementation of methods used in ["Automated cell type discovery and classification through knowledge transfer"] (https://pmc.ncbi.nlm.nih.gov/articles/PMC5447237/), a automated method for analyzing CyTOF data. This repository modified the [original ACDC repository](https://bitbucket.org/dudleylab/acdc/src/master/) in a way that it supports modern Python versions. Also, this modified code allows for paralellization and also allows for multiple samples and metadata and has additional features.


### Setup ###
* Dependencies
    * This package will work for Python versions >= 3.7
    * pandas, numpy, scipy, scikit-learn, seaborn, matplotlib,fcsy,flowio
    * third-party package [Phenograph](https://github.com/jacoblevine/PhenoGraph)


* Installation
Ensure you are in the directory containing the setup.py file. Then run to install the ACDC package locally.: 

```bash
pip install -e .
```

Alternatively, you can do it in an environment:

```bash
python3 -m venv acdc_env
source acdc_env/bin/activate  # On Windows: acdc_env\Scripts\activate
pip install -e .
```

### Tutorials ###
The original results from ACDC publication coming from the BMMC dataset were also reproduced and ACDC was run on new datasets as well. Notebooks for all datasets and preprocessing is included in the notebooks/ directory.