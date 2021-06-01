.. EPCY documentation master file, created by
   sphinx-quickstart on Fri Nov 20 10:21:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EPCY: Evaluation of Predictive CapabilitY
=========================================

EPCY is a method used to rank genes (features) according to their
potential as predictive (bio)markers, using quantitative data (like
gene expression).

**Installation**

Using PyPI

.. code:: bash

   pip install epcy
   epcy -h


**Overview**

Similarly to Differential Expression analyses, EPCY take as input two
tabulated files:

* a table which describes the comparative design of the experience
* a matrix which contains quantitative data for each gene (feature) and sample

Using these two input files, EPCY will evaluate the predictive
capacity of each gene individually and return predictive scores, along
with their confidence intervals.

To guarantee the reliability of predictive scores, EPCY uses a leave-one-out
cross validation to train multiple Kernel Density Estimation (`KDE`_)
classifiers and evaluate their performances on unseen samples
(see `method <https://epcy.readthedocs.io/en/latest/method.html>`_ for more
details).

**Background**

EPCY is a product of the `Leucegene project <leucegene.ca>`_ and has
been developed and tested specifically to analyse RNA-seq data of acute
myeloid leukemia (AML) patients. However, the method implemented in
EPCY is generic and should work on different quantitative data.

**Citing**

While we are finalizing work on the official paper, more details can be found
`in a poster presented at ISMB ECCB 2019 <https://f1000research.com/posters/8-1349>`_:

Audemard E, Sauvé L and Lemieux S. EPCY: Evaluation of Predictive
CapabilitY for ranking biomarker gene candidates [version 1; not peer reviewed].
*F1000Research 2019*, 8(ISCB Comm J):1349(poster)

.. _KDE: https://en.wikipedia.org/wiki/Kernel_density_estimation

.. toctree::
   :maxdepth: 2
   :caption: Documentation index

   basic_usage
   method
   predictive_capability_columns
   rna_tools
   cutoff
   small_dataset
