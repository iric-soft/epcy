.. EPCY documentation master file, created by
   sphinx-quickstart on Fri Nov 20 10:21:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EPCY: Evaluation of Predictive CapabilitY
=========================================

Evaluation of Predictive CapabilitY (EPCY) is a method used to rank genes
(features) according to their potential as predictive (bio)markers,
using quantitative data (like gene expression).

Similarly to Differential Expression analyses, EPCY take in input two
tabulated files:

* a table which describe the comparative design of the experience
* a matrix which contains quantitative data for each gene (feature) and sample

Using these two input files, EPCY will evaluate the predictive capacity of each
gene individually and will return predictive scores with their confidence
intervals, without using null-hypothesis testing.

To guarantee the reliability of predictive scores, EPCY uses a leave-one-out
cross validation to train multiple Kernel Density Estimation (KDE) classifiers
and evaluate their performances on unseen samples.

.. image:: images/method.png
   :width: 400px
   :alt: Overview of EPCY
   :align: center


EPCY is part of the `Leucegene project <leucegene.ca>`_ and has been developed
and tested to analyse RNAseq data. However, the method implemented in EPCY is
generic and should works on different quantitative data.

We working on a paper, in the meantime more details can be found
`in a poster presented at ISMB ECCB 2019 <https://f1000research.com/posters/8-1349>`_:

Audemard E, Sauvé L and Lemieux S. EPCY: Evaluation of Predictive
CapabilitY for ranking biomarker gene candidates [version 1; not peer reviewed].
*F1000Research 2019*, 8(ISCB Comm J):1349(poster)

**Installation**

PyPI install

.. code:: bash

   pip install epcy
   epcy -h


.. toctree::
   :maxdepth: 2
   :caption: Documentation index

   basic_usage
   rna_tools
   predictive_capability_columns
   small_dataset
