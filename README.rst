=============================================================================
EPCY :  Evaluation of Predictive CapabilitY for ranking biomarker candidates
=============================================================================

+------------------------------------------------------------------+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| .. image:: https://img.shields.io/badge/python-3.6-blue.svg      | .. image:: https://travis-ci.org/iric-soft/epcy.svg?branch=master | .. image:: https://codecov.io/gh/iric-soft/epcy/branch/master/graph/badge.svg |
|    :target: https://www.python.org/downloads/release/python-362/ |    :target: https://travis-ci.org/iric-soft/epcy                  |    :target: https://codecov.io/gh/iric-soft/epcy/                             |
+------------------------------------------------------------------+-------------------------------------------------------------------+-------------------------------------------------------------------------------+


-------
Citing:
-------
* EPCY: Evaluation of Predictive CapabilitY for ranking biomarker gene candidates. Poster at ISMB ECCB 2019: https://f1000research.com/posters/8-1349

-------------
Introduction:
-------------

This tool was developed to Evaluate Predictive CapabilitY of each feature to become a biomarker candidates.

-------------
Requirements:
-------------

* python3
* (Optional) virtualenv

--------
Install:
--------

.. code:: shell

  python3 -m venv $HOME/.virtualenvs/epcy
  source $HOME/.virtualenvs/epcy/bin/activate
  cd [your_epcy_folder]
  CFLAGS=-std=c99 pip3 install numpy==1.17.0
  python3 setup.py install
  epcy -h

------
Usage:
------

General:
--------

From source:
****************

.. code:: shell

  cd [your_epcy_folder]
  python3 -m epcy -h

After setup install:
********************

.. code:: shell

  epcy -h

Small example:
--------------

.. code:: shell

  cd [your_epcy_folder]
  # Run epcy using default parameter on normalized expression matrix (like TPM, FPKM or RPKM)
  epcy pred_rna -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup
  # Run epcy using default parameter on readcount not normalized add --cpm
  epcy pred_rna -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup --cpm
  # Run epcy without filter, this time you will have predictive analysis for all features (no NA in the output)
  epcy pred_rna -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/no_filter_subgroup -l 0
  # Run epcy on the second design (column subgroup2) describe in ./data/small_for_test/design.tsv
  epcy pred_rna -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup2 --subgroup subgroup2

  # Run epcy on any normalized quantification data
  epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup

Output:
-------

Prediction\_capability.xls: the main output which contain the evaluation of each features (genes, proteins, ...). It's a tabulated files 9 columns:

 * Default columns:

   - id: the id of each feature.
   - l2fc: log2 Fold change.
   - kernel\_mcc: Matthews Correlation Coefficient (`MCC`_) compute by a predictor using `KDE`_.
   - mean\_query: mean(values) of samples specify as Query in design.tsv
   - mean\_ref: mean(values) of samples specify as Ref in design.ts

 * Using --normal:

   - normal\_mcc: `MCC`_ compute a predictor using `normal`_ distributions.

 * Using --auc --utest:

   - auc: Area Under the Curve
   - u\_pv: pvalue compute by a `MannWhitney`_ rank test

 * Using --ttest:

   - t\_pv: pvalue compute by `ttest\_ind`_


Using --full a secondary output file (subgroup\_predicted.xls) specify for each features if the sample as been correctly predict. Build an heatmap with this output could help you to explore your data.

 * Legend:

   - 1: true positive
   - 2: false negative
   - 3: false positive
   - 4: true negative

   .. _MCC: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
   .. _KDE: https://en.wikipedia.org/wiki/Kernel_density_estimation
   .. _normal: https://en.wikipedia.org/wiki/Normal_distribution
   .. _MannWhitney: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
   .. _ttest\_ind: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
   .. _contingency: https://en.wikipedia.org/wiki/Confusion_matrix
