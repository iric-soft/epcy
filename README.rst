=============================================================================
EPCY :  Evaluation of Predictive CapabilitY for ranking biomarker candidates
=============================================================================

.. image:: https://travis-ci.org/ericloud/epsy.svg?branch=master
  :target: https://travis-ci.org/ericloud/epsy

---------
Contents:
---------

-------------
Introduction:
-------------

This tool was developed to Evaluate Predictive CapabilitY of each feature to become a biomarker candidates.

-------------
Requirements:
-------------

python3

--------
Install:
--------

.. code:: shell

  $ virtualenv $HOME/.virtualenvs/epcy
  $ source $HOME/.virtualenvs/epcy/bin/activate
  $ cd [your_epcy_folder]
  $ python3 setup.py install
  $ epcy -h

------
Usage:
------

General:
--------

From source:
****************

.. code:: shell

  $ cd [your_epcy_folder]
  $ python3 -m epcy -h

After setup install:
********************

.. code:: shell

  $ epcy -h

Small example:
--------------

.. code:: shell

  $ cd [your_epcy_folder]
  # Run epcy using default parameter
  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup
  # Run epcy without filter, this time you will have predictive analysis for all features (no NA in the output)
  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/no_filter_subgroup -l 0
  # Run epcy on the second design (column subgroup2) describe in ./data/small_for_test/design.tsv
  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup2 --subgroup subgroup2

Output:
-------

EPCY's output have 2 files :
 * prediction\_capability.xls: the main output which contain the evaluation of each features (genes, proteins, ...). It's a tabulated files 9 columns:
 
   - ID: the ID of each feature.
   - L2FC: log2 Fold change.
   - KENEL\_MCC: Matthews Correlation Coefficient (`MCC`_) compute by a predictorusing `KDE`_. 
   - NORMAL\_MCC: `MCC`_ compute a predictor using 2 `normal`_ distributions.
   - AUC: Area Under the Curve
   - MEAN\_QUERY: mean(values) of samples specify as Query in design.tsv
   - MEAN\_REF: mean(values) of samples specify as Ref in design.ts
   - U\_PV: pvalue compute by a `MannWhitney`_ rank test
   - T\_PV: pvalue compute by `ttest\_ind`_
   

 * subgroup\_predicted.xls: This secondary output specify for each features if the sample as been correctly predict to feed the `contingency`_ table use to compute KERNEL\_MCC. Build an heatmap with this output could help you to explore your data.
 
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
 