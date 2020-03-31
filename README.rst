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

Generic case:
--------------

* EPCY is design to work on any quantitative data, provided that values of each feature are comparable between each samples (normalized).
* To run a comparative analysis, `epcy pred` need two tabulated files:

  * A `matrix`_ of quantitative normalized data for each samples (column) with an "ID" column to identify each feature.
  * A `design`_ table which describe the comparison.

.. _matrix: https://github.com/iric-soft/epcy/blob/master/data/small_for_test/exp_matrix.tsv
.. _design: https://github.com/iric-soft/epcy/blob/master/data/small_for_test/design.tsv

.. code:: shell

  # Run epcy on any normalized quantification data
  epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup
  # If your data require a log2 transforamtion, add --log
  epcy pred --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup

* Result will be saved in prediction\_capability.xls file, which is detail below.
* You can personalize the design file using **--subgroup --query**

.. code:: shell

  epcy pred_rna -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup2 --subgroup subgroup2 --query A


Working on RNA sequencing readcounts:
-------------------------------------

* To run EPCY on readcounts not mormalized use **pred_rna** tool as follow:

.. code:: shell

  # To run on read count not normalized, add --cpm --log
  epcy pred_rna --cpm --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup

Working on kallisto quantification:
-----------------------------------

* EPCY allow to work directly on kallisto quantificaion using h5 files, to have access to bootstrapped samples. To do so, a `kallisto` column need to be add to the design file (to specify the directory path where to find *abundant.h5* file for each sample) and **epcy pred_rna** need to run as follow:

.. code:: shell

  # To run on kallisto quantification, add --kall (+ --cpm --log)
  epcy pred_rna --kal --cpm --log -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5/
  # !!! Take care kallisto quantification is on transcript not on gene

* To run on gene level, a gff3 file of the genome annotation is needed, to have the correspondence between transcript and gene. This file can be download on `ensembl`_

.. _ensembl: https://useast.ensembl.org/info/data/ftp/index.html

.. code:: shell

  # To run on kallisto quantification and gene level, add --gene --anno [file.gff] (+ --kall --cpm --log)
  epcy pred_rna --kal --cpm --log --gene --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5/

* kallisto quantification allow to work on TPM:

.. code:: shell

  # work on TPM, replace --cpm by --tpm
  epcy pred_rna --kal --tpm --log --gene --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5/


-------
Output:
-------

predictive\_capability.xls
---------------------------

This file is the main output which contain the evaluation of each features (genes, proteins, ...). It's a tabulated files 9 columns:

* Default columns:

  - id: the id of each feature.
  - l2fc: log2 Fold change.
  - kernel\_mcc: Matthews Correlation Coefficient (`MCC`_) compute by a predictor using `KDE`_.
  - kernel\_mcc\_low, kernel\_mcc\_high: boundaries of confidence interval (90%).
  - mean\_query: mean(values) of samples specify as Query in design.tsv
  - mean\_ref: mean(values) of samples specify as Ref in design.ts
  - bw\_query: Estimate bandwidth used by `KDE`_, to calculate the density of query samples
  - bw\_ref: Estimate bandwidth used by `KDE`_, to calculate the density of ref samples

* Using --normal:

  - normal\_mcc: `MCC`_ compute a predictor using `normal`_ distributions.

* Using --auc --utest:

  - auc: Area Under the Curve
  - u\_pv: pvalue compute by a `MannWhitney`_ rank test

* Using --ttest:

  - t\_pv: pvalue compute by `ttest\_ind`_


subgroup\_predicted.xls
-----------------------

Using --full a secondary output file (subgroup\_predicted.xls) specify for each features if the sample as been correctly predicted. Build an heatmap with this output could help you to explore your data.
More details coming soon.

  .. _MCC: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
  .. _KDE: https://en.wikipedia.org/wiki/Kernel_density_estimation
  .. _normal: https://en.wikipedia.org/wiki/Normal_distribution
  .. _MannWhitney: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
  .. _ttest\_ind: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
  .. _contingency: https://en.wikipedia.org/wiki/Confusion_matrix

--------
Bagging:
--------

To improve the stability and accuracy of MCC computed, you can add n `bagging`_ (using `-b n`)

.. code:: shell

  #Take care, it's take n time more longer!!!, use multiprocess (-t) seems a good idea :).
  epcy pred_rna -b 4 -t 4 --cpm --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup


.. _bagging: https://en.wikipedia.org/wiki/Bootstrap_aggregating
