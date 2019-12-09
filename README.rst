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

* Result will be saved in prediction\_capability.xls file, which is detail here.
* It is possible to change `subgroup` column name and the name of each subgroup using `--subgroup` `--query`

Working on RNA sequencing readcounts:
-------------------------------------

* EPCY allow to work directly on readcounts not mormalized, using `epcy pred_rna` as follow

.. code:: shell

  # To run on read count not normalized, add --cpm --log
  epcy pred_rna --cpm --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/default_subgroup

Working on kallisto quantification:
-----------------------------------

* EPCY allow to work directly kallisto quantificaion using h5 files, to have access to bootstrapped samples. To do so, a `kallisto` column need to be add to the design file (to specify the directory path where to find `abundant.h5` file for each sample) and `epcy pred_rna` need to run as follow:

.. code:: shell

  # To run on kallisto quantification, add --kall (+ --cpm --log)
  epcy pred_rna --kal --cpm --log -d [design.tsv] -o [output_directory]
  # !!! Take care kallisto quantification is on transcript not on gene

* To run on gene level, a gff3 file of the genome annotation is needed, to have the correspondence between transcript and gene. This file can be download on `ensenbl`_

.. code:: shell

  # To run on kallisto quantification and gene level, add --gene --anno [file.gff] (+ --kall --cpm --log)
  epcy pred_rna --kal --cpm --log --gene --anno [gff.file] -d [design.tsv] -o [output_directory]

* kallisto quantification allow to work on TPM:

.. code:: shell

  # work on TPM, replace --cpm by --tpm
  epcy pred_rna --kal --tpm --log --gene --anno [gff.file] -d [design.tsv] -o [output_directory]



.. _ensembl: https://useast.ensembl.org/info/data/ftp/index.html
