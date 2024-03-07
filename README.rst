=============================================================================
EPCY :  Evaluation of Predictive CapabilitY for ranking biomarker candidates
=============================================================================

+------------------------------------------------------------+------------------------------------------------------------------+
| .. image:: https://zenodo.org/badge/197271057.svg          | .. image:: https://img.shields.io/badge/python-3.11.5-blue.svg      |
|    :target: https://zenodo.org/doi/10.5281/zenodo.10407905 |    :target: https://www.python.org/downloads/release/python-3115/|
+------------------------------------------------------------+------------------------------------------------------------------+


-------
Citing:
-------
* Manuscript in preparation
* EPCY: Evaluation of Predictive CapabilitY for ranking biomarker gene candidates. Poster at ISMB ECCB 2019: https://f1000research.com/posters/8-1349

-------------
Introduction:
-------------

This tool was developed to Evaluate Predictive CapabilitY of each gene (feature) to become a predictive (bio)marker candidates.
Documentation is available `via Read the Docs <https://epcy.readthedocs.io/>`_.

-------------
Requirements:
-------------

* python >= 3.11.5

--------
Install:
--------

Using pypi:
-----------

.. code:: shell

  pip install epcy

From source:
------------
.. code:: shell

  python3 -m venv $HOME/.virtualenvs/epcy
  source $HOME/.virtualenvs/epcy/bin/activate
  pip install pip setuptools --upgrade
  pip install wheel
  cd [your_epcy_folder]
  pip install -e .
  epcy -h

------
Usage:
------

General:
--------

After install:
**************

.. code:: shell

  epcy -h

From source:
************

.. code:: shell

  cd [your_epcy_folder]
  python3 -m epcy -h

Generic case:
-------------

* EPCY is design to work on any quantitative data, provided that values of each feature are comparable between each samples (normalized).
* To run a comparative analysis, `epcy pred` need two tabulated files:

  * A `matrix`_ of quantitative normalized data for each samples (column) with an "ID" column to identify each feature.
  * A `design`_ table which describe the comparison.

.. _matrix: https://github.com/iric-soft/epcy/blob/master/data/small_for_test/normalized_matrix.tsv
.. _design: https://github.com/iric-soft/epcy/blob/master/data/small_for_test/design.tsv

.. code:: shell

  # Run epcy on any normalized quantification data
  epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/log_normalized_matrix.tsv -o ./data/small_for_test/EPCY_output

  # If your data are normalized, but require a log2 transforamtion, add --log
  epcy pred --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/normalized_matrix.tsv -o ./data/small_for_test/EPCY_output

  # If your data are not normalized and require a log2 transforamtion, add --norm --log
  epcy pred --norm --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/matrix.tsv -o ./data/small_for_test/EPCY_output

  # Different runs might show small variations.
  # To ensure reproducibility set a random seed, using --randomseed
  epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/normalized_matrix.tsv -o ./data/small_for_test/EPCY_output --randomseed 42
  epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/normalized_matrix.tsv -o ./data/small_for_test/EPCY_output2 --randomseed 42
  diff ./data/small_for_test/EPCY_output/predictive_capability.tsv ./data/small_for_test/EPCY_output2/predictive_capability.tsv


More documentation is available `via Read the Docs <https://epcy.readthedocs.io/>`_.
