
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

  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup
  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup -l 0
  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup2 --subgroup subgroup2
