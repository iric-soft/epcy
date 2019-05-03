=============================================================================
EPSY :  Evaluation of Predictive CapabilitY for ranking biomarker candidates
=============================================================================

.. image:: https://travis-ci.org/ericloud/epsy.svg?branch=master
  :target: https://travis-ci.org/ericloud/epsy

---------
Contents:
---------

-------------
Introduction:
-------------

This tool was developed to made differential analysis on large cohort, from `Kallisto`_ quantification.

.. _Kallisto: https://pachterlab.github.io/kallisto/

-------------
Requirements:
-------------

python3

--------
Install:
--------

.. code:: shell

  $ python3 setup.py install
  $ epsy -h

------
Usage:
------

General:
--------

From source:
****************

.. code:: shell

  $ cd [your_kt_folder]
  $ python3 -m epsy -h

After setup install:
********************

.. code:: shell

  $ epsy -h

Small example:
--------------

.. code:: shell

  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup
  $ epcy pred -d ./data/small_for_test/design.tsv -m ./data/small_for_test/exp_matrix.tsv -o ./data/small_for_test/subgroup2 --subgroup subgroup2
