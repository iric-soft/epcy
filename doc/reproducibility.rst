Reproducibility
===============

EPCY draw a random value to assign a class according probabilities learn by
the KDE classifier, to fill a contingency table (see algorithme section).
This means that different runs of EPCY can produce different results.

However, EPCY output is relatively stable, as each predictive score returned
is a mean on several predictive scores (by default 100), to minimize the
variance between runs. Nevertheless, different runs may have small variations.
To ensure the reproducibility, we add a parameter to fix the random seed,
using **-\-randomseed**.

Here an example on the dataset used for the tutorial (see, How to use EPCY).

.. code:: bash

   epcy pred --randomseed 42 --log -t 4 -m cpm.xls  -d design.txt --subgroup AML --query inv16 -o ./30_inv16_vs_70/
