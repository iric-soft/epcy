Details of predictive capability columns
========================================

By default, *epcy pred* and *epcy pred_rna* will retun an output file
*predictive_capability.xls* of 9 columns:

* Default columns:

  - **id**: the id of each feature.
  - **l2fc**: log2 fold change.
  - **kernel\_mcc**: Matthews Correlation Coefficient (`MCC`_) compute by a predictor using `KDE`_.
  - **kernel\_mcc\_low**, **kernel\_mcc\_high**: boundaries of confidence interval (90%).
  - **mean\_query**: mean(values) of samples specify as Query in design.tsv
  - **mean\_ref**: mean(values) of samples specify as Ref in design.ts
  - **bw\_query**: Estimate bandwidth used by `KDE`_, to calculate the density of query samples
  - **bw\_ref**: Estimate bandwidth used by `KDE`_, to calculate the density of ref samples

However, epcy can expand or modify this default output using several options,
see:

.. code:: bash

   epcy pred -h


Predictive scores
-----------------

By default we decide to return `MCC`_ scores. However, it's possible to compute
other predictive scores, in case they are more suitable for your needs. Use
these parameters will add new columns to the default output, as:

* -\-ppv:

  - **kernel\_ppv**: Positive Predictive Value (`PPV`_) compute by a predictor using `KDE`_.
  - **kernel\_ppv\_low**, **kernel\_ppv\_high**: boundaries of confidence interval (90%).

* -\-npv:

  - **kernel\_npv**: Negative Predictive Value (`NPV`_) compute by a predictor using `KDE`_.
  - **kernel\_npv\_low**, **kernel\_ppv\_high**: boundaries of confidence interval (90%).


* -\-auc:

  - **auc**: Area Under the Curve

Statistical test
----------------

EPCY is able to perform statistical tests, using:

* -\-auc -\-utest:

  - **u\_pv**: pvalue compute by a `MannWhitney`_ rank test

* -\-ttest:

  - **t\_pv**: pvalue compute by `ttest\_ind`_


Normal distribution
-------------------

EPCY is able to evaluate predictive capacity using a normal distribution
(in place of `KDE`_), using:

* -\-normal:

  - **normal\_mcc**: `MCC`_ compute a predictor using `normal`_ distributions.

.. _KDE: https://en.wikipedia.org/wiki/Kernel_density_estimation
.. _MCC: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
.. _PPV: https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
.. _NPV: https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
.. _normal: https://en.wikipedia.org/wiki/Normal_distribution
.. _MannWhitney: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
.. _ttest\_ind: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
