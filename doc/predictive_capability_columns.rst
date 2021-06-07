Details of predictive capability columns
========================================

By default, *epcy pred* and *epcy pred_rna* will return an output file
*predictive_capability.xls* of 9 columns:

* Default columns:

  - **id**: the id of each feature.
  - **l2fc**: log2 fold change.
  - **kernel\_mcc**: Matthews Correlation Coefficient (`MCC`_) compute by a predictor using `KDE`_.
  - **kernel\_mcc\_low**, **kernel\_mcc\_high**: lower and upper bounds of confidence interval (90%).
  - **mean\_query**: average values of this feature for samples in the subgroup of interest defined using the --query parameter.
  - **mean\_ref**: average values of this feature for samples in the reference group (not in the query subset).
  - **bw\_query**: estimated bandwidth used by `KDE`_, to calculate the density of query samples.
  - **bw\_ref**: estimated bandwidth used by `KDE`_, to calculate the density of ref samples.

However, epcy can expand or modify this default output using several
options, see:

.. code:: bash

   epcy pred -h


Predictive scores
-----------------

By default we decide to return `MCC`_ scores. However, it's possible to compute
other predictive scores, in case they are more suitable for your needs. Using the following
parameters will add new columns to the default output, as:

* -\-ppv:

  - **kernel\_ppv**: Positive Predictive value (`PPV`_, precision) compute by a predictor using `KDE`_.
  - **kernel\_ppv\_low**, **kernel\_ppv\_high**: boundaries of confidence interval (90%).

* -\-npv:

  - **kernel\_npv**: Negative Predictive value (`NPV`_) compute by a predictor using `KDE`_.
  - **kernel\_npv\_low**, **kernel\_ppv\_high**: boundaries of confidence interval (90%).

* -\-tpr:

  - **kernel\_tpr**: True Positive Rate value (`TPR`_, sensitivity) compute by a predictor using `KDE`_.
  - **kernel\_tpr\_low**, **kernel\_tpr\_high**: boundaries of confidence interval (90%).

* -\-tnr:

  - **kernel\_tnr**: True Negative Rate value (`TNR`_, specificity) compute by a predictor using `KDE`_.
  - **kernel\_tnr\_low**, **kernel\_tnr\_high**: boundaries of confidence interval (90%).

* -\-fnr:

  - **kernel\_fnr**: False Negative Rate value (`FNR`_, miss rate) compute by a predictor using `KDE`_.
  - **kernel\_fnr\_low**, **kernel\_fnr\_high**: boundaries of confidence interval (90%).

* -\-fpr:

  - **kernel\_fpr**: False Positive Rate Rate value (`FPR`_, fall-out) compute by a predictor using `KDE`_.
  - **kernel\_fpr\_low**, **kernel\_fpr\_high**: boundaries of confidence interval (90%).

* -\-fdr:

  - **kernel\_fdr**: False Discovery Rate value (`FDR`_) compute by a predictor using `KDE`_.
  - **kernel\_fdr\_low**, **kernel\_fdr\_high**: boundaries of confidence interval (90%).

* -\-for:

  - **kernel\_for**: False Omission Rate value (`FOR`_) compute by a predictor using `KDE`_.
  - **kernel\_for\_low**, **kernel\_for\_high**: boundaries of confidence interval (90%).

* -\-ts:

  - **kernel\_ts**: Threat Score value (`TS`_, critical sucess index) compute by a predictor using `KDE`_.
  - **kernel\_ts\_low**, **kernel\_ts\_high**: boundaries of confidence interval (90%).

* -\-acc:

  - **kernel\_acc**: Accuracy value (`ACC`_) compute by a predictor using `KDE`_.
  - **kernel\_acc\_low**, **kernel\_acc\_high**: boundaries of confidence interval (90%).

* -\-f1:

  - **kernel\_f1**: F1 score value (`F1`_) compute by a predictor using `KDE`_.
  - **kernel\_f1\_low**, **kernel\_f1\_high**: boundaries of confidence interval (90%).

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

The type of classifier used to evaluate the predictive score of each gene
(feature), is a parameter of EPCY. By default, EPCY will use a `KDE`_
classifier. However, it is possible to replace the `KDE`_ classifier by
a normal classifier, using -\-normal.

Using the normal classifier, all predictive scores (listed above) remain
available. However, the column name of each predictive score, will be changed
to start with *normal* instead of *kernel* (**normal\_mcc** vs **kernel\_mcc**),
to be consistant.

Missing values
--------------

On some dataset (as in proteomics or single-cell), quantitative matrix can have
some missing values (*nan*). In that case there are different alternatives to manage
these missing values within EPCY:
* Impute missing values before running EPCY.
* Replace missing values by a constant, using -\-replacena.
* For each gene (or feature), remove samples with missing values.

If you choose to remove samples with missing values, EPCY will return
a *predictive_capability.xls* with two new columns,
**sample_query** and **sample_ref**, to report for each gene (feature),
the number of query and reference samples used (without missing values).

If you have downloaded the source code or data on `git`_,
you can test these procedures using:

 .. code:: bash

    epcy pred --norm --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/matrix.tsv -o ./data/small_for_test/using_na
    epcy pred --replacena 0 --norm --log -d ./data/small_for_test/design.tsv -m ./data/small_for_test/matrix.tsv -o ./data/small_for_test/replace_na

.. _KDE: https://en.wikipedia.org/wiki/Kernel_density_estimation
.. _MCC: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
.. _PPV: https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
.. _NPV: https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
.. _TPR: https://en.wikipedia.org/wiki/Sensitivity_and_specificity
.. _TNR: https://en.wikipedia.org/wiki/Sensitivity_and_specificity
.. _FNR: https://en.wikipedia.org/wiki/Type_I_and_type_II_errors#False_positive_and_false_negative_rates
.. _FPR: https://en.wikipedia.org/wiki/False_positive_rate
.. _FDR: https://en.wikipedia.org/wiki/False_discovery_rate
.. _FOR: https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
.. _ACC: https://en.wikipedia.org/wiki/Accuracy_and_precision
.. _TS: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
.. _F1: https://en.wikipedia.org/wiki/F-score
.. _normal: https://en.wikipedia.org/wiki/Normal_distribution
.. _MannWhitney: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
.. _ttest\_ind: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
.. _git: https://github.com/iric-soft/epcy/tree/master/data/small_for_test
