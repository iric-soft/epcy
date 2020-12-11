How to use EPCY
===============

We have implemented EPCY in python3 as a software, which allow you to perform
predictive analyses using bash command line.

If it's not already done, you need to install EPCY

.. code:: bash

   pip install epcy
   epcy -h

Input files
-----------

EPCY is design to work on severals quantitative data (like genes expression),
provided that data are normalized: quantitative values of each genes (features)
need to be comparable between each samples.

To run EPCY, we need two **tabulated** files, as input:

  * A **matrix of quantitative data**, normalized for each samples (in columns)
    with an *ID* column to identify each genes (features), in first.

  .. list-table:: Example of a quantitative matrix
       :widths: 30 20 20 20
       :header-rows: 1

       * - ID
         - Sample1
         - ...
         - SampleX
       * - Gene1
         - 10
         - ...
         - 20
       * - ...
         - ...
         - ...
         - ...
       * - GeneY
         - 12
         - ...
         - 16

  * A **design table** which describe the comparison. This table is composed
    by a first column *Sample*, followed by at least one column to assign
    observed labels for each samples. A new column is needed for each
    conditions.

  .. list-table:: Example of a design table
     :widths: 30 20 20
     :header-rows: 1

     * - Sample
       - Condition1
       - Gender
     * - Sample1
       - Query
       - M
     * - ...
       - ...
       - ...
     * - SampleX
       - Ref
       - F


Download input files
--------------------

For the tutorial, we propose to download a part of the
Leucegene cohort composed by 100 Acute Myeloid Leukemia (AML) individual
samples. To reduce the execution time, we are going to analyze only genes that
code for a protein (19,892 genes).
These input file are available in EPCY source code, using *git*:

.. code:: bash

   git clone
   cd epcy/data/leucegene3

If you take a look of *design.txt*, you can see an *AML_subtype* column which
classify each samples in one of these 3 subtypes of AML:
*t15_17*, *inv16* and *other*.

  .. list-table:: design.txt
     :widths: 30 20
     :header-rows: 1

     * - Sample
       - AML_subtype
     * - 01H001
       - Other
     * - 01H002
       - t15_17
     * - ...
       - ...
     * - 14H103
       - inv16
     * - 14H133
       - t15_17

We will start by analyse *t15_17* samples versus all other samples (*inv16* and
*other*). On a macbook pro 2 GHz Dual-Core Intel Core i5, this analysis take
8 min, using 4 thread.

Run your first EPCY analysis
----------------------------

EPCY is divide in severals tools, which can be listed using:

.. code:: bash

   epcy -h

Among all these tools, *epcy pred* is the one which allow to run a default
comparative predictive analysis, as follow:

.. code:: bash

   epcy pred --log -t 4 -m cpm.xls  -d design.txt --subgroup AML_subtype --query t15_17 -o ./30_t15_17_vs_70/

Which:
  * **-\-log**: specify that quantitative data need to to be log transformed
    before to be analyzed.
  * **-t 4**: to allow to run run on 4 threads and divide the execution time by 4.
  * **-m cpm.xls**: to get the path of quantitative matrix file.
  * **-d design.txt**: to get the path of the design table.
  * **-\-subgroup AML_subtype**: to specify the condition column we want use.
  * **-\-query t15_17**: to specify which type of samples we want to compare
                         to each other samples.
  * **-o ./30_t15_17_vs_70/**: to specify the output directory.

More information can be found, using *epcy pred -h*.

If everything is correct you should see appear 4 lines, similar to:

.. code:: bash

   16:31:24: Read design and matrix features
   16:31:34: Start epcy analysis of 19892 features
   16:39:48: Save epcy results
   16:39:49: End

Results
-------

**predictive_capability.xls** is the main output of EPCY analysis, which contain the evaluation
of each genes (features). It's a tabulated files of 9 columns:

* **id**: the id of each feature.
* **l2fc**: log2 fold change.
* **kernel\_mcc**: Matthews Correlation Coefficient (`MCC`_) compute by a predictor using `KDE`_.
* **kernel\_mcc\_low**, **kernel\_mcc\_high**: boundaries of confidence interval (90%).
* **mean\_query**: mean(values) of samples specify as Query in design.tsv
* **mean\_ref**: mean(values) of samples specify as Ref in design.ts
* **bw\_query**: Estimate bandwidth used by `KDE`_, to calculate the density of query samples
* **bw\_ref**: Estimate bandwidth used by `KDE`_, to calculate the density of ref samples

  .. _MCC: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
  .. _KDE: https://en.wikipedia.org/wiki/Kernel_density_estimation

It remains to order this table *kernel_mcc*, to rank most predictive genes.

.. list-table:: ./30_t15_17_vs_70/predictive_capability.xls ordered on kernel_mcc
   :widths: 30 10 15 20 20 15 15 15 15
   :header-rows: 1

   * - id
     - l2fc
     - kernel_mcc
     - kernel_mcc_low
     - kernel_mcc_high
     - mean_query
     - mean_ref
     - bw_query
     - bw_ref
   * - ENSG00000089820.15
     - -4.30
     - 0.96
     - 0.51
     - 0.97
     - 4.23
     - 8.53
     - 0.43
     - 0.22
   * - ENSG00000173531.15
     - 3.23
     - 0.92
     - 0.58
     - 0.97
     - 6.22
     - 2.99
     - 0.52
     - 0.21
   * - ENSG00000168004.9
     - 3.64
     - 0.91
     - 0.82
     - 0.95
     - 3.90
     - 0.26
     - 0.29
     - 0.10
   * - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...

Some details on the design table
--------------------------------

As mentioned before, *design.txt* classify samples in 3 different
subtypes (*t15_17*, *inv16* and *other*). Similarly as we did for *t15_17*, we
can analyse *inv16* samples vs all others samples (*t15_17* and
*other*), using the command below:

.. code:: bash

   epcy pred --log -t 4 -m cpm.xls  -d design.txt --subgroup AML --query inv16 -o ./30_inv16_vs_70/


Moreover, it is possible to add a column in **design.txt** for each conditions
you want to compare. Indeed, with the design table given as example
(in introduction), we could make an analyse on **Gender**,
using *-\-subgroup Gender -\-query M -o ./gender*.

Also, if some annotation are unknown for some samples, we can removed these
samples from the analysis, using **None** in cells which correspond.

  .. list-table:: Example where the AML subtype of sampleX is unknown and
                  need to be removed from the analysis.
     :widths: 30 20 20
     :header-rows: 1

     * - Sample
       - AML_subtype
       - Gender
     * - Sample1
       - t15_17
       - M
     * - ...
       - ...
       - ...
     * - SampleX
       - None
       - F

Using all these variations, you should be able to create a unique *design.txt*
to performed all *1 vs all* analyses linked to a dataset.
