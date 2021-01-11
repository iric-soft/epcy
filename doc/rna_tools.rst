Working on RNA quantification
=============================

EPCY is a generic method, which can be used to analyse several types of
quantification data like statistical tests, compare to Differential Expression
methods specifically design to analyse RNA quantification data. However, EPCY
have specific tools, *pred_rna* and *profile_rna*, to simplify the analysis
of RNA data.

bulk RNA
--------

In the first steps section, we saw how to run EPCY on normalized quantitative
data. However, EPCY can work directly on read counts, using *pred_rna* tool:

.. code:: bash

   # The commande seen in first steps sections
   epcy pred --log -t 4 -m cpm.xls  -d design.txt --subgroup AML --query t15_17 -o ./30_t15_17_vs_70/ --randomseed 42
   # Equivalent analysis using read counts quantification and pred_rna
   epcy pred_rna --cpm --log -t 4 -m readcounts.xls  -d design.txt --subgroup AML --query t15_17 -o ./30_t15_17_vs_70_readcounts/ --randomseed 42

Similarly, you can use *profile_rna* to plot a learned KDE with his
gene expression, directly from read counts:

.. code:: bash

  # The commande seen in first steps sections
  epcy profile --log -m cpm.xls -d design.txt --subgroup AML --query t15_17 -o ./30_t15_17_vs_70/figures/ --ids ENSG00000162493.16 ENSG00000227268.4

  # Equivalent commande using read counts and profile_rna
  epcy profile_rna --cpm --log -m readcounts.xls -d design.txt --subgroup AML --query t15_17 -o ./30_t15_17_vs_70_readcounts/figures/ --ids ENSG00000162493.16 ENSG00000227268.4


.. image:: images/profile_readcounts.png
   :width: 400px
   :alt: readcounts gene profiles
   :align: center

More help can be found using:

.. code:: bash

  epcy pred_rna -h
  epcy profile_rna -h

Kallisto quantification
-----------------------

EPCY allow to work directly on kallisto quantificaion using h5 files and have
access to bootstrapped samples. To do so, a `kallisto` column need to be added
to the design file (to specify the directory path of *abundant.h5*
file for each sample).

Using data available on `git`_ and **epcy pred_rna**, you can run EPCY
as follow:

.. code:: shell

  # To run on kallisto quantification, add --kall --cpm --log
  epcy pred_rna --kal --cpm --log -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5/trans/
  # Take care:
  # - kallisto quantification is on transcript, not on gene
  # - On real (complet) dataset, it is recommended to add some threads (-t)

To run on gene level, a gff3 file of the genome annotation need to be specify,
to have the correspondence between transcript and gene. This file can be
download on `ensembl`_ and added in the command, as follow:

.. code:: shell

  # To run on kallisto quantification and gene level, add --gene --anno [file.gff]
  epcy pred_rna --kal --cpm --log --gene --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5/gene/ --randomseed 42
  epcy profile_rna --kal --cpm --log --gene --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5/gene/figures --ids ENSG00000100345
  # If you prefer analyse your data on tpm, replace --cpm by --tpm

To take account of inferential variance (introduce by `sleuth`_), EPCY can use
bootstrapped samples, using -\-bs:

.. code:: shell

  epcy pred_rna --kal --cpm --log --gene --bs 10 --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5_bs/gene/ --randomseed 42
  epcy profile_rna --kal --cpm --log --gene --bs 10 --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5_bs/gene/figures --ids ENSG00000100345

When reading all kallisto files is time consuming, you can use *epcy kal2mat*
tool, to create a quantification matrix file and use EPCY, as before:

.. code:: shell

  # Without bootstrapped samples
  epcy kal2mat --gene --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5_mat/gene/
  epcy pred_rna --cpm --log -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -m ./data/small_leucegene/5_inv16_vs_5_mat/gene/readcounts.xls -o ./data/small_leucegene/5_inv16_vs_5_mat/gene/ --randomseed 42
  epcy profile_rna --cpm --log -m ./data/small_leucegene/5_inv16_vs_5_mat/gene/readcounts.xls -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5_mat/gene/figures --ids ENSG00000100345

  # With bootstrapped samples
  epcy kal2mat --gene --bs 10 --anno ./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3 -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5_mat_bs/gene/
  epcy pred_rna --bs 10 --cpm --log -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -m ./data/small_leucegene/5_inv16_vs_5_mat_bs/gene/readcounts.xls -o ./data/small_leucegene/5_inv16_vs_5_mat_bs/gene/ --randomseed 42
  epcy profile_rna --bs 10 --cpm --log -m ./data/small_leucegene/5_inv16_vs_5_mat_bs/gene/readcounts.xls -d ./data/small_leucegene/5_inv16_vs_5/design.tsv -o ./data/small_leucegene/5_inv16_vs_5_mat_bs/gene/figures --ids ENSG00000100345

Single-cell
-----------

Several developments are planned, to get EPCY more user friendly on
single-cell data (to manage parse matrix, run on GPU, ...). Pending, you can
analyses your single-cell data with *epcy pred* and *epcy profile*, on
normalized expression data (like in first steps).

On read counts (not normalized), you can use *epcy pred_rna* and
*epcy profile_rna* with -\-cpmed (in place of -\-cpm) to normalized read
counts according to median depth of the dataset.

.. code:: shell

  epcy pred_rna --cpmed --log ...


.. _git: https://github.com/iric-soft/epcy/tree/master/data/small_leucegene/5_inv16_vs_5/
.. _ensembl: https://useast.ensembl.org/info/data/ftp/index.html
.. _sleuth: https://www.nature.com/articles/nmeth.4324?WT.feed_name=subjects_gene-expression#Sec1
