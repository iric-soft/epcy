import sys
import os
import time

import numpy as np

from .. utils import Classifier as uc

#from guppy import hpy
#hp = hpy()

def cut_version(id_str):
    """Remove the version ensembl gene/transcript ids
    """
    index_version = id_str.find('.')
    if index_version != -1:
        id_str = id_str[0:index_version]

    return(id_str)


def compute_pred(args, num_pred, list_ids, data, design):
    sys.stderr.write(time.strftime('%X') + ": Start epcy analysis of " +
                     str(num_pred) + " features\n")

#    print(hp.heap())
    c = uc.Classifier(args, design, data, list_ids)
    del design
    del data
    c.run()
    c.pred2csv()

    sys.stderr.write(time.strftime('%X') + ": End\n")
    return(None)
