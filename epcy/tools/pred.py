import time
import sys

import numpy as np
import pandas as pd

from .. utils import other as uo
from .. utils import readers as ur


def main_pred(args, argparser):

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix " +
                     "features\n")
    (design, data, list_ids) = ur.read_design_matrix(args)

    num_pred = data.shape[0]

    all_classifier = uo.compute_pred(args, num_pred, list_ids, data, design)

    #uo.save_pred_res(args, all_classifier)
