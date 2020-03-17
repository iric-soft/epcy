from .common import *

from .get_rna_norm_part import *
from .get_kallisto_part import *
from .get_gene_part import *
from .profile import *


def get_argparser_profile_rna(parser):

    get_argparser_profile(parser)
    get_argparser_rna_norm_part(parser)
    get_argparser_kallisto_part(parser)
    get_argparser_gene_part(parser)
