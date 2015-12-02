#!/usr/local/bin/python3
__author__ = 'Ambrose J. Carr'

import argparse
import os
import seqc
from copy import copy
import sys
import json


if __name__ == "__main__":
    parser = seqc.core.create_parser()
    kwargs = seqc.core.parse_args(parser)
    seqc.log.setup_logger()
    try:
        # log command line arguments for debugging
        arg_copy = copy(kwargs)
        del arg_copy['func']  # function is not serializable
        seqc.log.info('SEQC version: %s' % seqc.__version__)
        seqc.log.info('SEQC working directory: %s' % os.getcwd())
        seqc.log.info('Passed command line arguments: %s' %
                      json.dumps(arg_copy, separators=(',', ': '), indent=4))

        func = kwargs['func']
        func(**kwargs)

    except:
        seqc.log.exception()
        raise
