#!/bin/env python

import os
import shutil


def ignore_test_and_tools(dir_, files):
    """Filter files to be moved by shutil.copytree. Ignore any hidden file and the
    test and tools directories, which are not needed by the remote instance.
    :param dir_: dummy variable, must be present to be passed to shutil.copytree()
    :param files: output of os.listdir(), files to be subjected to filtering
    :return list: list of files that should be filtered, and not copied.
    """
    return [
        f
        for f in files
        if (f == "test" or f == "test-data" or f == "__pycache__" or f.startswith("."))
    ]


setup_dir = os.path.dirname(os.path.realpath(__file__))
seqc_dir = os.path.expanduser("~/.seqc/seqc")

print("setup_dir: {}".format(setup_dir))
print("seqc_dir: {}".format(seqc_dir))

# delete the existing one
if os.path.isdir(seqc_dir):
    shutil.rmtree(seqc_dir)

# copy the SEQC files in the working directory to ~/.seqc/seqc
shutil.copytree(setup_dir, seqc_dir, ignore=ignore_test_and_tools)

# create .tar.gz of ~/.seqc/seqc/*
shutil.make_archive(base_name=seqc_dir, format="gztar", root_dir=seqc_dir)
