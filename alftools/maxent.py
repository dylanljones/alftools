# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import re
import shutil
import logging
from .utils import ALF_DIR, call

logger = logging.getLogger(__name__)


def copy_parameters(src_dir, dst_dir):
    dst = os.path.join(dst_dir, "parameters")
    src = os.path.join(src_dir, "parameters")
    if os.path.exists(dst):
        os.remove(dst)
    logger.debug("Copying parameters to %s", dst_dir)
    shutil.copy(src, dst)


def prepare_green_maxent(root):
    # Copy parameter file to directories
    regex = re.compile(r"Green_\d+\.\d{2}_\d+\.\d{2}")
    dirs = list()
    logger.info("Preparing directory %s for MaxEnt", root)
    for name in os.listdir(root):
        if regex.match(name) is not None:
            green_dir = os.path.join(root, name)
            copy_parameters(root, green_dir)
            dirs.append(green_dir)
    return dirs


def run_maxent(directory, verbose=False):
    logger.info("Running Max_SAC in %s", directory)
    cmd = os.path.join(ALF_DIR, "Analysis", "Max_SAC.out")
    call(cmd, directory, verbose)
