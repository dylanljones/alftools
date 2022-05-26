# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import shutil
import logging
from .utils import ALF_DIR, call

logger = logging.getLogger(__name__)


def copy_parameters(src_dir, dst_dir, overwrite=False):
    dst = os.path.join(dst_dir, "parameters")
    src = os.path.join(src_dir, "parameters")
    if os.path.exists(dst):
        if overwrite:
            os.remove(dst)
        else:
            return
    logger.debug("Copying parameters to %s", dst_dir)
    shutil.copy(src, dst)


def run_maxent(directory, name, verbose=False):
    path = os.path.join(directory, name)
    # Prepare directory
    src = os.path.join(directory, "parameters")
    dst = os.path.join(path, "parameters")
    if not os.path.exists(dst):
        logger.debug("Copying parameters to %s", path)
        shutil.copy(src, dst)

    # Run MaxEnt
    logger.info("Running Max_SAC in %s", path)
    cmd = os.path.join(ALF_DIR, "Analysis", "Max_SAC.out")
    call(cmd, path, verbose)
