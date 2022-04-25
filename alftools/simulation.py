# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import shutil
import logging
import subprocess
from .config import conf
from .parameters import Parameters

logger = logging.getLogger(__name__)

ALF_DIR = conf["alf_dir"]


def get_simulation_path(name, root=""):
    if not root:
        root = os.path.dirname(os.path.dirname(ALF_DIR))
    return os.path.join(root, name)


def init_simulation_directory(name, root="", overwrite=False):
    """Initializes an empty ALF simulation directory.

    Equivalent to the terminal command

    .. code-block:: commandline

        cp -r $ALF_DIR/Scripts_and_Parameters_files/Start ./<Name> && cd ./<Name>

    which creates the initial simulation directory storing the default simulation
    parameters. After the simulation is run the results will be saved in this
    directory as well.

    Parameters
    ----------
    name : str
        The name of the ALF output directory that will be created.
    root : str, optional
        The directory where the ALF output directory will be created. by default, the
        root directory of ALF itself is used.
    overwrite : bool, optional
        If True any existing directories with the same name will be deleted before
        initializing the new directoy.

    Returns
    -------
    out_dir : str
        The path of the newly created ALF output directory.
    """
    if not root:
        root = os.path.dirname(os.path.dirname(ALF_DIR))
    out_dir = os.path.join(root, name)
    if overwrite and os.path.exists(out_dir):
        logger.info("Removing directory %s", out_dir)
        shutil.rmtree(out_dir)
    logger.info("Creating initial simulation directory: %s", out_dir)
    src_dir = os.path.join(ALF_DIR, "Scripts_and_Parameters_files", "Start")
    shutil.copytree(src_dir, out_dir)
    return out_dir


def init_simulation(name, root="", overwrite=False):
    out_dir = init_simulation_directory(name, root, overwrite)
    params = Parameters(out_dir)
    return out_dir, params


def load_simulation(name, root=""):
    out_dir = get_simulation_path(name, root)
    params = Parameters(out_dir)
    return out_dir, params


def run_simulation(out_dir=""):
    print("Starting simulation")
    print("-" * 90)
    alf_prog = os.path.join(ALF_DIR, "Prog", "ALF.out")
    subprocess.call(
        [
            alf_prog,
        ],
        cwd=out_dir,
    )
    print("-" * 90)
    print("Finished!")
