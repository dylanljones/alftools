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


def _check_initialized(directory):
    if not os.path.exists(directory):
        return False
    files = list(os.listdir(directory))
    expected = ["parameters", "seeds", "out_to_in.sh"]
    for name in expected:
        if name not in files:
            return False
    return True


def init_simulation(directory, start_dir="", overwrite=False):
    """Initializes an empty ALF simulation directory.

    Equivalent to the terminal command
    .. code-block:: commandline

        cp -r $ALF_DIR/Scripts_and_Parameters_files/<Start> ./<Name>

    which creates the initial simulation directory storing the default simulation
    parameters. After the simulation is run the results will be saved in this
    directory as well.

    Parameters
    ----------
    directory : str
        The path of the ALF output directory that will be created.
    start_dir : str, optional
        The start directory to copy. The default is the `Start` directory provided by
        ALF.
    overwrite : bool, optional
        If True any existing directories with the same name will be deleted before
        initializing the new directoy.
    """
    if not start_dir:
        start_dir = "Start"

    out_dir = os.path.abspath(os.path.normpath(directory))
    if _check_initialized(out_dir):
        if overwrite:
            logger.info("Removing directory %s", out_dir)
            shutil.rmtree(out_dir)
        else:
            logger.info("Directory already exists: %s", out_dir)
            return
    logger.info("Creating initial simulation directory: %s", out_dir)
    src_dir = os.path.join(ALF_DIR, "Scripts_and_Parameters_files", start_dir)
    shutil.copytree(src_dir, out_dir)


def out_to_in(directory, verbose=True):
    """Converts the output configuration to an input configuration.

    Equivalent to the terminal command
    .. code-block:: commandline

        ./out_to_in.sh

    Parameters
    ----------
    directory : str
        The path of an existing ALF simualtion directory.
    verbose : bool, optional
        If True, print the output of the command.
    """
    logger.info("Running 'out_to_in.sh' in '%s'", directory)
    # Run `out_to_in` Script
    call("./out_to_in.sh", cwd=directory, verbose=verbose)


def run_simulation(directory, verbose=True):
    """Runs a new or continued ALF simulation.

    Equivalent to the terminal command
    .. code-block:: commandline

        $ALF_DIR/Prog/ALF.out

    If an existing simulation is continued the method `out_to_in` has to be called
    before running this method.

    Parameters
    ----------
    directory : str
        The path of an initialized ALF simualtion directory.
    verbose : bool, optional
        If True, print the output of the command.
    """
    logger.info("Running simulation in '%s'", directory)
    call(os.path.join(ALF_DIR, "Prog", "ALF.out"), cwd=directory, verbose=verbose)
