# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import shutil
import logging
from typing import Union
from .utils import conf, call
from .parameters import Parameters
from .analysis import (
    run_analysis,
    read_data_latt,
    read_data_mat,
    read_green_eq,
    read_green_tau,
    read_greenmat_eq,
    read_greenmat_tau,
)

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
        The start directory to copy. Can be a name of the directories provided by ALF
        or a path to a custom start directory. The default is the `Start` directory
        provided by ALF.
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
    base, name = os.path.split(start_dir)
    if not base:
        src_dir = os.path.join(conf["ALF_DIR"], "Scripts_and_Parameters_files", start_dir)
    else:
        src_dir = start_dir
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
    cmd = os.path.join(conf["ALF_DIR"], "Prog", "ALF.out")
    call(cmd, cwd=directory, verbose=verbose)


class Simulation:
    def __init__(self, directory):
        self.directory = directory
        self.parameters: Union[Parameters, None] = None
        try:
            self.load_parameters()
        except FileNotFoundError:
            pass

    def load_parameters(self):
        self.parameters = Parameters(self.directory)

    def init(self, start_dir="", overwrite=False):
        init_simulation(self.directory, start_dir, overwrite)
        self.load_parameters()

    def run(self, verbose=True):
        run_simulation(self.directory, verbose)

    def out_to_in(self, verbose=True):
        out_to_in(self.directory, verbose)

    def analyze(self, files="*", verbose=True):
        run_analysis(self.directory, files, verbose)

    def update_params(self, data, save=False):
        for sec, items in data.items():
            for key, val in items.items():
                self.parameters.set("var_" + sec, key, val)
        if save:
            self.parameters.save()

    def listdir(self):
        return os.listdir(self.directory)

    def join(self, *args):
        return os.path.join(self.directory, *args)

    def read_obs_latt(self, obs_name, nrebin=None, nskip=None, subtract_back=True):
        return read_data_latt(self.directory, obs_name, nrebin, nskip, subtract_back)

    def read_obs_mat(self, obs_name, nrebin=None, nskip=None, subtract_back=True):
        return read_data_mat(self.directory, obs_name, nrebin, nskip, subtract_back)

    def read_green_eq(self, iorb=0, nrebin=None, nskip=None):
        return read_green_eq(self.directory, iorb, nrebin, nskip)

    def read_green_tau(self, iorb=0, total=False, nrebin=None, nskip=None):
        return read_green_tau(self.directory, iorb, total, nrebin, nskip)

    def read_greenmat_eq(self, iorb=0, nrebin=None, nskip=None):
        return read_greenmat_eq(self.directory, iorb, nrebin, nskip)

    def read_greenmat_tau(self, iorb=0, total=False, nrebin=None, nskip=None):
        return read_greenmat_tau(self.directory, iorb, total, nrebin, nskip)
