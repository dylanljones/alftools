# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import re
import logging
import itertools
import numpy as np
from .utils import (
    ALF_DIR,
    ComplexParseError,
    BinHeaderError,
    ParseError,
    strings_to_numbers,
    csv_to_complex,
    call,
)

logger = logging.getLogger(__name__)


def contains_analysis(directory):
    for name in os.listdir(directory):
        path = os.path.join(directory, name)
        if os.path.isdir(path):
            return True
    return False


def get_kspace_dirs(directory, obs_name):
    regex = re.compile(obs_name + r"_\d+\.\d{2}_\d+\.\d{2}")
    names = list()
    for name in os.listdir(directory):
        if regex.match(name):
            names.append(name)
    return names


def run_analysis(directory, files="*", verbose=True):
    """Runs the ALF analysis program.

    Equivalent to the terminal command
    .. code-block:: commandline

        $ALF_DIR/Analysis/ana.out <files>

    Parameters
    ----------
    directory : str
        The path of an initialized ALF simualtion directory.
    files : str, optional
        The output files to analyze. By default, all files are analyzed.
    verbose : bool, optional
        If True, print the output of the command.
    """
    if os.path.exists(os.path.join(directory, "data.h5")):
        cmd = os.path.join(ALF_DIR, "Analysis", "ana_hdf5.out")
        if files != "*":
            cmd += " " + files
    else:
        cmd = os.path.join(ALF_DIR, "Analysis", "ana.out") + " " + files
    logger.info("Running Analysis in %s: %s", directory, cmd)
    call(cmd, cwd=os.path.normpath(directory), verbose=verbose)


def read_info(directory: str, name: str, key_map: dict = None):
    """Parses an ALF output info file

    Parameters
    ----------
    directory : str
        The path of the ALF simulation output directory.
    name : str
        The name of the info file. Must end with `_info`.
    key_map: dict, optional
        A dictionary for renaming additional keys of the info file. By default,
        'number of orbitals' is renamed to 'norbs' and 'unit cells' to 'ncells'.

    Returns
    -------
    info : dict
        The data of the info file
    """
    path = os.path.join(directory, name)
    if not name.endswith("_info"):
        raise ValueError(f"File {path} is not a ALF info file!")

    # Initialize key map for renaming
    default_key_map = {
        "number of orbitals": "norbs",
        "unit cells": "ncells",
    }
    if key_map is not None:
        default_key_map.update(key_map)

    # Read info file contents
    with open(path, "r") as fh:
        data = fh.read()

    # Parse info data to a dict with numerical values (if possible)
    info = dict()
    for line in data.splitlines():
        line = line.strip()
        if not line.startswith("="):
            try:
                if ": " in line:
                    name, valuestr = line.split(": ")
                    values = [s.strip() for s in valuestr.split()]
                    key = name.strip().lower()
                    value = strings_to_numbers(values)
                else:
                    key = line.strip().lower()
                    value = ""

                # Rename key and set value
                if key in default_key_map:
                    key = default_key_map[key]
                info[key] = value

            except ValueError as e:
                raise ParseError(f"Couldn't parse info-line '{line}': {e}")

    return info


def rebin(x, nrebin):
    """Combine each N_rebin bins into one bin.

    If the number of bins (=N0) is not an integer multiple of N_rebin,
    the last N0 modulo N_rebin bins get dismissed.
    """
    if nrebin == 1:
        return x
    n0 = len(x)
    n = n0 // nrebin
    shape = (n,) + x.shape[1:]
    y = np.empty(shape, dtype=x.dtype)
    for i in range(n):
        y[i] = np.mean(x[i * nrebin : (i + 1) * nrebin], axis=0)
    return y


def jack(x, nrebin, nskip):
    """Create jackknife bins out of input bins after skipping and rebinning.

    Parameters
    ----------
    x : array-like object
        Input bins. Bins run over first index.
    nskip : int, default=par.N_skip()
        Number of bins to skip.
    nrebin : int, default=par.N_rebin()
        Number of bins to recombine into one.

    Returns
    -------
    jack_bins : np.ndarray
        Jackknife bins after skipping and rebinning.
    """
    if nskip != 0:
        x = x[nskip:]
    x = rebin(x, nrebin)
    n = len(x)
    y = (np.sum(x, axis=0) - x) / (n - 1)
    return y


def subtract_background(values, backs, ncells):
    n = 0  # latt.invlistk[0, 0]
    norbs = values.shape[1]
    ntau = values.shape[3]
    for orb1, orb2 in itertools.product(range(norbs), repeat=2):
        for tau in range(ntau):
            values[:, orb2, orb1, tau, n] -= ncells * backs[:, orb2] * backs[:, orb1]


def read_data_tau(directory, obs_name, nrebin=0, nskip=0, subtract_back=True):
    info = read_info(directory, obs_name + "_info")
    ntau = info["ntau"]
    ncells = info["ncells"]
    norbs = info["norbs"]
    dtau = info["dtau"]

    # Read lines of the output file
    path = os.path.join(directory, obs_name)
    with open(path, "r") as fh:
        data = fh.read()
    lines = data.splitlines()

    # Sanity check of line numbers
    num_bins0 = len(lines) / (1 + norbs + ncells + ncells * ntau * norbs**2)
    num_bins = int(round(num_bins0))
    if abs(num_bins0 - num_bins) > 1e-10:
        raise ParseError(
            f"Error in reading data: File '{path}', line number does not fit!"
            "Did you forget to clear the output-dir before re-running the simulation?"
        )

    # Initialize output arrays for the data, signs and backgrounds
    values = np.zeros((num_bins, norbs, norbs, ntau, ncells), dtype=np.complex128)
    signs = np.zeros(num_bins, dtype=np.int8)
    backs = np.zeros((num_bins, norbs), dtype=np.complex128)

    # Parse output file
    i = 0
    for ibin in range(num_bins):
        # Parse header line
        header = lines[i].split()
        signs[ibin] = float(header[0])
        # Check bin parameters
        if int(header[1]) != norbs:
            raise BinHeaderError("norbs", obs_name, ibin)
        if int(header[2]) != ncells:
            raise BinHeaderError("ncells", obs_name, ibin)
        if int(header[3]) != ntau:
            raise BinHeaderError("ntau", obs_name, ibin)
        if float(header[4]) != dtau:
            raise BinHeaderError("dtau", obs_name, ibin)
        i += 1

        # First `norbs` lines for background
        for iorb in range(norbs):
            backs[ibin, iorb] = csv_to_complex(lines[i])
            i += 1

        # Parse main data
        for icell in range(ncells):
            # What is this line??? (without braces)
            # line = lines[i]
            # val1, val2 = [_string_to_number(s) for s in line.split()]
            # print(val1, val2)
            i += 1
            for itau in range(ntau):
                for orb1, orb2 in itertools.product(range(norbs), repeat=2):
                    line = lines[i]
                    try:
                        values[ibin, orb1, orb2, itau, icell] = csv_to_complex(line)
                    except ValueError:
                        raise ComplexParseError(lines[i])
                    i += 1

    if nrebin and nskip:
        values = jack(values, nrebin, nskip)
        backs = jack(backs, nrebin, nskip)
        signs = jack(signs, nrebin, nskip)

    if subtract_back:
        subtract_background(values, backs, ncells)

    tau = np.arange(info["ntau"]) * info["dtau"]
    return info, tau, values, signs


def read_green_tau(directory, total=True):
    info, tau, gf_tau, signs = read_data_tau(directory, "Green_tau")
    # Remove orbs
    gf_tau = gf_tau[:, 0, 0]
    # Normalize data
    gf_tau = -gf_tau / (info["ncells"] * info["norbs"])
    if total:
        gf_tau = np.sum(gf_tau, axis=-1)
    return info, tau, gf_tau.real


def read_mean_tau(directory, obs_name):
    folder = os.path.join(directory, obs_name)
    file = os.path.join(folder, "g_dat")
    with open(file, "r") as fh:
        data = fh.read()

    lines = data.splitlines(keepends=False)

    header = lines.pop(0)
    strings = header.split()
    ntau = int(strings[0])
    # nbins = int(strings[1])
    # beta = float(strings[2])
    # norbs = int(strings[3])

    # Read tau, <mean[Tr(X(tau))]> and  <error[Tr(X(tau))]>
    tau = np.zeros(ntau, dtype=np.float64)
    mean = np.zeros(ntau, dtype=np.float64)
    err = np.zeros(ntau, dtype=np.float64)
    for i in range(ntau):
        strings = lines.pop(0).split()
        tau[i] = float(strings[0])
        mean[i] = float(strings[1])
        err[i] = float(strings[2])

    # Read covariance? matrix
    values = np.zeros((ntau * ntau), dtype=np.float64)
    for i, line in enumerate(lines):
        values[i] = float(line)
    mat = values.reshape((ntau, ntau))

    return tau, mean, err, mat


def read_greens_kspace(root):
    path = os.path.join(root, "Green")
    data = np.loadtxt(path, usecols=(1, 2, 3))
    return data.T
