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
from .parameters import Parameters
from .utils import (
    conf,
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
        cmd = os.path.join(conf["ALF_DIR"], "Analysis", "ana_hdf5.out")
        if files != "*":
            cmd += " " + files
    else:
        cmd = os.path.join(conf["ALF_DIR"], "Analysis", "ana.out") + " " + files
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


def jacknife_bins(x, nrebin=1, nskip=0):
    """Create jackknife bins out of input bins after skipping and rebinning.

    Parameters
    ----------
    x : (Nbin, ...) np.ndarray
        Input bins. Bins run over first index.
    nskip : int
        Number of bins to skip.
    nrebin : int
        Number of bins to combine into one bin.

    Returns
    -------
    jack_bins : (M, ...) np.ndarray
        Jackknife bins after skipping and rebinning.
    """
    if nskip != 0:
        # Skip first nskip bins
        x = x[nskip:]
    if nrebin > 1:
        # Rebin each nrebin into one bin
        nbins = len(x) // nrebin
        bins = np.empty((nbins,) + x.shape[1:], dtype=x.dtype)
        for i in range(nbins):
            bins[i] = np.mean(x[i * nrebin : (i + 1) * nrebin], axis=0)  # noqa
    else:
        bins = x
    return (np.sum(bins, axis=0) - bins) / (len(bins) - 1)


def error(jacks):
    """Calculates the errors of the given jackknife bins.

    Parameters
    ----------
    jacks : (Nbin, ....) np.ndarray
        The input Jackknife bins.

    Returns
    -------
    errors : (...) np.ndarray
        The error values.
    """
    # n = len(jacks)
    return np.sqrt(np.var(jacks, axis=0) * len(jacks))


def mean(jacks):
    """Calculates the mean of the given jackknife bins.

    Parameters
    ----------
    jacks : (Nbin, ....) np.ndarray
        The input Jackknife bins.

    Returns
    -------
    mean : (...) np.ndarray
        The mean of the input bins.
    """
    return np.mean(jacks.real, axis=0)


def subtract_background(values, backs, ncells):
    n = 0  # latt.invlistk[0, 0]
    norbs = values.shape[1]
    ntau = values.shape[3]
    for orb1, orb2 in itertools.product(range(norbs), repeat=2):
        for tau in range(ntau):
            values[:, orb2, orb1, tau, n] -= ncells * backs[:, orb2] * backs[:, orb1]


def subtract_background_mat(values, backs, ncells):
    n = 0  # latt.invlistk[0, 0]
    norbs = values.shape[1]
    ntau = values.shape[3]
    for orb1, orb2 in itertools.product(range(norbs), repeat=2):
        for tau in range(ntau):
            values[:, orb2, orb1, tau, n, n] -= ncells * backs[:, orb2] * backs[:, orb1]


def read_data_latt(directory, obs_name, nrebin=None, nskip=None, subtract_back=True):
    """Read the data file of `Obser_Latt` observable types from ALF.

    The `Obser_Latt` observables include equal-time and time-displaced observables,
    for example the Green's function :math:`G(0)` and :math:`G(\tau, 0)`.

    Parameters
    ----------
    directory : str
        The simulation root directory.
    obs_name : str
        The name of the observable
    nrebin : int, optional
        The number of rebinning to use. By default, the value from the parameter file
        is used.
    nskip : int, optional
        The number of bins to skip. By default, the value from the parameter file
        is used.
    subtract_back : bool, optional
        If True, the background will be subtracted from the data.

    Returns
    -------
    tau : (L,) np.ndarray
        The imaginary time values. If the observable is an equal time observable, L=1.
    values : (Nbin, L, Norb, Norb, N) np.ndarray
        The data of the output file. If the observable is an equal time observable, L=1.
    signs : (Nbin, ) np.ndarray
        The signs in each bin.
    info : dict
        The data of the observable info file.
    """
    params = Parameters(directory)
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
        if len(header) > 3:
            # obs_tau
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
            # K-Point (without braces)
            # kx, ky = [_string_to_number(s) for s in lines[i].split()]
            i += 1
            for itau in range(ntau):
                for orb1, orb2 in itertools.product(range(norbs), repeat=2):
                    line = lines[i]
                    try:
                        values[ibin, orb1, orb2, itau, icell] = csv_to_complex(line)
                    except ValueError:
                        raise ComplexParseError(lines[i])
                    i += 1

    # Get default values from parameter file
    if nrebin is None:
        nrebin = params.get_errors("n_rebin")
    if nskip is None:
        nskip = params.get_errors("n_skip")
    # Process jacknife bins
    if nrebin and nskip:
        values = jacknife_bins(values, nrebin, nskip)
        backs = jacknife_bins(backs, nrebin, nskip)
        signs = jacknife_bins(signs, nrebin, nskip)
    if subtract_back:
        subtract_background(values, backs, ncells)

    tau = np.arange(info["ntau"]) * info["dtau"]
    return tau, values, signs, info


def read_data_mat(directory, obs_name, nrebin=None, nskip=None, subtract_back=True):
    """Read the data file of `Obser_mat` observable types from ALF.

    The `Obser_mat` observables include equal-time and time-displaced observables in
    matrix form, for example the Green's function :math:`G(0)` and :math:`G(\tau, 0)`.

    Parameters
    ----------
    directory : str
        The simulation root directory.
    obs_name : str
        The name of the observable
    nrebin : int, optional
        The number of rebinning to use. By default, the value from the parameter file
        is used.
    nskip : int, optional
        The number of bins to skip. By default, the value from the parameter file
        is used.
    subtract_back : bool, optional
        If True, the background will be subtracted from the data.

    Returns
    -------
    tau : (L,) np.ndarray
        The imaginary time values. If the observable is an equal time observable, L=1.
    values : (Nbin, L, Norb, Norb, N, N) float np.ndarray
        The data of the output file. If the observable is an equal time observable, L=1.
    signs : (Nbin, ) np.ndarray
        The signs in each bin.
    info : dict
        The data of the observable info file.
    """
    params = Parameters(directory)
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
    num_elements = ncells**2 * norbs**2 * (1 + ntau)
    num_bins0 = len(lines) / (1 + norbs + num_elements)
    num_bins = int(round(num_bins0))
    if abs(num_bins0 - num_bins) > 1e-10:
        raise ParseError(
            f"Error in reading data: File '{path}', line number does not fit!"
            "Did you forget to clear the output-dir before re-running the simulation?"
        )

    # Initialize output arrays for the data, signs and backgrounds
    values = np.zeros(
        (num_bins, norbs, norbs, ntau, ncells, ncells), dtype=np.complex128
    )
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
        if len(header) > 3:
            # obs_tau
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
            for jcell in range(ncells):
                line = lines[i]
                i1, i2 = [int(s) - 1 for s in line.split()]
                i += 1
                for itau in range(ntau):
                    for orb1, orb2 in itertools.product(range(norbs), repeat=2):
                        line = lines[i]
                        try:
                            val = csv_to_complex(line)
                            values[ibin, orb1, orb2, itau, i1, i2] = val
                        except ValueError:
                            raise ComplexParseError(lines[i])
                        i += 1

    # Get default values from parameter file
    if nrebin is None:
        nrebin = params.get_errors("n_rebin")
    if nskip is None:
        nskip = params.get_errors("n_skip")
    # Process jacknife bins
    if nrebin and nskip:
        values = jacknife_bins(values, nrebin, nskip)
        backs = jacknife_bins(backs, nrebin, nskip)
        signs = jacknife_bins(signs, nrebin, nskip)

    if subtract_back:
        subtract_background_mat(values, backs, ncells)

    tau = np.arange(info["ntau"]) * info["dtau"]
    return tau, values, signs, info


# == Green's functions =================================================================


def read_green_eq(directory, iorb=0, nrebin=None, nskip=None):
    """Read the equal time Green's function data.

    Parameters
    ----------
    directory : str
        The simulation root directory.
    iorb : int, optional
        The orbital index. If None the Green's function for all orbitals are returned.
    nrebin : int, optional
        The number of rebinning to use. By default, the value from the parameter file
        is used.
    nskip : int, optional
        The number of bins to skip. By default, the value from the parameter file
        is used.

    Returns
    -------
    gf_eq : np.ndarray
        The equal time Green's function. If an orbital index is given the shape of the
        resulting array is (Nbin, N), otherwise it is (Nbin, Norb, Norb, N).
    """
    tau, gf_eq, signs, info = read_data_latt(directory, "Green_eq", nrebin, nskip)
    # Remove orbs and tau axis
    gf_eq = gf_eq[:, :, :, 0]
    if iorb is not None:
        gf_eq = gf_eq[:, iorb, iorb]
    # Normalize data
    gf_eq = -gf_eq / (info["ncells"] * info["norbs"])
    errs = error(gf_eq)
    return gf_eq.real, errs, info


def read_green_tau(directory, iorb=0, total=False, nrebin=None, nskip=None):
    """Read the time displaced Green's function data.

    Parameters
    ----------
    directory : str
        The simulation root directory.
    iorb : int, optional
        The orbital index. If None the Green's function for all orbitals are returned.
    total : bool, optional
        If True, the Green's function is summed over all lattice sites.
    nrebin : int, optional
        The number of rebinning to use. By default, the value from the parameter file
        is used.
    nskip : int, optional
        The number of bins to skip. By default, the value from the parameter file
        is used.

    Returns
    -------
    gf_tau : np.ndarray
        The time displaced Green's function.
    """
    tau, gf_tau, signs, info = read_data_latt(directory, "Green_tau", nrebin, nskip)
    # Remove orbs
    if iorb is not None:
        gf_tau = gf_tau[:, iorb, iorb]
    # Normalize data
    gf_tau = -gf_tau / (info["ncells"] * info["norbs"])
    if total:
        gf_tau = np.sum(gf_tau, axis=-1)

    errs = error(gf_tau)
    return tau, gf_tau.real, errs, info


def read_greenmat_eq(directory, iorb=0, nrebin=None, nskip=None):
    """Read the equal time Green's function data in matrix representation.

    Parameters
    ----------
    directory : str
        The simulation root directory.
    iorb : int, optional
        The orbital index. If None the Green's function for all orbitals are returned.
    nrebin : int, optional
        The number of rebinning to use. By default, the value from the parameter file
        is used.
    nskip : int, optional
        The number of bins to skip. By default, the value from the parameter file
        is used.

    Returns
    -------
    gf_eq : np.ndarray
        The equal time Green's function. If an orbital index is given the shape of the
        resulting array is (Nbin, N, N), otherwise it is (Nbin, Norb, Norb, N, N).
    """
    tau, gf_eq, signs, info = read_data_mat(directory, "Greenmat_eq", nrebin, nskip)
    # Remove orbs and tau axis
    gf_eq = gf_eq[:, :, :, 0]
    if iorb is not None:
        gf_eq = gf_eq[:, iorb, iorb]
    # Normalize data
    gf_eq = -gf_eq / (info["ncells"] * info["norbs"])

    errs = error(gf_eq)
    return gf_eq.real, errs, info


def read_greenmat_tau(directory, iorb=0, total=False, nrebin=None, nskip=None):
    """Read the time displaced Green's function data in matrix representation.

    Parameters
    ----------
    directory : str
        The simulation root directory.
    iorb : int, optional
        The orbital index. If None the Green's function for all orbitals are returned.
    total : bool, optional
        If True, the Green's function is summed over all lattice sites.
    nrebin : int, optional
        The number of rebinning to use. By default, the value from the parameter file
        is used.
    nskip : int, optional
        The number of bins to skip. By default, the value from the parameter file
        is used.

    Returns
    -------
    gf_tau : np.ndarray
        The time displaced Green's function.
    """
    tau, gf_tau, signs, info = read_data_mat(directory, "Greenmat_tau", nrebin, nskip)
    # Remove orbs
    if iorb is not None:
        gf_tau = gf_tau[:, iorb, iorb]
    # Normalize data
    gf_tau = -gf_tau / (info["ncells"] * info["norbs"])
    if total:
        gf_tau = np.trace(gf_tau, axis1=-2, axis2=-1)
    errs = error(gf_tau)
    return tau, gf_tau.real, errs, info


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
    mn = np.zeros(ntau, dtype=np.float64)
    err = np.zeros(ntau, dtype=np.float64)
    for i in range(ntau):
        strings = lines.pop(0).split()
        tau[i] = float(strings[0])
        mn[i] = float(strings[1])
        err[i] = float(strings[2])

    # Read covariance? matrix
    values = np.zeros((ntau * ntau), dtype=np.float64)
    for i, line in enumerate(lines):
        values[i] = float(line)
    mat = values.reshape((ntau, ntau))

    return tau, mn, err, mat


def read_greens_kspace(root):
    path = os.path.join(root, "Green")
    data = np.loadtxt(path, usecols=(1, 2, 3))
    return data.T
