# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import numpy as np


def read_error_tau(directory, name):
    folder = os.path.join(directory, name)
    file = os.path.join(folder, list(os.listdir(folder))[0])
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
