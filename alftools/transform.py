# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import numpy as np
import gftool as gt


def f_tau2iw(f_tau, beta):
    """Transforms a function from imaginary times `τ` to Matsubara frequencies `iω_n`.

    Parameters
    ----------
    f_tau : (..., N) np.ndarray
        The input function .math:`f(τ)` at imaginary times .math:`τ ∈ [0, β]`.
    beta : float
        The inverse temperature .math:`β = 1 / k_B T`.

    Returns
    -------
    iw : ((N + 1) / 2, ) np.ndarray
        The non-negative fermionic Matsubara frequencies .math:`iω_n`.
    f_iw : (..., (N + 1) / 2) np.ndarray
        The Fourier transform of `f_tau` for non-negative fermionic Matsubara
        frequencies .math:`iω_n`.
    """
    num_times = f_tau.shape[-1]
    n_points = range(num_times // 2)
    iws = gt.matsubara_frequencies(n_points, beta=beta)
    f_iw = gt.fourier.tau2iw(f_tau, beta=beta)
    return iws, f_iw


def f_iw2z(f_iw, iw, num=1000, eta=1e-6):
    """Transforms a function from Matsubara frequencies `iω_n` to real frequencies `z`.

    Uses the Pade analytical continuation algorithm.

    Parameters
    ----------
    f_iw : (..., N) np.ndarray
        The Fourier transform of `f_tau` for non-negative fermionic Matsubara
        frequencies .math:`iω_n`.
    iw : (..., N) np.ndarray
        The non-negative fermionic Matsubara frequencies .math:`iω_n`.
    num : int, optional
        The number of points M at which the continuated function is evalued.
        The default is 1000.
    eta : float
        The complex broadening used for the real frequencies .math:`z = ω + i η`.
        The default is `1e-6`.

    Returns
    -------
    z : (M, ) np.ndarray
        The real frequencies with a complex broadening term  .math:`z = ω + i η`.
    f_z : (..., M) np.ndarray
        The continuated input function evaluated at the frequenzies `z`.
    """
    pade = gt.polepade.continuation(iw, f_iw, degree=-1, moments=[1])
    print(f"[{pade.zeros.size}/{pade.poles.size}]")
    z = np.linspace(-5, 5, num=num) + 1j * eta
    f_z = pade.eval_polefct(z)
    return z, f_z


def f_tau2z(f_tau, beta, num=1000, eta=1e-6):
    iw, f_iw = f_tau2iw(f_tau, beta)
    return f_iw2z(f_iw, iw, num, eta)
