# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import numpy as np
import lattpy as lp


LATTICE_BASIS = {
    "Chain": lp.Lattice.chain,
    "Square": lp.Lattice.square,
    "Cube": lp.Lattice.sc
}


def lattice_from_info(info):
    l1 = info["l1"]
    l2 = info["l2"]
    a1 = info["a1"]
    a2 = info["a2"]
    ncelss = info["ncells"]
    size = l1[0], l2[1]

    vecs = np.array([a1, a2])
    latt = lp.Lattice(vecs)
    latt.add_atom()
    latt.add_connections()
    latt.build(size, primitive=True)
    assert latt.num_cells == ncelss
    return latt


def lattice_from_params(params):
    var_latt = params["var_lattice"]
    latt_type = var_latt["lattice_type"]
    size = var_latt["l1"], var_latt["l2"]

    latt = LATTICE_BASIS[latt_type]()
    latt.add_atom()
    latt.add_connections()
    latt.build(size, primitive=True)

    return latt
