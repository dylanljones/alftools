# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import f90nml
import shutil
from collections import abc


class Parameters(abc.Mapping):
    def __init__(self, root_or_path):
        if os.path.isdir(root_or_path):
            root_or_path = os.path.join(root_or_path, "parameters")
        self._path = root_or_path
        self._nml = f90nml.read(self._path)

    def __len__(self):
        return len(self._nml)

    def __iter__(self):
        return iter(self._nml)

    def __getitem__(self, key):
        return dict(self._nml[key])

    def keys(self):
        return self._nml.keys()

    def values(self):
        return self._nml.values()

    def items(self):
        return self._nml.items()

    def get(self, section, key="", default=None):
        try:
            sec = self._nml.get(section)
            if key:
                return sec[key]
            return dict(sec)
        except KeyError:
            return default

    def set(self, section, key, value):
        self._nml[section][key] = value

    def save(self, name=""):
        if name:
            path = os.path.join(os.path.dirname(self._path), name)
            self._path = path
            f90nml.write(self._nml, path)
        else:
            f90nml.patch(self._path, self._nml, self._path + "_temp")
            shutil.move(self._path + "_temp", self._path)

    def __str__(self):
        return str(self._nml)

    # -- Getters -----------------------------------------------------------------------

    def get_lattice(self, key, default=None):
        return self.get("var_lattice", key, default)

    def get_model_generic(self, key, default=None):
        return self.get("var_model_generic", key, default)

    def get_qmc(self, key, default=None):
        return self.get("var_qmc", key, default)

    def get_errors(self, key, default=None):
        return self.get("var_errors", key, default)

    def get_temp(self, key, default=None):
        return self.get("var_temp", key, default)

    def get_hubbard(self, key, default=None):
        return self.get("var_hubbard", key, default)

    def get_hubbard_plain_vanilla(self, key, default=None):
        return self.get("var_hubbard_plain_vanilla", key, default)

    def get_tv(self, key, default=None):
        return self.get("var_tv", key, default)

    def get_kondo(self, key, default=None):
        return self.get("var_kondo", key, default)

    def get_lrc(self, key, default=None):
        return self.get("var_lrc", key, default)

    def get_z2_matter(self, key, default=None):
        return self.get("var_z2_matter", key, default)

    # -- Setters -----------------------------------------------------------------------

    def set_lattice(self, key, value):
        self.set("var_lattice", key, value)

    def set_model_generic(self, key, value):
        self.set("var_model_generic", key, value)

    def set_qmc(self, key, value):
        self.set("var_qmc", key, value)

    def set_errors(self, key, value):
        self.set("var_errors", key, value)

    def set_temp(self, key, value):
        self.set("var_temp", key, value)

    def set_hubbard(self, key, value):
        self.set("var_hubbard", key, value)

    def set_hubbard_plain_vanilla(self, key, value):
        self.set("var_hubbard_plain_vanilla", key, value)

    def set_tv(self, key, value):
        self.set("var_tv", key, value)

    def set_kondo(self, key, value):
        self.set("var_kondo", key, value)

    def set_lrc(self, key, value):
        self.set("var_lrc", key, value)

    def set_z2_matter(self, key, value):
        self.set("var_z2_matter", key, value)
