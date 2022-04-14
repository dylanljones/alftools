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
            f90nml.patch(self._path, self._nml, self._path + '_temp')
            shutil.move(self._path + '_temp', self._path)

    def __str__(self):
        return str(self._nml)
