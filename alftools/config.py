# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
from configparser import ConfigParser

ALF_DIR = os.environ["ALF_DIR"]


def _read_config(path=""):
    if not path:
        path = os.path.join(os.getcwd(), ".alftools.cfg")
    parser = ConfigParser()
    parser.read(path)
    return parser


def init_config():
    parsed = _read_config()
    paths = parsed["paths"]
    config = dict()
    config["alf_dir"] = paths.get("alf_dir", fallback="") or os.environ["ALF_DIR"]
    config["out"] = paths.get("out_dir", fallback="")
    return config
