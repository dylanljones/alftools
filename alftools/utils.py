# coding: utf-8
#
# This code is part of alftools.
#
# Copyright (c) 2022, Dylan Jones

import os
import logging
import numpy as np

ALF_DIR = os.environ["ALF_DIR"]


logger = logging.getLogger(__name__.split(".")[0])

# Logging format
frmt = "[%(asctime)s] %(name)s:%(levelname)-8s - %(message)s"
formatter = logging.Formatter(frmt, datefmt="%H:%M:%S")

# Set up console logger
sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)
sh.setFormatter(formatter)
logger.addHandler(sh)

# Set logging level
logger.setLevel(logging.INFO)
logging.root.setLevel(logging.NOTSET)


class ParseError(ValueError):
    pass


class ComplexParseError(ParseError):
    def __init__(self, string):
        super().__init__(f"complex() arg is a malformed string: {string}")


def call(cmd, cwd=None, verbose=False):
    old_cwd = ""
    if cwd:
        old_cwd = os.getcwd()
        os.chdir(cwd)
    if not verbose:
        cmd += " > /dev/null"
    os.system(cmd)
    if cwd:
        os.chdir(old_cwd)


def string_to_number(string):
    """Converts a string to an integer, float or a complex if possible.

    Parameters
    ----------
    string : str
        The input string to convert.

    Returns
    -------
    value : int or float or complex or str
        The converted input string `string` if the conversion to either type was
        possible. Otherwise, the original input string `string` is returned.

    Examples
    --------
    >>> string_to_number("1")
    1

    >>> string_to_number("1.2")
    1.2

    >>> string_to_number("1.1 + 2.2j")
    (1.1+2.2j)

    >>> string_to_number("test")
    'test'
    """
    # Try to parse as int
    try:
        return int(string)
    except ValueError:
        pass
    # Try to parse as float
    try:
        return float(string)
    except ValueError:
        pass
    # Try to parse as complex
    try:
        return complex(string.replace(" ", ""))
    except ValueError:
        pass
    return string


def strings_to_numbers(strings):
    """Converts strings to an integer, float, complex or array of values if possible.

    If a single string is passed it is converted to a scalar value. Otherwise, the
    strings are converted to a numpy array if possible.

    Parameters
    ----------
    strings : Iterable of str
        The input strings to convert.

    Returns
    -------
    values : int or float or complex or np.ndarray or Iterable of str
        The converted input strings `strings` if the conversion to either type was
        possible. Otherwise, the original input strings `s` are returned.

    Examples
    --------
    >>> strings_to_numbers(["1"])
    1

    >>> strings_to_numbers(["1.2"])
    1.2

    >>> strings_to_numbers(["1.1 + 2.2j"])
    (1.1+2.2j)

    >>> strings_to_numbers(["1", "2"])
    array([1, 2])

    >>> strings_to_numbers(["1.1", "2.2"])
    array([1.1, 2.2])
    """
    if len(strings) == 1:
        return string_to_number(strings[0])
    else:
        values = [string_to_number(s) for s in strings]
        if not any(isinstance(x, str) for x in values):
            return np.array(values)
        return values


def csv_to_complex(s):
    s_frmt = s.replace(" ", "").replace(",", "+").replace("+-", "-").replace(")", "j)")
    try:
        return complex(s_frmt)
    except ValueError:
        raise ComplexParseError(s_frmt)
