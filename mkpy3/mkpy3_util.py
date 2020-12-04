#!/usr/bin/env python3

# file://mkpy3_util.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute


def mkpy3_util_str2bool(v):
    """
Utility function for argparse.
    """
    import argparse
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    # pass:if
# pass:def


def mkpy3_util_accept_str_or_int(v):
    """
Utility function for argparse.
    """
    import argparse
    if isinstance(v, int):
        return str(v)
    elif isinstance(v, str):
        return v
    else:
        raise argparse.ArgumentTypeError('str or int value expected.')
    # pass:if
# pass:def


def mkpy3_util_check_file_exists(filename, overwrite):
    """
Utility function.
    """
    import os
    import sys
    assert(isinstance(filename, str))
    assert(isinstance(overwrite, bool))
    msg = 'Requested output file already exists (overwrite=False):\n'
    if (not overwrite):
        if (os.path.isfile(filename)):
            print('\n***** ERROR *****\n\n%s' % (msg))
            print("new_filename='%s'\n" % filename)
            sys.exit(1)
        # pass:if
    # pass:if
# pass:def


if (__name__ == '__main__'):
    pass
# pass#if

# EOF
