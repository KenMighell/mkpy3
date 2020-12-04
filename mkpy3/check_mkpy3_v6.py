#!/usr/bin/env python3

# file://check_mkpy3.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def check_mkpy3():
    """
Purpose: Unit tests for mkpy3  (https://github.com/KenMighell/mkpy3)

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    import os
    import glob
    import subprocess
    #
    # ASCII color codes
    CEND = '\33[0m'
    # CRED = '\33[31m'  # red
    CRED = '\33[1;31m'  # bold red
    # CGREEN = '\33[32m'  # green
    CGREEN = '\33[1;32m'  # bold green
    #
    cwd = os.getcwd()
    print(cwd, ' =$PWD')
    good = 0
    bad = 0
    filel = sorted(glob.glob('./mkpy3_*.py'))
    filel += sorted(glob.glob('./xmkpy3_*.py'))
    sz = len(filel)
    print()
    print(sz, 'files to check')
    print()
    for k, name in enumerate(filel, start=1):
        print('[%d] %s : ' % (k, name.strip()), end='')
        result = subprocess.run(['python3', name], capture_output=True)
        returncode = result.returncode
        if (returncode != 0):
            bad += 1
            print(CRED + '***** ERROR FOUND *****' + CEND)
            print()
        else:
            good += 1
            print(CGREEN + 'OK' + CEND)
        # pass:if
    # pass:for
    print()
    if (good == sz):
        msg_ = '\nAll %d files PASS  :-)' % good
        print(CGREEN + msg_ + CEND)
    else:
        msg_ = '%d of %d files FAIL  8=X' % (bad, sz)
        print(CRED + msg_ + CEND)
    # pass:if
    print()
# pass:def


# =============================================================================


if (__name__ == '__main__'):
    check_mkpy3()
# pass:if

# EOF
