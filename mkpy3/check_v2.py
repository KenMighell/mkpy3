#!/usr/bin/env python3

# file:///check_v2.py

version_ = '2020AUG13T1033 0.02'

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute

if (__name__ == '__main__'): 
    """
Purpose: Unit tests for mkpy3  (https://github.com/KenMighell/mkpy3)

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    import os
    import glob
    import subprocess
    #
    cwd = os.getcwd()
    print(cwd, ' =$PWD')
    good = 0
    bad = 0
    print()
    filel = sorted(glob.glob('./mkpy3*.py'))
    for name in filel:
        result = subprocess.run(['python3', name], capture_output=True)
        returncode = result.returncode
        if (returncode != 0):
            bad += 1
            print('\r',name, ': ***** ERROR FOUND *****\n')
        else:
            good += 1
            print('\r',name,': OK',end='')
            print('                                                        \r',
              end='')
        pass#if
    pass#for
    print()
    if (good == len(filel)):
        print('\nAll %d files PASS  :-)' % good)
    else:
        print('%d of %d files FAIL  8=X' % (bad,good))
    pass#if
    print()
pass#if
#EOF
