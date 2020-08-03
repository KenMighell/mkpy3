#!/usr/bin/env python3
import os
import subprocess
import glob

verbose = False

cwd = os.getcwd()
print(cwd, ' =$PWD')
for name in glob.glob('./mkpy3*.py'):
    result = subprocess.run(['python3', name], capture_output=True)
    returncode = result.returncode
    if (returncode != 0):
        print(name, ': ***** ERROR FOUND *****')
    else:
        if (verbose): print(name, ': OK')
    pass#if
pass#for

#EOF

