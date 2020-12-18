#!/bin/bash
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#
# example:
#
#$ python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#Uploading distributions to https://test.pypi.org/legacy/
#Enter your username: __token__
#Enter your password: 
#Uploading mkpy3-0.2b1-py3-none-any.whl
#100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 50.3k/50.3k [00:02<00:00, 22.7kB/s]
#Uploading mkpy3-0.2b1.tar.gz
#100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 38.2k/38.2k [00:01<00:00, 27.7kB/s]
#
#View at:
#https://test.pypi.org/project/mkpy3/0.2b1/
#
#EOF