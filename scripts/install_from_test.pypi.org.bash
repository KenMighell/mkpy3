#!/bin/bash
url='https://test.pypi.org/simple/'
pkg='mkpy3'
pip uninstall $pkg
pip install --index-url $url --extra-index-url $url $pkg
#EOF
