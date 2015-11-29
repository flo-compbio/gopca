#!/bin/bash

sphinx-apidoc -e -o source/api ../gopca
rm source/api/modules.rst
rm source/api/gopca.rst
rm source/api/gopca.scripts.rst
rm source/api/gopca.plotting.rst
