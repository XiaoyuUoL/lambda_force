#!/bin/bash

# This file will start automatically in your docker run. You can assume the presence of a
# molecule.yml and a calculator.yml in the work directory.

# Make sure that after this script finishes a result.yml exists.
# The workdir_bundle.tar.gz will also be staged out for debugging purposes, if you create it.
cd /opt/binary_folder
python main.py

# A good way to pack all files smaller than e.g 500k for stageout is:
#cd /opt
#find result_folder -type f -size -500k -print0 | xargs -0 tar czf workdir_bundle.tar.gz