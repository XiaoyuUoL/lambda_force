#!/bin/bash

export SOFTWARE=/opt/software_folder

# orca
export ORCA=$SOFTWARE/orca_5_0_4_linux_x86-64_shared_openmpi411
export PATH=$ORCA:$PATH
export LD_LIBRARY_PATH=$ORCA:$LD_LIBRARY_PATH

$ORCA/orca S0freq.inp > S0freq.log