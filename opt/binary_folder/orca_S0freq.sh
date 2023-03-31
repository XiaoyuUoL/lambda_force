#!/bin/bash

export SOFTWARE=/opt/software_folder

# openmpi
export OPENMPI_DIR=$SOFTWARE/openmpi-4.1.5
export PATH=$OPENMPI_DIR/bin:$PATH
export LD_LIBRARY_PATH=$OPENMPI_DIR/lib:$LD_LIBRARY_PATH

# orca
export ORCA=$SOFTWARE/orca_5_0_4_linux_x86-64_shared_openmpi411
export PATH=$ORCA:$PATH
export LD_LIBRARY_PATH=$ORCA:$LD_LIBRARY_PATH

$ORCA/orca S0freq.inp > S0freq.log