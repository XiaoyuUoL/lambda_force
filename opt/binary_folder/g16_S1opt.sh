#!/bin/bash

export SOFTWARE=/opt/software_folder

# g16
export G16=$SOFTWARE
export g16root=$G16
export GAUSS_EXEDIR=$G16/g16:$G16/g16/bsd
export GAUSS_SCRDIR=/tmp
export PATH=$G16/g16:$PATH
export LD_LIBRARY_PATH=PATH=$G16/g16:$LD_LIBRARY_PATH

g16 S1opt.gjf
formchk S1opt.chk > /dev/null