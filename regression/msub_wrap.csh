#!/bin/tcsh

printenv

setenv host `echo $HOST | sed -e 's/[.].*//g'`
echo " "
switch ( "${host}" )
case tu-*:
    echo "/opt/MOAB/bin/msub -A access $HOME/draco/regression/cdash-release-tu.msub"
    /opt/MOAB/bin/msub -A access $HOME/draco/regression/cdash-release-tu.msub
    echo "/opt/MOAB/bin/msub -A access $HOME/draco/regression/cdash-debug-tu.msub"
    /opt/MOAB/bin/msub -A access $HOME/draco/regression/cdash-debug-tu.msub 
    breaksw
case yr-*:
    echo "/opt/MOAB/bin/msub $HOME/draco/regression/cdash-release-yr.msub"
    /opt/MOAB/bin/msub $HOME/draco/regression/cdash-release-yr.msub 
    echo "/opt/MOAB/bin/msub $HOME/draco/regression/cdash-debug-yr.msub"
    /opt/MOAB/bin/msub $HOME/draco/regression/cdash-debug-yr.msub 
    breaksw
default:
    echo "I don't know how to run regression on host = ${host}."
    breaksw
endsw


