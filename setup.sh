#!/bin/bash

#================================================================================================
# Script for setting up standalone environment for accessing the python libraries and scripts.
#
# Usage:
# cd HLTausAnalysis
# source setup.sh
#
#================================================================================================

echo $HLTAUSANALYSIS_BASE
if [ "x$HLTAUSANALYSIS_BASE" != "x" ]; then
    echo "Standalone environment already loaded"
    return
fi

LOCATION="" 
if [ -n "$CMSSW_BASE" ] ; then
    LOCATION="CMSSW"
    echo "=== LOCATION is $LOCATION. Known problems with Python when 'cmsenv' is set (Perhaps PYTHONPATH?)"
    read -p "Press any key to continue: "
    echo "Continuing ..."
fi

# Detect LXPLUS or MAC OS X (Darwin)
if [ -z "$LOCATION" ] ; then
    ho=`hostname`
    echo $ho
    if [[ `hostname` =~ "lxplus"* ]] ; then
        LOCATION="lxplus"
    elif [[ `hostname` =~ "super"* ]] || [[ $ho =~ "Super"* ]] ; then
        LOCATION="mac"
    elif [[ `hostname` =~ *".cern.ch" ]] ; then
        LOCATION="lxbatch"
    fi
fi
echo "Location is $LOCATION"

# Detect lxplus and jade
if [ "x$LOCATION" = "x" ]; then
    case "$HOSTNAME" in
        lxplus*) LOCATION="lxplus";;
        jade*) LOCATION="jade";;
    esac
fi

if [ "x$LOCATION" = "xlxplus" ]; then
    echo "Sourcing lxplus environments for gcc 4.8 and ROOT 6.02"
    source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/setup.sh
    pushd /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.00/x86_64-slc6-gcc48-opt/root >/dev/null
    source bin/thisroot.sh
    popd >/dev/null
fi

LD_LIBRARY_PATH_APPEND=""
if [[ "$LOCATION" = "lxplus" ]] || [[ "$LOCATION" = "lxbatch" ]]; then
    echo "=== Sourcing $LOCATION environments (Hand-picked from CMSSW_7_6_5)"
    echo "To update:"
    echo "1) create a developer area (cmsrel):"
    echo "cmsrel CMSWW_X_Y_Z"
    echo "2) source the CMSSW environment (cmsenv):"
    echo "cd  CMSWW_X_Y_Z/src/"
    echo "cmcsev"
    echo "3) Look the new paths with 'scram tool list' and 'scram tool info'"
    echo "[See setup.csh for more details]"

    # scram tool info gcc-cxxcompiler (Look for GCC_CXXCOMPILER_BASE)
    GCC_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3

    # scram tool info root_interface (Look for ROOT_INTERFACE_BASE)
    ROOTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/lcg/root/6.02.12-kpegke4

    # scram tool info xrootd (Look for XROOTD_BASE)
    XROOTD_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/xrootd/4.0.4-kpegke2

    # scram tool info xz (Look for XZ_BASE)
    XZ_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/xz/5.2.1

    # scram tool info python
    PYTHON_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/python/2.7.6-kpegke

    # Set the run-time shared library loader (ld.so) an extra set of directories to look for when searching for shared libraries.
    LD_LIBRARY_PATH_APPEND=$ROOTSYS/lib:$GCC_BASE/lib64:$GCC_BASE/lib:$XROOTD_BASE/lib:$XZ_BASE/lib:$PYTHON_BASE/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/libjpg/8b-cms/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/libpng/1.6.16/lib
    
    # Tell the shell which directories to search for executable files
    PATH=$ROOTSYS/bin:$GCC_BASE/bin:$XROOTD_BASE/bin:$PATH
    export PATH

    if [[ -n "$PYTHONPATH" ]]; then
	PYTHONPATH="$ROOTSYS/lib:$PYTHONPATH"; export PYTHONPATH 
    else
	PYTHONPATH="$ROOTSYS/lib" ; export PYTHONPATH 
    fi
    LD_LIBRARY_PATH_APPEND=$GCC_BASE/lib64:$GCC_BASE/lib:$XROOTD_BASE/lib:$XZ_BASE/lib:$PYTHON_BASE/lib

    export PATH=$ROOTSYS/bin:$GCC_BASE/bin:$XROOTD_BASE/bin:$PATH 

    pushd $ROOTSYS >/dev/null
    source bin/thisroot.sh
    popd >/dev/null
fi

if [[ "$LOCATION" = "mac" ]]; then
    echo "=== Setting ROOT version 5-34-00-patches"
    source /sw/Codes/root/v6-06-08/bld-root/bin/thisroot.sh
    if [[ -n "$PYTHONPATH" ]]; then 
	PYTHONPATH="$ROOTSYS/lib:$PYTHONPATH" 
        export PYTHONPATH 
    else
        PYTHONPATH="$ROOTSYS/lib"
        export PYTHONPATH 
    fi
fi

echo "=== Appending to LD_LIBRARY_PATH"
export HLTAUSANALYSIS_BASE=$PWD
LD_LIBRARY_PATH_APPEND=$HLTAUSANALYSIS_BASE/NtupleAnalysis/lib:$LD_LIBRARY_PATH_APPEND

if [ "x$LD_LIBRARY_PATH" = "x" ]; then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_APPEND
else
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_APPEND:$LD_LIBRARY_PATH
fi

echo "=== Creating symbolic links and hidden directories for $LOCATION"
export PPATHPREFIX=.python
if [ "x$LOCATION" = "xCMSSW" ]; then
    if [ ! -e $CMSSW_BASE/python/HLTausAnalysis/NtupleAnalysis ]; then
        ln -s $HLTAUSANALYSIS_BASE/NtupleAnalysis/python $CMSSW_BASE/python/HLTausAnalysis/NtupleAnalysis
    fi
    export PPATHPREFIX="$CMSSW_BASE/python"
fi

# Need to create the following also on lxplus for limit calculation
if [ ! -e $PPATHPREFIX/HLTausAnalysis ]; then
    mkdir -p $PPATHPREFIX/HLTausAnalysis
    touch $PPATHPREFIX/HLTausAnalysis/__init__.py
fi
for DIR in NtupleAnalysis; do
    if [ ! -e $PPATHPREFIX/HLTausAnalysis/$DIR ]; then
        ln -s $HLTAUSANALYSIS_BASE/$DIR/python $PPATHPREFIX/HLTausAnalysis/$DIR
        touch $PPATHPREFIX/HLTausAnalysis/$DIR/__init__.py
        for d in $PPATHPREFIX/HLTausAnalysis/$DIR/*; do
            if [ -d $d ]; then
                touch $d/__init__.py
            fi
        done
    fi
done
for DIR in `ls NtupleAnalysis/src` ; do
    if [[ ! -e $PPATHPREFIX/HLTausAnalysis/$DIR ]] && [[ -e $HLTAUSANALYSIS_BASE/NtupleAnalysis/src/$DIR/python ]]; then
        ln -s $HLTAUSANALYSIS_BASE/NtupleAnalysis/src/$DIR/python $PPATHPREFIX/HLTausAnalysis/$DIR
        touch $PPATHPREFIX/HLTausAnalysis/$DIR/__init__.py
        for d in $PPATHPREFIX/HLTausAnalysis/$DIR/*; do
            if [ -d $d ]; then
                touch $d/__init__.py
            fi
        done
    fi
done

if [ "x$PYTHONPATH" = "x" ]; then
    export PYTHONPATH=$PWD/.python
else
    export PYTHONPATH=$PWD/.python:$PYTHONPATH
fi

#echo "=== Setting PATH variable"
PATH="${HLTAUSANALYSIS_BASE}/NtupleAnalysis/scripts:${PATH}" ; export PATH

echo "=== The environment variables set are:"
echo "LOCATION is $LOCATION"
echo "HLTAUSANALYSIS_BASE is $HLTAUSANALYSIS_BASE"
echo "PATHPREFIX is $PATHPREFIX"
echo "ROOTSYS is $ROOTSYS"
echo "LD_LIBRARY_PATH is $LD_LIBRARY_PATH"
echo "PYTHONPATH is $PYTHONPATH"
echo "PATH is $PATH"
