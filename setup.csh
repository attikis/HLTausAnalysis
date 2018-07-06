#================================================================================================
# Script for setting up standalone environment for accessing the python libraries and scripts.
#
# Usage:
# cd HLTausAnalysis
# source setup.csh
#
# Note:
# tested so far LOCATION="" and LOCATION="jade"
#================================================================================================
if ( $?HLTAUSANALYSIS_BASE ) then
    echo "Standalone environment already loaded"
    exit
endif

set LOCATION=""
if ( $?CMSSW_BASE ) then
    set LOCATION="CMSSW"
    echo "=== LOCATION is $LOCATION. Known problems with Python when 'cmsenv' is set. Press any key to continue: "
    set proceed=$<
    echo "Continuing ..."
endif

# Detect LXPLUS or MAC OS X (Darwin)
if ( $LOCATION == "" ) then
    if (`hostname` =~ "lxplus"* ) then
        set LOCATION="lxplus"
    else if (`hostname` =~ "Mac"* ) then
	set LOCATION="mac"
    else if (`hostname` =~ *".cern.ch" ) then #Example: p06109780e53561.cern.ch
	set LOCATION="lxbatch"
    endif
endif


# Set the HLTausAnalysis base directory
setenv HLTAUSANALYSIS_BASE $PWD


set LD_LIBRARY_PATH_APPEND=""
if ( $LOCATION == "lxplus" || $LOCATION == "lxbatch" ) then
    echo "=== Sourcing $LOCATION environments (Hand-picked from CMSSW_7_6_5). To update these follow the instutions below:"
    echo "1) create a developer area (cmsrel):"
    echo "cmsrel CMSWW_X_Y_Z"
    echo "2) source the CMSSW environment (cmsenv):"
    echo "cd  CMSWW_X_Y_Z/src/"
    echo "cmcsev"
    echo "3) Look the new paths with 'scram tool list' and 'scram tool info'"
    echo "[See setup.csh for more details]"
       
    # scram tool info gcc-cxxcompiler (Look for GCC_CXXCOMPILER_BASE)
    set GCC_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3
    #echo "=== GCC_BASE=${GCC_BASE}"

    # scram tool info root_interface (Look for ROOT_INTERFACE_BASE)
    setenv ROOTSYS /cvmfs/cms.cern.ch/slc6_amd64_gcc493/lcg/root/6.02.12-kpegke4
    #echo "=== ROOTSYS=$ROOTSYS"

    # scram tool info xrootd (Look for XROOTD_BASE)
    set XROOTD_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/xrootd/4.0.4-kpegke2
    #echo "=== XROOTD_BASE=${XROOTD_BASE}"

    # scram tool info xz (Look for XZ_BASE)
    set XZ_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/xz/5.2.1
    #echo "=== XZ_BASE=${XZ_BASE}"

    # scram tool info python
    set PYTHON_BASE=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/python/2.7.6-kpegke
    #echo "=== PYTHON_BASE=${PYTHON_BASE}"

    # Set the run-time shared library loader (ld.so) an extra set of directories to look for when searching for shared libraries.
    set LD_LIBRARY_PATH_APPEND=$ROOTSYS/lib:$GCC_BASE/lib64:$GCC_BASE/lib:$XROOTD_BASE/lib:$XZ_BASE/lib:$PYTHON_BASE/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/libjpg/8b-cms/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/libpng/1.6.16/lib
    
    # Tell the shell which directories to search for executable files
    setenv PATH $ROOTSYS/bin:$GCC_BASE/bin:$XROOTD_BASE/bin:$PATH

    # $?VAR expands to 1 (true) if $VAR is set (to anything, even the empty string), 0 (false) if it isn't.
    if ($?PYTHONPATH) then
	setenv PYTHONPATH "$ROOTSYS/lib:$PYTHONPATH"
    else
        setenv PYTHONPATH "$ROOTSYS/lib"
    endif
endif

# Working locally?
if ( $LOCATION == "mac" ) then
    echo "=== Using default ROOT version"
    #source /Users/$USER/ROOT/v5-34-00-patches/bin/thisroot.csh
    #setenv ROOTSYS /Users/$USER/ROOT/v5-34-00-patches/
    #source /Users/$USER/root/v6-10-08/bin/thisroot.csh
    #RootSetup.sh v6-11-03
    source /Users/$USER/root/v6-11-03/bin/thisroot.csh
    if ($?PYTHONPATH) then
	setenv PYTHONPATH "$ROOTSYS/lib:$PYTHONPATH"
    else
	setenv PYTHONPATH "$ROOTSYS/lib"
    endif
endif


echo
echo "=== Appending to LD_LIBRARY_PATH"
set LD_LIBRARY_PATH_APPEND="$HLTAUSANALYSIS_BASE/NtupleAnalysis/lib:${LD_LIBRARY_PATH_APPEND}"
if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH_APPEND}"
else
    setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH_APPEND}:${LD_LIBRARY_PATH}"
endif

echo "=== Creating symbolic links and hidden directories for $LOCATION"
set PATHPREFIX=.python

if ( $LOCATION == "CMSSW" ) then    
    if ( ! $?CMSSW_BASE || ! -e $CMSSW_BASE/python/HLTausAnalysis/NtupleAnalysis ) then
        ln -s $HLTAUSANALYSIS_BASE/NtupleAnalysis/python $CMSSW_BASE/python/HLTausAnalysis/NtupleAnalysis
    endif

else
    if ( ! -e $PATHPREFIX/HLTausAnalysis ) then
	echo "=== Creating $PATHPREFIX directory under `pwd`. Creating __init__.py"
        mkdir -p $PATHPREFIX/HLTausAnalysis
        touch $PATHPREFIX/HLTausAnalysis/__init__.py
    endif

    echo "=== Loop over directories under NtupleAnalysis/ and HeavyChHLTausToTauNu/"
    foreach DIR ( NtupleAnalysis )
	# echo "DIR=$DIR"

	set LINK_NAME=$PATHPREFIX/HLTausAnalysis/$DIR
	set TARGET=$HLTAUSANALYSIS_BASE/$DIR/python
	set PYINIT=$LINK_NAME/__init__.py

	# If $PATHPREFIX/HLTausAnalysis/$DIR does not exist
        if ( ! -e $PATHPREFIX/HLTausAnalysis/$DIR ) then

            echo "Linking $TARGET with $LINK_NAME"
	    ln -s $TARGET $LINK_NAME

	    echo "Creating $PYINIT"
            touch $PYINIT
	    
            foreach d ( $PATHPREFIX/HLTausAnalysis/$DIR/* )
                if ( -d $d ) then
		    echo "Creating $d/__init__.py"
                    touch $d/__init__.py
                endif
            end
        endif
    end


    echo "=== Loop over directories under NtupleAnalysis/src"
    foreach DIR ( `ls NtupleAnalysis/src` )
	# echo "DIR=$DIR"

	# NOTE: Remove last "/" from directory name. The "/" at the end causes the linking to FAIL for some shells
	set DIR=`echo $DIR | sed 's/\(.*\)\//\1 /'`
	set LINK_NAME=$PATHPREFIX/HLTausAnalysis/$DIR
	set TARGET=$HLTAUSANALYSIS_BASE/NtupleAnalysis/src/$DIR/python
	set PYINIT=$LINK_NAME/__init__.py

	# If $LINK_NAME does not exist and $TARGET exists
        if ( ! -e $LINK_NAME && -e $HLTAUSANALYSIS_BASE/NtupleAnalysis/src/$DIR/python ) then

            echo "Linking $TARGET with $LINK_NAME"
            ln -s $TARGET $LINK_NAME

            # echo "Creating $PYINIT"
            touch $PYINIT

            foreach d ( $PATHPREFIX/HLTausAnalysis/$DIR/* )
                if ( -d $d ) then
		    echo "Creating $d/__init__.py"
                    touch $d/__init__.py
                endif
            end
        endif
    end

    # Set PYTHONPATH (-z string is null, that is, has zero length)
    if ( -z PYTHONPATH ) then
        setenv PYTHONPATH ${PWD}/${PATHPREFIX} #NOTE: Double quotes will NOT WORK for some shells!!!
    else
        setenv PYTHONPATH ${PWD}/${PATHPREFIX}:${PYTHONPATH}  #NOTE: Double quotes will NOT WORK for some shells!!!
    endif

endif


#echo "=== Setting PATH variable"
setenv PATH "${HLTAUSANALYSIS_BASE}/NtupleAnalysis/scripts:${PATH}"

# Get the python version
set PYTHON_VERSION=`python -c 'import sys ;version=sys.version_info[:3]; print("{0}.{1}.{2}".format(*version))'`

echo
echo "=== The environment variables set are:"
echo "LOCATION is $LOCATION"
echo "PYTHONPATH is $PYTHONPATH"
echo "HLTAUSANALYSIS_BASE is $HLTAUSANALYSIS_BASE"
echo "PATHPREFIX is $PATHPREFIX"
echo "ROOTSYS is $ROOTSYS"
echo "LD_LIBRARY_PATH is $LD_LIBRARY_PATH"
echo "PYTHONPATH is $PYTHONPATH"
echo "PATH is $PATH"
echo "PYTHON_VERSION is $PYTHON_VERSION"
echo

if ( $PYTHON_VERSION =~ *"2.6"* ) then
    echo "=== WARNING! Requires python 2.7 and later. The python version used is:"
    python -V
    echo

    # Python versions on LXPLUS [https://cern.service-now.com/service-portal/article.do?n=KB0000730]
    if ( $LOCATION == "lxplus" || $LOCATION == "lxbatch" ) then
	echo "=== Printing list of available python versions:"
	scl -l | grep python
	echo

	echo "=== To enable python27 for SHELL=$SHELL copy/paste the following:"
	if ( $SHELL == "/bin/tcsh") then
	    echo "scl enable python27 csh"
	    # scl enable python27 csh
	    echo
	else if ( $SHELL == "/bin/bash") then
	echo "scl enable python27 bash"
	#scl enable python27 bash
	echo
	endif 
    else
	echo "Unsupported SHELL == $SHELL. EXIT shell script"
	exit 1
    endif
else
    echo "=== The python version used is:"
    python -V
    echo 
endif


