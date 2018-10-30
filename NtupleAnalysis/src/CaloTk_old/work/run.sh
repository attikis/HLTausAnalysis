#!/bin/bash
#====================================================================================================
# DESCRIPTION:
# Written to accommodate the need to run all datasets for a given analyser (CalTk) in one go, instead 
# of the default one-by-one dataset. The idea is to be able to fire-and-forget, with all datasets
# run on automatically and once done a pseudo-multicrab directory also created by running the dedicated
# python script.
#
#
# EXAMPLES:
# ./run.sh <multicrab_path>
# ./run.sh /eos/user/m/mtoumazo/multicrab_HLTaus_v1015_20180710T1650
# ./run.sh /eos/user/m/mtoumazo/multicrab_HLTaus_v1015_20180710T1650 -1 test
#
# LAST USED:
#
#
#
# USEFUL LINKS:
# https://root.cern.ch/root/roottalk/roottalk01/4044.html
# https://root.cern.ch/root/htmldoc/guides/primer/ROOTPrimer.html#n-tuples-in-root# (7.2.4 Processing N-tuples Spanning over Several Files)
#
#====================================================================================================

#====================================================================================================
# Ensure all script arguments are passed from command line
#====================================================================================================
if [ $# -eq 0 ]
then
    echo "=== You must give at least 1 argument (the path to a valid multicrab directory). For example:"
    echo "./run.sh /eos/user/m/mtoumazo/multicrab_HLTaus_v1015_20180710T1650"
    echo
    exit 1
fi

if [ -z "$2" ]
  then
    echo "=== MAXEVENTS argument not found. Using default value of -1"
    MAXEVENTS=-1
else
    # echo "=== MAXEVENTS is ${2}"
    MAXEVENTS=${2}
fi

# if [ -z "$3" ]
#   then
#     echo "=== OUTPUTDIR argument not found. Using default value of multicrab_CaloTkSkim_v92X_<time>"
#     STIME=`date '+%Hh%Mm%Ss_%d%h%Y'`
#     OUTPUTDIR="multicrab_CaloTkSkim_v92X_${STIME}"
# else
#     STIME=`date '+%Hh%Mm%Ss_%d%h%Y'`
#     OUTPUTDIR="${3}_${STIME}"
# fi

#====================================================================================================
# Define Variables
#====================================================================================================
MCRABDIR=${1}
CWD=`pwd`

# Save the submit/start time for future use
if [ -d ${MCRABDIR} ]; then
    echo "=== Directory ${MCRABDIR} found"

    # For-loop: All directories in multicrab dir
    for d in "${MCRABDIR}"/*/; do
	DATASET=`basename "${d}"`
	root -l -b -q "run.cc(\"${MCRABDIR}\", \"${DATASET}\", \"\", ${MAXEVENTS})" &
	# echo "${MCRABDIR} ${DATASET} ${MAXEVENTS}"
	sleep 2
    done
else
    echo "=== Multicrab directory ${MCRABDIR} not found. EXIT"
    exit
fi

# echo "=== Creating pseudo-multicrab directory"
# # For-loop: All root files!
# for f in "${CWD}"/*.root; do
#     FILE=`basename ${f}`
#     echo "${FILE}"
#     hltausPseudoMulticrab.py -f ${FILE} --dir ${OUTPUTDIR}
#     sleep 1
#     done