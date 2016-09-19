#!/bin/bash

if [[ $1 == -h || $1 == --help ]]
then 
  echo Instructions:
  echo "Invoke the script with MACRONAME MODE"
  echo "   MACRONAME is the name of the new macro to create"
  echo "   MODE can be either:"
  echo "      \"reco\" to analyze ntuples containing only reco information"
  echo "               (Real data or MC)"
  echo "       \"MC\"  to analyze ntuples containing reco and generator"
  echo "               information (MC only)"
  echo "EXAMPLE: . creatNewMacro MyMacro reco"
  echo ""
  echo "The script will create a directory called MACRONAME"
  echo "and 3 files in it:"
  echo "   MACRONAME.h is the header file of the new nacro"
  echo "   MACRONAME.C is the file with the Loop() function"
  echo "       (this is the only file that the user normally"
  echo "        needs to modify)"
  echo "   runMACRONAME.cc is the file to be called to execute the macro"
  echo ""
  echo "MACRO EXECUTION"
  echo "==============="
  echo "Macro is executed by running the file runMACRONAME.cc from within ROOT:"
  echo "   root [0] .x runMACRONAME.cc(\"SAMPLENAME\",\"TEXT\",NEvts)"
  echo "The command arguments are:"
  echo "   SAMPLENAME is the name of the sample to analyze"
  echo "   TEXT is a free text (this is available from within the macro class,"
  echo "        in the variable \"std::string text\")"
  echo "   NEvts is the number of events to analyze (if NEvts = -1 => all evts)"
  echo ""
  echo "INPUT FILES"
  echo "==========="
  echo "Input sample must be defined in the form of a text file, containing"
  echo "a list of ROOT ntuple files. The file must be located in the fileLists" 
  echo "directory and its name must follow this convention:"
  echo "   fileList_SAMPLENAME.txt"
  echo "The string SAMPLENAME is the same that is used in the macro invocation."
  echo "Example files containing lists of ntuples in local or remote locations"
  echo "are already in the fileLists directory. It is suggested to be as "
  echo "descriptive as possible in the naming of the file lists. "
  echo "In particular, the sample name should always begin with the string"
  echo "MC_\" if it is a simulated one, or with the string \"Data_\" if it"
  echo "contains real data."
  echo ""
  echo "OUTPUT FILES"
  echo "============"
  echo "Output file for histograms, canvases, etc, is automatically managed"
  echo "Every time the macro is invoked, a new output file is created"
  echo "with the name: MACRONAME_Histograms_SAMPLENAME_TEXT.root"
  echo "The file is overwritten if it already existed in the working directory."
  echo ""
  echo "HISTOGRAM BOOKING AND PLOTTING"
  echo "=============================="
  echo "A utility class to book and plot, with good-looking styles consistent"
  echo "across all the macros, 1D and 2D histograms and graphs, W I L L   B E"
  echo "provided in utilities/histoPlotter.C[h]."
  echo ""
  echo "EVENT AND PHYSICS OBJECT SELECTION"
  echo "=================================="
  echo "A utility class to define and apply selections in a consistent way"
  echo "across all the macros is provided in utilities/OjectSelector.C[h]."
  echo "The class will be augmented in the near future, however some" 
  echo "functions are already available" 
  echo ""
  echo "PHYSICAL AND MATHEMATICAL CONSTANTS"
  echo "==================================="
  echo "A file in which all the constants definitions can be collected is "
  echo "provided in utilities/constants.h."
  echo "Constants defined in this file can be visible to all the macros."
  echo "Within a macro, a constant is used via  \"constants::CONSTANTNAME\"."
  echo "To avoid to use \"constants::\" before the constant name, use instead"
  echo "   using namespace constants;"
  echo "in MACRONAME.C, just before the beginning of the Loop() function."
  return
fi

if [[ $2 == '' ]]
then 
    echo Script accepts exactly two arguments! Write \"--help\" for help.
    return
fi

if [[ $2 == 'reco' ||  $2 == 'MC' ]]
then 
    echo Creating tree analyzer for $2
else 
    echo Second argument must be \"reco\", \"MC\"!
    return
fi

macroName=$1
mode=$2

echo Name of the macro to be created: $macroName
echo

if [[ -e $macroName ]]
then 
    echo ERROR: directory $macroName exists!
    echo Aborting macro creation!
    return
fi

echo Creating directory $macroName
mkdir $macroName

if [[ $mode == "MC" ]]
then
    echo Creating $macroName/$macroName.h...
    if [[ -e $macroName/$macroName.h ]]
    then 
	echo ERROR: file $macroName/$macroName.h exists!
	echo Aborting macro creation!
	return
    fi
    sed "s#MACRONAMEMC#$macroName#g" Templates/template_MacroMC.h > $macroName/$macroName.h
    echo Creating $macroName/$macroName.C...
    if [[ -e $macroName/$macroName.C ]]
    then 
	echo ERROR: file $macroName/$macroName.C exists!
	echo Aborting macro creation!
	return
    fi

    sed "s#MACRONAMEMC#$macroName#g" Templates/template_MacroMC.C > $macroName/$macroName.C
    echo Creating $macroName/run$macroName.cc...
    if [[ -e $macroName/run$macroName.cc ]]
    then
	echo ERROR: file $macroName/run$macroName.cc exists!
	echo Aborting macro creation!
	return
    fi

    sed "s#MACRONAMEMC#$macroName#g" Templates/template_runMacroMC.cc > $macroName/run$macroName.cc
    echo File creation completed! 
    echo You can find the macro skeleton in directory ./$macroName


elif [[ $mode == "reco" ]]
then
    echo Creating $macroName/$macroName.h...
    if [[ -e $macroName/$macroName.h ]]
    then 
	echo ERROR: file $macroName/$macroName.h exists!
	echo Aborting macro creation!
	return
    fi

    sed "s#MACRONAMERECO#$macroName#g" Templates/template_MacroReco.h > $macroName/$macroName.h
    echo Creating $macroName/$macroName.C...
    if [[ -e $macroName/$macroName.C ]]
    then 
	echo ERROR: file $macroName/$macroName.C exists!
	echo Aborting macro creation!
	return
    fi

    sed "s#MACRONAMERECO#$macroName#g" Templates/template_MacroReco.C > $macroName/$macroName.C
    echo Creating $macroName/run$macroName.cc...
    if [[ -e $macroName/run$macroName.cc ]]
    then 
	echo ERROR: file $macroName/run$macroName.cc exists!
	echo Aborting macro creation!
	return
    fi

    sed "s#MACRONAMERECO#$macroName#g" Templates/template_runMacroReco.cc > $macroName/run$macroName.cc
    echo File creation completed! You can find the macro skeleton in directory ./$macroName

else 
    echo ERROR: unrecognized mode! Exiting.
    return
fi
