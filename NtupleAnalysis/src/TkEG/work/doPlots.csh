#!/bin/csh                                                                                                                                                                          
# SingleNeutrino_L1TPU140
# SingleNeutrino_L1TPU200
# SingleTau_L1TnoPU
# SingleTau_L1TPU200
# TT_TuneCUETP8M2T4_14TeV_L1TnoPU
# TT_TuneCUETP8M2T4_14TeV_L1TPU140
# TT_TuneCUETP8M2T4_14TeV_L1TPU200
# ChargedHiggs200_14TeV_L1TnoPU
# ChargedHiggs200_14TeV_L1TPU140
# ChargedHiggs200_14TeV_L1TPU200
# ChargedHiggs500_14TeV_L1TnoPU
# ChargedHiggs500_14TeV_L1TPU140
# ChargedHiggs500_14TeV_L1TPU200
# ChargedHiggs1000_14TeV_L1TnoPU
# ChargedHiggs1000_14TeV_L1TPU140
# ChargedHiggs1000_14TeV_L1TPU200
# GluGluHToTauTau_14TeV_L1TnoPU
# GluGluHToTauTau_14TeV_L1TPU140
# GluGluHToTauTau_14TeV_L1TPU200

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 1) then
    echo "=== You must give exactly 1 argument:"
    echo "1=PSEUDO_MCRAB_DIR"
    echo "\n=== For example:"
    echo "./doPlots.csh multicrab_TkTau_v92X_SeedPt5_SigRMax0p15_IsoRMax0p35_VtxIso0p5_RelIso0p5_22h48m02s_09Oct2018"
    echo
    exit 1
endif

#================================================================================================
# Define variables
#================================================================================================
set INITIAL = `echo $USER | cut -c1-1`
set PSEUDO_MCRAB_DIR = ${1}

#"SingleNeutrino_L1TPU140|SingleNeutrino_L1TPU200|SingleTau_L1TnoPU|SingleTau_L1TPU200|TT_TuneCUETP8M2T4_14TeV_L1TPU140|TT_TuneCUETP8M2T4_14TeV_L1TPU200|ChargedHiggs200_14TeV_L1TnoPU|ChargedHiggs200_14TeV_L1TPU140"
#"ChargedHiggs200_14TeV_L1TPU200|ChargedHiggs500_14TeV_L1TnoPU|ChargedHiggs500_14TeV_L1TPU140|ChargedHiggs500_14TeV_L1TPU200|ChargedHiggs1000_14TeV_L1TnoPU|ChargedHiggs1000_14TeV_L1TPU140|ChargedHiggs1000_14TeV_L1TPU200"
#"GluGluHToTauTau_14TeV_L1TPU140|GluGluHToTauTau_14TeV_L1TPU200"

./plotCounters.py -i 'SingleTau_L1TnoPU|SingleNeutrino_L1TPU200|Glu' -n -m $PSEUDO_MCRAB_DIR --url
./plotTH1.py -n -i "SingleNeutrino|GluGluHToTauTau_14TeV_L1TPU200|SingleTau_L1TPU200|GluGluHToTauTau_14TeV_L1TPU140" -m $PSEUDO_MCRAB_DIR --url
./plotRateVsEff.py -e "SingleE" -m $PSEUDO_MCRAB_DIR --url
./plotTH2.py -e "SingleE" --logZ --normalizeToOne -m $PSEUDO_MCRAB_DIR --url
./plotResolutions.py -i "GluGluHToTauTau_14TeV_L1TPU200" -n -m $PSEUDO_MCRAB_DIR --url
./plotResolutionsPU.py -i "GluGluHToTauTau" -n -m $PSEUDO_MCRAB_DIR  --url
./plotComparison.py -i "SingleTau_L1TPU200" -n -m $PSEUDO_MCRAB_DIR --url

#./plotTkTau.py -n -e "SingleE|Charged|TT" -m $PSEUDO_MCRAB_DIR --url
#./plotRateVsEff.py -e "SingleE" -m $PSEUDO_MCRAB_DIR --url
#./plotTH2.py -e "SingleE" --logZ --normalizeToOne -m $PSEUDO_MCRAB_DIR --url
#./plotResolutions.py -e "SingleE" -n -m $PSEUDO_MCRAB_DIR --url

