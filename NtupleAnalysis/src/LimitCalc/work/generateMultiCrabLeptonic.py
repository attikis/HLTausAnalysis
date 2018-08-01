#!/usr/bin/env python

import HiggsAnalysis.LimitCalc.LandSTools as lands

lepType = False
lhcType = False
lhcTypeAsymptotic = False

#lepType = True
lhcType = True
#lhcTypeAsymptotic = True

mutau = False
etau = False
emu = False

#mutau = True
#etau = True
#emu = True


crabScheduler = "arc"
crabOptions = {
#        "GRID": [
#            "ce_white_list = jade-cms.hip.fi",
#            "ce_white_list = korundi.grid.helsinki.fi",
#            ]
    }

massPoints = lands.allMassPoints
#massPoints = ["155", "160"]
#massPoints = ["160"]

def main(opts):
    postfix = "taujets"


    if opts.mutau:
        generate(opts, lands.mutauDatacardPattern, "mutau",
                 lepRmax={"default": "0.05",
                          "160": "0.3",
                          },
                 lhcRmax="0.5",
                 lhcvR={"default": ("0.02", "0.2"),
                        "150": ("0.03", "0.3"),
                        "155": ("0.03", "0.4"),
                        "160": ("0.03", "0.5"),
                        },
                 lhcScanRmin={"default": None,
                              "160": "0.02",
                              },
                 lhcScanRmax={"default": None,
                              "160": "0.5",
                              },
                 lhcToysCLsb={"default": 600,
                              "160": 300
                              },
                 lhcToysCLb={"default": 300,
                             "160": 150
                             },
                 lhcJobs = {"default": 5,
                            "160": 20
                            },
                 lhcAsyRmax="0.2",
                 lhcAsyOptsExp={"default": lands.lhcAsymptoticOptionsExpected, # for migrad
                                "160": lands.lhcAsymptoticOptionsObserved, # for minos instead of migrad
                                }
                 )

    if opts.etau:
        generate(opts, lands.etauDatacardPattern,  "etau",
                 lepRmax={"default": "0.05",
                          "160": "0.3"
                          },
                 lhcRmax="0.5",
                 lhcvR={"default": ("0.02", "0.25"),
                        "140": ("0.03", "0.3"),
                        "150": ("0.04", "0.4"),
                        "155": ("0.04", "0.4"),
                        "160": ("0.05", "0.6"),
                        },
                 lhcScanRmin={"default": None,
                              "140": "0.02",
                              "155": "0.02",
                              "160": "0.02",
                              },
                 lhcScanRmax={"default": None,
                              "140": "0.4",
                              "155": "0.5",
                              "160": "0.5",
                              },
                 lhcToysCLsb=600, lhcToysCLb=300, lhcPostfix="jobs10_sb_600_b300",
                 lhcJobs = {"default": 20,
                            "155": 40,
                            "160": 40,
                            },
                 lhcAsyRmax="0.2",
                 lhcAsyOptsExp={"default": lands.lhcAsymptoticOptionsExpected, # for migrad
                                "155": lands.lhcAsymptoticOptionsObserved, # for minos
                                "160": lands.lhcAsymptoticOptionsObserved, # for minos
                                }
                 )

    if opts.emu:
        generate(opts, lands.emuDatacardPattern,   "emu",   lepRmax="0.3")

def generate(opts, datacard, postfix,
             lepRmax=None,
             lhcRmax=None, lhcvR=None, lhcToysCLsb=1200, lhcToysCLb=600, lhcJobs=5, lhcPostfix="jobs5_sb1200_b600", lhcScanRmin=None, lhcScanRmax=None,
             lhcAsyRmax=None, lhcAsyvR=None, lhcAsyOptsObs=None, lhcAsyOptsExp=None):
    if opts.lepType:
        lands.generateMultiCrab(
            opts,
            massPoints = massPoints, datacardPatterns = [datacard], rootfilePatterns = [],
            clsType = lands.LEPType(toysPerJob=1000, rMax=lepRmax),
            numberOfJobs = 1,
            postfix = postfix+"_lep_toys1k",
            crabScheduler=crabScheduler, crabOptions=crabOptions)
    if opts.lhcType:
        lands.generateMultiCrab(
            opts,
            massPoints = massPoints, datacardPatterns = [datacard], rootfilePatterns = [],
            clsType = lands.LHCType(toysCLsb=lhcToysCLsb, toysCLb=lhcToysCLb, rMax=lhcRmax, vR=lhcvR, scanRmin=lhcScanRmin, scanRmax=lhcScanRmax),
            numberOfJobs = lhcJobs,
            postfix = postfix+"_lhc_"+lhcPostfix,
            crabScheduler=crabScheduler, crabOptions=crabOptions)
    if opts.lhcTypeAsymptotic:
        lands.produceLHCAsymptotic(
            opts,
            massPoints = massPoints, datacardPatterns = [datacard], rootfilePatterns = [],
            clsType = lands.LHCTypeAsymptotic(rMax=lhcAsyRmax, vR=lhcAsyvR, optionsObserved=lhcAsyOptsObs, optionsExpected=lhcAsyOptsExp),
            postfix = postfix+"_lhcasy"
            )
    

if __name__ == "__main__":
    parser = lands.createOptionParser(lepType, lhcType, lhcTypeAsymptotic)
    parser.add_option("--etau", dest="etau", default=etau, action="store_true",
                      help="Generate for e+tau final state (default %s)" % str(etau))
    parser.add_option("--mutau", dest="mutau", default=mutau, action="store_true",
                      help="Generate for mu+tau final state (default %s)" % str(mutau))
    parser.add_option("--emu", dest="emu", default=emu, action="store_true",
                      help="Generate for e+mu final state (default %s)" % str(emu))

    opts = lands.parseOptionParser(parser)
    main(opts)
