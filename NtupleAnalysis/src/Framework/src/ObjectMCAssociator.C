#ifndef ObjectMCAssociator_cxx
#define ObjectMCAssociator_cxx

// User
#include "../interface/ObjectMCAssociator.h"

//#define myDEBUG

short int ObjectMCAssociator::FindTauMCProvenance(const unsigned short iTau)
{
  int idGenTau = GenP_PdgId.at(iTau);
  if (abs(idGenTau) != 15 ) {
    std::cout << " E R R O R ! The Method FindTauMCProvenance can be called"
              << " only for generated muons \n";
    std::cout << " The id of the tau passed as argument is instead "<< idGenTau
	      << std::endl;
    exit(1);
  }
  short int tauSign = idGenTau > 0 ? 1 : -1;
  return tauSign;
}

#endif //ObjectMCAssociator_cxx
