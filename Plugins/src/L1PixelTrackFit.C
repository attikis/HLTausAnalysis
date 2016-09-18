#ifndef L1PixelTrackFit_cxx
#define L1PixelTrackFit_cxx

#include "../interface/L1PixelTrackFit.h"


//****************************************************************************
L1PixelTrackFit::L1PixelTrackFit(double Bz)
//****************************************************************************
{

  // Set the magnetic field to be used for pT calculation.
  // Assume that the B-field is uniform inside the magnet solenoid (z-component)
  mMagneticFieldStrength = Bz;
  
}


//****************************************************************************
L1PixelTrackFit::~L1PixelTrackFit()
//****************************************************************************
{

  _InitVars();
  
}


//****************************************************************************
void L1PixelTrackFit::_InitVars(void)
//****************************************************************************
{

  pixHitL1.clear();
  pixHitL2.clear();
  pixHitL3.clear();
  pixHitL4.clear();
  pixHitD1.clear();
  pixHitD2.clear();
  pixHitD3.clear();
  rinvfit   = 0.0;
  phi0fit   = 0.0;
  d0fit     = 0.0;
  tfit      = 0.0;
  z0fit     = 0.0;
  chisqfit  = 0.0;
  sigmarinv = 0.0;
  sigmaphi0 = 0.0;
  sigmad0   = 0.0;
  sigmat    = 0.0;
  sigmaz0   = 0.0;
  pixHits.clear();
  candidatePixHits.clear();
  cand_sorted = false;

  return;
  
}


//****************************************************************************
void L1PixelTrackFit::_SetCandPixHitContainers(vector<double> candPixHits_X,
					       vector<double> candPixHits_Y,
					       vector<double> candPixHits_Z,
					       vector< int >  candPixHits_Type)
//****************************************************************************
{
  
  const int nHits = candPixHits_Type.size();
  
  // For-loop: All Pixel Hits
  for(int i = 0; i < nHits; i++){
    
    // Create 3D space point for current pixel hit
    ROOT::Math::XYZVector candPixHit_XYZ( candPixHits_X.at(i), candPixHits_Y.at(i), candPixHits_Z.at(i) );
    candidatePixHits.push_back( candPixHit_XYZ );

    // Fill the correct container according to the pixel hit type. Layer: +1, +2, +3, +4. Disk : -1, -2, -3
    if (candPixHits_Type.at(i) == +1)      pixHitL1.push_back(candPixHit_XYZ);
    else if (candPixHits_Type.at(i) == +2) pixHitL2.push_back(candPixHit_XYZ);
    else if (candPixHits_Type.at(i) == +3) pixHitL3.push_back(candPixHit_XYZ);
    else if (candPixHits_Type.at(i) == +4) pixHitL4.push_back(candPixHit_XYZ);
    else if (candPixHits_Type.at(i) == -1) pixHitD1.push_back(candPixHit_XYZ);
    else if (candPixHits_Type.at(i) == -2) pixHitD2.push_back(candPixHit_XYZ);
    else if (candPixHits_Type.at(i) == -3) pixHitD3.push_back(candPixHit_XYZ);
    else{
      std::cout << "E R R O R ! L1PixelTrackFit::L1PixelTrackFit(...) - Invalid pixel hit type \"" << candPixHits_Type.at(i) << "\". EXIT" << std::endl;
      exit(1);
    }
    
  }// for(int i = 0; i < nHits; i++){
  
  return;
}



//****************************************************************************
void L1PixelTrackFit::Multifit(double rinv,
			       double phi0,
			       double d0,
			       double t,
			       double z0,
			       std::vector<ROOT::Math::XYZVector> hitL1,
			       std::vector<ROOT::Math::XYZVector> hitL2,
			       std::vector<ROOT::Math::XYZVector> hitL3,
			       std::vector<ROOT::Math::XYZVector> hitL4,
			       std::vector<ROOT::Math::XYZVector> hitD1,
			       std::vector<ROOT::Math::XYZVector> hitD2,
			       std::vector<ROOT::Math::XYZVector> hitD3,		
			       bool& success,
			       double& invrfinal,
			       double& phi0final,
			       double& d0final,
			       double& tfinal,
			       double& z0final,
			       double& chisqfinal,
			       int& nhit,
			       double &sigmarinv,
			       double &sigmaphi0,
			       double &sigmad0,
			       double &sigmat,
			       double &sigmaz0,
			       std::vector<ROOT::Math::XYZVector> &pixHits,
			       unsigned int cfg_NPixelHits_Min)
//****************************************************************************
{

  // Boolean for a sucessful track fit to pixel hits
  success = false;

  // Create an array of vectors to hold the pixel hits (descending order wrt occupancy)
  std::vector<ROOT::Math::XYZVector> hits[7];
  bool barrel[7];
  
  barrel[0] = true;
  barrel[1] = true;
  barrel[2] = true;
  barrel[3] = true;
  barrel[4] = false;
  barrel[5] = false;
  barrel[6] = false;

  hits[0] = hitL1;
  hits[1] = hitL2;
  hits[2] = hitL3;
  hits[3] = hitL4;
  hits[4] = hitD1;
  hits[5] = hitD2;
  hits[6] = hitD3;

  
  // Sort the number of hits per layer
  bool more=false;
  do {
    more=false;
     for(int i=0;i<6;i++) {
       if (hits[i].size()<hits[i+1].size()) {
	 more=true;
	 std::vector<ROOT::Math::XYZVector> tmp=hits[i];
	 hits[i]=hits[i+1];
	 hits[i+1]=tmp;
	 bool tmpb=barrel[i];
	 barrel[i]=barrel[i+1];
	 barrel[i+1]=tmpb;
       }
     }
  } while(more);
  

  bool bPrintBeforeSort = false;
  bool bPrintAfterSort  = false;
  // Print hit containers (before sorting)
  if(bPrintBeforeSort){
    
    std::cout  << "*** Before sort:" << std::endl;
    for (int i=0;i<7;i++) {
      int nHits = (int) hits[i].size();
      
      std::cout  << "hits[" << i << "]:";
      if(barrel[i]) std::cout << " (" << nHits << " barrel hits)" << std::endl;
      else std::cout << " (" << nHits << " disk hits)" << std::endl;
      for (int j=0; j < nHits ; j++){
	std::cout << "\t"<< j << ") r = " << hits[i].at(j).Rho() << ", z = " << hits[i].at(j).z() << ", Phi = " << hits[i].at(j).Phi() << std::endl; 
      }
      std::cout << "" << std::endl;
    }
  }
  
  // Print hit containers (after sorting)
  if(bPrintAfterSort){
    std::cout  << "*** After sort:" << std::endl;
    for (int i=0;i<7;i++) {
      int nHits = (int) hits[i].size();
      std::cout  << "hits[" <<i<<"]:"; 
      if(barrel[i]) std::cout << " (" << nHits << " barrel hits)" << std::endl;
      else std::cout << " (" << nHits << " disk hits)"<< std::endl;
      for (int j=0; j < nHits ; j++){
	std::cout << "\t"<< j << ") r = " << hits[i].at(j).Rho() << ", z = " << hits[i].at(j).z() << ", Phi = " << hits[i].at(j).Phi() << std::endl; 
      }
      std::cout << "" << std::endl;
    }
  }
  
  // Now start fitting
  double bestChisqdof    = 1e30;
  int i0best             = -1;
  int i1best             = -1;
  int i2best             = -1;
  int i3best             = -1;
  int i0temp, i1temp, i2temp, i3temp;
  unsigned int nhitsbest = 0;
  double rinvbest        = 9999.9;
  double phi0best        = 9999.9;
  double d0best          = 9999.9;
  double tbest           = 9999.9;
  double z0best          = 9999.9;
  double chisqbest       = 9999.9;

  
  // For-loop over sorted (CANDIDATE) pixel hit containers (most occupied Layer/Disk at hits[0])
  // Try all possible permutations. The permutation that gives the best chi-square value is the chosen one
  for(unsigned int i0=0; i0<hits[0].size()+1; i0++) {
    for(unsigned int i1=0; i1<hits[1].size()+1; i1++) {
      for(unsigned int i2=0; i2<hits[2].size()+1; i2++) {
	for(unsigned int i3=0; i3<hits[3].size()+1; i3++) {
	  
	  // Count the number of pixel hits
	  unsigned int npixel=0;
	  if (i0!=hits[0].size()) npixel++;
	  if (i1!=hits[1].size()) npixel++;
	  if (i2!=hits[2].size()) npixel++;
	  if (i3!=hits[3].size()) npixel++;
	     
	  // Require at least 3 pixel hits
          if (npixel < cfg_NPixelHits_Min) continue;
	  
	  // Initialise Hit Vectors (passed for fitting)
	  std::vector<ROOT::Math::XYZVector> fithits;
	  std::vector<bool> fitbarrel;
	  
	  // Fill first the most active layer
	  if (i0<hits[0].size()) {
	    i0temp = i0;
	    fithits.push_back(hits[0].at(i0));
	    fitbarrel.push_back(barrel[0]);
	  }
	  else i0temp = -1; 
	  // Fill the second most active layer
	  if (i1<hits[1].size()) {
	    i1temp = i1;
	    fithits.push_back(hits[1].at(i1));
	    fitbarrel.push_back(barrel[1]);
	  }
	  else i1temp = -1; 
	  // Fill the third most active layer
	  if (i2<hits[2].size()) {
	    i2temp = i2;
	    fithits.push_back(hits[2].at(i2));
	    fitbarrel.push_back(barrel[2]);
	  }
	  else i2temp = -1;
	  // Fill the fourth most active layer
	  if (i3<hits[3].size()) {
	    i3temp = i3;
	    fithits.push_back(hits[3].at(i3));
	    fitbarrel.push_back(barrel[3]);
	  }
	  else i3temp = -1;

	  // Perform the actual track fit to the given pixel hit permutation. Return the TTPixelTrack track & fit parameters.
	  double rinvfit, phi0fit, d0fit, tfit, z0fit, chisqfit;
	  TrackFit(rinv, phi0, d0, t, z0,
		   rinvfit, phi0fit, d0fit, tfit, z0fit, chisqfit,
		   sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0,
		   fithits, fitbarrel);
	  double chisqdof = chisqfit/(2.0*npixel-5.0);

	  SaveHitPermutationCandidate(i0temp, i1temp, i2temp, i3temp, chisqfit, chisqdof, rinvfit, phi0fit, d0fit, tfit, z0fit);	  

	  // If the reduced chi-square of this combination is better, select the current pixel combination for pixel track re-fit
	  if (chisqdof < bestChisqdof) {
	  // if ( fabs(chisqdof-1.0)  < fabs(bestChisqdof-1.0)) {
	      
	    bestChisqdof = chisqdof;
	    nhitsbest    = npixel;
	    i0best       = i0temp;
	    i1best       = i1temp;
	    i2best       = i2temp;
	    i3best       = i3temp;
	    rinvbest     = rinvfit;
	    phi0best     = phi0fit;
	    d0best       = d0fit;
	    tbest        = tfit;
	    z0best       = z0fit;
	    chisqbest    = chisqfit;

	  } // if (chisqdof<bestChisqdof) {
	  
	} // for(unsigned int i3=0; i3<hits[3].size()+1; i3++) {
      } // for(unsigned int i2=0; i2<hits[2].size()+1; i2++) {
    } // for(unsigned int i1=0; i1<hits[1].size()+1; i1++) {
  } // for(unsigned int i0=0; i0<hits[0].size()+1; i0++) {


  if (bestChisqdof<1e29) {
    if (0) {
      std::cout << i0best<<i1best<<i2best<<i3best<<std::endl;
    }
    
    success    = true;

    invrfinal  = rinvbest;
    phi0final  = phi0best;
    d0final    = d0best;
    tfinal     = tbest;
    z0final    = z0best;
    chisqfinal = chisqbest;
    nhit       = nhitsbest;

    // SaveHitPermutationCandidate(i0best, i1best, i2best, i3best, chisqfinal, bestChisqdof, invrfinal, phi0final, d0final, tfinal, z0final, "Selected"); //duplicate!
    if (0) PrintHitPermutationCandidate();
    
    // Save hits in array of std::vector
    if (0) std::cout << "*** Fitting: " << i0best << ", " << i1best << ", " << i2best << ", " << i3best << std::endl;
    if (i0best > -1 && i0best < (int) hits[0].size() ) pixHits.push_back(hits[0].at(i0best) );
    if (i1best > -1 && i1best < (int) hits[1].size() ) pixHits.push_back(hits[1].at(i1best) );
    if (i2best > -1 && i2best < (int) hits[2].size() ) pixHits.push_back(hits[2].at(i2best) );
    if (i3best > -1 && i3best < (int) hits[3].size() ) pixHits.push_back(hits[3].at(i3best) );

  }
  
  return;
}


//****************************************************************************
TTPixelTrack L1PixelTrackFit::FitPixelTrack(double rinv,
					    double phi0,
					    double d0,
					    double t,
					    double z0,
					    vector<double> candPixHits_X,
					    vector<double> candPixHits_Y,
					    vector<double> candPixHits_Z,
					    vector< int >  candPixHits_Type)
//****************************************************************************
{
  
  // Make sure the pixel hit containers are empty
  _InitVars();
  
  // Create an array of vectors to hold the pixel hits (descending order wrt occupancy)
  _SetCandPixHitContainers(candPixHits_X, candPixHits_Y, candPixHits_Z, candPixHits_Type);

  // Create a TTPixelTrack
  TTPixelTrack aTrack;  
  bool success = false;
  double invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit;
  double sigmainvr, sigmaphi0, sigmad0, sigmat, sigmaz0;
  int nhit;
  
  // Anders "hack"
  if (fabs(d0) < 10.0){ 

    // XENIOS FIXME ATTIKIS
    // 1) Default-Start
    // Perform the TTrack fit to the candidate pixel hits. Select the best (chi2-driven) pixel hit permutation
     Multifit(rinv, phi0, d0, t, z0,
     	     pixHitL1, pixHitL2, pixHitL3, pixHitL4, pixHitD1, pixHitD2, pixHitD3,
     	     success, invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit, nhit,
     	     sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, 3);
    // 1) Default-End
         

    // 2) AtLeast4Not-Start
    // Multifit(rinv, phi0, d0, t, z0,
    //   	     pixHitL1, pixHitL2, pixHitL3, pixHitL4, pixHitD1, pixHitD2, pixHitD3,
    //    	     success, invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit, nhit,
    //    	     sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, 4);
    //
    //  if (success) { // rejects 4-hit tracks
    //    Multifit(rinv, phi0, d0, t, z0,
    //  	       pixHitL1, pixHitL2, pixHitL3, pixHitL4, pixHitD1, pixHitD2, pixHitD3,
    //  	       success, invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit, nhit,
    // 		sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, 1e2);
    //   }
    //  else{
    //    Multifit(rinv, phi0, d0, t, z0,
    //  	       pixHitL1, pixHitL2, pixHitL3, pixHitL4, pixHitD1, pixHitD2, pixHitD3,
    //  	       success, invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit, nhit,
    //  	       sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, 3);
    //  }
     // 2) AtLeast4Not-End
    
    
    // 3) Priority4-Start
    // Multifit(rinv, phi0, d0, t, z0,
    //   	     pixHitL1, pixHitL2, pixHitL3, pixHitL4, pixHitD1, pixHitD2, pixHitD3,
    //   	     success, invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit, nhit,
    //   	     sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, 4);

    // // If no 4 hits fit possible re-do enabling 3 hits fits
    // if(!success)
    //   { 
    // 	Multifit(rinv, phi0, d0, t, z0,
    // 		 pixHitL1, pixHitL2, pixHitL3, pixHitL4, pixHitD1, pixHitD2, pixHitD3,
    // 		 success, invrfit, phi0fit, d0fit, tfit, z0fit, chisqfit, nhit,
    // 		 sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, 3);
    // }
    // 3) Priority4-Start    
    
  } // if (fabs(d0) < 10.0){
  
  
  // Create the Pixel Track. Initialise it only if the fit is successful
  if (!success) return aTrack;

  // Initialise the Pixel Track
  double pt = 0.299792 * mMagneticFieldStrength / ( 100*fabs(invrfit) );
  ROOT::Math::XYZVector thePOCA(-d0fit * sin(phi0fit), d0fit * cos(phi0fit), z0fit );
  TVector3 theMomentum;
  // TVector3 theMomentum(GlobalVector::Cylindrical(pt, phi0fit, pt*tfit) ); // r, phi, z
  theMomentum.SetPtEtaPhi(pt, asinh(tfit), phi0fit ); // t = sinh(eta)
  aTrack.init(theMomentum, thePOCA, invrfit, chisqfit, nhit, sigmainvr, sigmaphi0, sigmad0, sigmat, sigmaz0, pixHits, candidatePixHits);

  return aTrack;
}
  

//****************************************************************************
void L1PixelTrackFit::TrackFit(double rinv,
			       double phi0,
			       double d0,
			       double t,
			       double z0,
			       double& rinvfit,
			       double& phi0fit,
			       double& d0fit,
			       double& tfit,
			       double& z0fit,
			       double& chisqfit,
			       double& sigmarinv,
			       double& sigmaphi0,
			       double& sigmad0,
			       double& sigmat,
			       double& sigmaz0, 
			       std::vector<ROOT::Math::XYZVector> fithits,
			       std::vector<bool> fitbarrel)
//****************************************************************************
{
  double D[5][8];
  double MinvDt[5][8];  

  CalculateDerivatives(rinv, phi0, t, z0, fithits, fitbarrel, D, MinvDt);
  LinearTrackFit(rinv, phi0, d0, t, z0, rinvfit, phi0fit, d0fit, tfit, z0fit, chisqfit, sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0, fithits, fitbarrel ,D, MinvDt);

  return;
}


//****************************************************************************
void L1PixelTrackFit::CalculateDerivatives(double rinv,
					   double phi0,
					   double t,
					   double z0,
					   std::vector<ROOT::Math::XYZVector> fithits,
					   std::vector<bool> fitbarrel,
					   double D[5][8],
					   double MinvDt[5][8])
//****************************************************************************
{

  // Create matrix for the 5 track parameters: rinv, phi0, t, z0, d0
  double M[5][10];
  unsigned int n = fithits.size();
  assert(n<=4);
  int j=0;
    
  // Note: D[5][8] is the derivative matrix. It has 5 rows and 8 columns.
  // 5 rows due to 5 track parameters (rinv, phi0, t, z0, d0)
  // 8 columns because each of the 4 hits has two paramaters it depends on:
  // a) barrel: rphi, z
  // b) disk  :   r, phi
  // Matrix looks like (asssuming 4 barrel hits):
  //   dphidrinv_1  dzdrinv1   dphidrinv_2  dzdrinv2   dphidrinv_3  dzdrinv3   dphidrinv_4  dzdrinv4
  //   dphidphi0_1  dzdphi0_1  dphidphi0_2  dzdphi0_2  dphidphi0_3  dzdphi0_3  dphidphi0_4  dzdphi0_4
  //   dphidt_1     dzdt_1     dphidt_2     dzdt_2     dphidt_3     dzdt_3     dphidt_4     dzdt_4
  //   dphidz0_1    dzdz0_1    dphidz0_2    dzdz0_2    dphidz0_3    dzdz0_3    dphidz0_4    dzdz0_4
  //   dphidd0_1    dzdd0_1    dphidd0_2    dzdd0_2    dphidd0_3    dzdd0_3    dphidd0_4    dzdd0_4

  
  // For-loop: Pixel Hits (3 or 4) and fill the Derivative Matrix D[5][8]
  for(unsigned int i=0; i < n; i++) {

    double ri   = fithits[i].Rho();
    double zi   = fithits[i].z();
    double phii = fithits[i].Phi();

    double sigmax = 0.007/sqrt(12.0);
    double sigmaz = 0.01/sqrt(12.0);
    
    // Handle the barrel hits (r is known for each layer)
    if (fitbarrel[i]){
            
      // Track Derivatives (rphi)
      double dphidrinv  = -0.5 * ri * ri/sqrt(1 - 0.25 * ri * ri * rinv * rinv);
      double drphidphi0 =  ri;
      double drphidt    =  0.0;
      double drphidz0   =  0.0;
      double drphidd0   = -1.0;

      // Track Derivatives (z)
      double dzdrinv = 0.0;
      double dzdphi0 = 0.0;
      double dzdt    = (2/rinv) * asin(0.5 * ri * rinv);
      double dzdz0   = 1.0;
      double dzdd0   = 0.0;

      // First we have the phi position
      D[0][j] = dphidrinv /sigmax;
      D[1][j] = drphidphi0/sigmax;
      D[2][j] = drphidt   /sigmax;
      D[3][j] = drphidz0  /sigmax;
      D[4][j] = drphidd0  /sigmax;
      j++;

      // Second the z position
      D[0][j]=0.0;
      D[1][j]=0.0;
      D[2][j]=(2/rinv)*asin(0.5*ri*rinv)/sigmaz;
      D[3][j]=1.0/sigmaz;
      D[4][j]=0.0;


      D[0][j] = 0.0;
      D[1][j] = 0.0;
      D[2][j] = dzdt /sigmaz;
      D[3][j] = dzdz0/sigmaz;
      D[4][j] = 0.0;
      j++;
    }
    else { // Handle the disk hits (z is known for each disk)
      
      // First we have the r position
      double r_track   = 2.0 * sin(0.5 * rinv * (zi-z0)/t)/rinv;
      double phi_track = phi0 - 0.5 * rinv * (zi-z0)/t; 
      
      double rmultiplier   = sin(phi_track - phii);
      double phimultiplier = r_track * cos(phi_track-phii);

      // Track Derivatives (r)      
      double drdrinv = -2.0 * sin(0.5 * rinv * (zi-z0)/t)/(rinv * rinv) + (zi-z0) * cos(0.5 * rinv * (zi-z0)/t)/(rinv * t);
      double drdphi0 =  0.0;
      double drdt    = -(zi-z0) * cos(0.5 * rinv * (zi-z0)/t)/(t * t);
      double drdz0   = -cos(0.5 * rinv * (zi-z0)/t)/t;
      double drdd0   = 0.0;
      
      // Track Derivatives (phi)
      double dphidrinv = -0.5*(zi-z0)/t;
      double dphidphi0 = 1.0;
      double dphidt    = 0.5*rinv*(zi-z0)/(t*t);
      double dphidz0   = 0.5*rinv/t;
      double dphidd0   = -1.0;
	
      D[0][j] = drdrinv/sigmaz;
      D[1][j] = drdphi0/sigmaz;
      D[2][j] = drdt   /sigmaz;
      D[3][j] = drdz0  /sigmaz;
      D[4][j] = drdd0  /sigmaz;
      j++;

      // Second the rphi position
      D[0][j] = (phimultiplier * dphidrinv + rmultiplier * drdrinv)/sigmax;
      D[1][j] = (phimultiplier * dphidphi0 + rmultiplier * drdphi0)/sigmax;
      D[2][j] = (phimultiplier * dphidt    + rmultiplier * drdt)   /sigmax;
      D[3][j] = (phimultiplier * dphidz0   + rmultiplier * drdz0)  /sigmax;
      D[4][j] = -1.0/sigmax;
      j++;
    }

    // std::cout << "Exact rinv derivative: " << i << " " << D[0][j-2] << " " << D[0][j-1] << std::endl;
    // std::cout << "Exact phi0 derivative: " << i << " " << D[1][j-2] << " " << D[1][j-1] << std::endl;
    // std::cout << "Exact t derivative   : " << i << " " << D[2][j-2] << " " << D[2][j-1] << std::endl;
    // std::cout << "Exact z0 derivative  : " << i << " " << D[3][j-2] << " " << D[3][j-1] << std::endl;
	
  }
    
  unsigned int npar=5;

  // For-loop: Track Parameters
  for(unsigned int i1 = 0; i1 < npar; i1++){
    for(unsigned int i2 = 0; i2 < npar; i2++){
      M[i1][i2]=0.0;
      // For-loop: Fit Hits (x2)
      for(unsigned int j=0; j < 2*n; j++){
	M[i1][i2] += D[i1][j] * D[i2][j];
      }
    }// for(unsigned int i2 = 0; i2 < npar; i2++){
  }// for(unsigned int i1 = 0; i1 < npar; i1++){
  
  // Approximate errors from L1 tracks
  M[0][0] += 1.0/pow(0.03*rinv,2);
  M[1][1] += 1.0/pow(0.0005,2);
  M[2][2] += 1.0/pow(0.0025,2);

  // Invert the Matrix 
  Invert(M, npar);

  // For-loop: Fit Hits (x2)
  for(unsigned int j=0; j < 2*n; j++) {
    // For-loop: Track Parameters
    for(unsigned int i1=0; i1 < npar; i1++) {
      MinvDt[i1][j]=0.0;
      // For-loop: Track Parameters
      for(unsigned int i2=0; i2 < npar; i2++) {
	MinvDt[i1][j]+=M[i1][i2+npar]*D[i2][j];
      }
    }
  }// for(unsigned int j=0; j < 2*n; j++) {

  return;
}



//****************************************************************************
void L1PixelTrackFit::Invert(double M[5][10],
			     unsigned int n)
//****************************************************************************
{ 
  
  assert(n<=5);
  
  unsigned int i,j,k;
  double ratio,a;
  
  for(i = 0; i < n; i++){
    for(j = n; j < 2*n; j++){
      if(i==(j-n))
	M[i][j] = 1.0;
      else
	M[i][j] = 0.0;
    }
  }
  
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if(i!=j){
	ratio = M[j][i]/M[i][i];
	for(k = 0; k < 2*n; k++){
	  M[j][k] -= ratio * M[i][k];
	}
      }
    }
  }
  
  for(i = 0; i < n; i++){
    a = M[i][i];
    for(j = 0; j < 2*n; j++){
      M[i][j] /= a;
    }
  }

  return;
}


//****************************************************************************
void L1PixelTrackFit::LinearTrackFit(double rinv,
				     double phi0,
				     double d0,
				     double t,
				     double z0,
				     double& rinvfit,
				     double& phi0fit,
				     double& d0fit,
				     double& tfit,
				     double& z0fit,
				     double& chisqfit,
				     double& sigmarinv,
				     double& sigmaphi0,
				     double& sigmad0,
				     double& sigmat,
				     double& sigmaz0, 
				     std::vector<ROOT::Math::XYZVector> fithits,
				     std::vector<bool> fitbarrel,
				     double D[5][8],
				     double MinvDt[5][8])
//****************************************************************************
{
  
  unsigned int n= fithits.size();
  
  //Next calculate the residuals
  double delta[40];
  double chisq=0;
  unsigned int j=0;

  for(unsigned int i=0;i<n;i++) {
    double ri=fithits[i].Rho();
    double zi=fithits[i].z();
    double phii=fithits[i].Phi();
    double sigmax=0.007/sqrt(12.0);
    double sigmaz=0.01/sqrt(12.0);
        
    if (fitbarrel[i]) {
      //we are dealing with a barrel stub
      static const double two_pi=8.0*atan(1.0);

      double deltaphi=phi0-asin(0.5*ri*rinv)+d0/ri-phii;
      if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
      if (deltaphi<-0.5*two_pi) deltaphi+=two_pi;
      assert(fabs(deltaphi)<0.1*two_pi);

      delta[j++]=-ri*deltaphi/sigmax;
      delta[j++]=(z0+(2.0/rinv)*t*asin(0.5*ri*rinv)-zi)/sigmaz;     
    }
    else {
      //we are dealing with a disk hit
      double r_track=2.0*sin(0.5*rinv*(zi-z0)/t)/rinv;
      double phi_track=phi0-0.5*rinv*(zi-z0)/t+d0/ri;
      
      double Delta=r_track*sin(phi_track-phii);

      delta[j++]=(r_track-ri)/sigmaz;
      delta[j++]=-Delta/sigmax;      
    }
    
    chisq+=(delta[j-2]*delta[j-2]+delta[j-1]*delta[j-1]);
  }
  
  double drinv=0.0;
  double dphi0=0.0;
  double dd0=0.0;
  double dt=0.0;
  double dz0=0.0;

  double drinv_cov=0.0;
  double dphi0_cov=0.0;
  double dd0_cov=0.0;
  double dt_cov=0.0;
  double dz0_cov=0.0;
   

  for(unsigned int j=0;j<2*n;j++) {
    drinv-=MinvDt[0][j]*delta[j];
    dphi0-=MinvDt[1][j]*delta[j];
    dt-=MinvDt[2][j]*delta[j];
    dz0-=MinvDt[3][j]*delta[j];
    dd0-=MinvDt[4][j]*delta[j];
    
    drinv_cov+=D[0][j]*delta[j];
    dphi0_cov+=D[1][j]*delta[j];
    dt_cov+=D[2][j]*delta[j];
    dz0_cov+=D[3][j]*delta[j];
    dd0_cov+=D[4][j]*delta[j];
  }
    
  double vpar[5];

  for(unsigned ipar=0;ipar<5;ipar++){
    vpar[ipar]=0.0;
    for(unsigned int j=0;j<2*n;j++) {
      vpar[ipar]+=MinvDt[ipar][j]*MinvDt[ipar][j];
    }
  }

  sigmarinv=sqrt(vpar[0]);
  sigmaphi0=sqrt(vpar[1]);
  sigmat=sqrt(vpar[2]);
  sigmaz0=sqrt(vpar[3]);
  sigmad0=sqrt(vpar[4]);

  double deltaChisq=drinv*drinv_cov+
    dphi0*dphi0_cov+
    dt*dt_cov+
    dz0*dz0_cov+
    dd0*dd0_cov;
  
  rinvfit = rinv+drinv;
  phi0fit = phi0+dphi0;
  tfit    = t+dt;
  z0fit   = z0+dz0;
  d0fit   = d0+dd0; 
   
  chisqfit = (chisq + deltaChisq);
  
  return;
}


//****************************************************************************
void L1PixelTrackFit::PrintProperties(void)
//****************************************************************************
{
  
  Table info("Bz | Var A | Var B | Var C", "Text");
  info.AddRowColumn(0, auxTools.ToString( mMagneticFieldStrength ) );
 // double rinvfit;
 // double phi0fit;
 // double d0fit;
 // double tfit;
 // double z0fit;
 // double chisqfit;
 // double sigmarinv;
 // double sigmaphi0;
 // double sigmad0;
 // double sigmat;
 // double sigmaz0;
  info.Print();

  return;
}


//****************************************************************************
void L1PixelTrackFit::SaveHitPermutationCandidate(int i0temp,
						  int i1temp,
						  int i2temp,
						  int i3temp,
						  double chiSq,
						  double chiSqDof,
						  double rinv,
						  double phi, 
						  double d0,
						  double t,
						  double z0,
						  string remarks)
//****************************************************************************
{
  
  int nHits = 0;
  if (i0temp > -1) nHits++;
  if (i1temp > -1) nHits++;
  if (i2temp > -1) nHits++;
  if (i3temp > -1) nHits++;

  cand_i0.push_back( i0temp );
  cand_i1.push_back( i1temp );
  cand_i2.push_back( i2temp );
  cand_i3.push_back( i3temp );
  cand_nHits.push_back(nHits);
  cand_chiSq.push_back( chiSq );
  cand_chiSqDof.push_back( chiSqDof );
  cand_rinv.push_back( rinv );
  cand_phi.push_back(  phi );
  cand_d0.push_back( d0 );
  cand_t.push_back(  t );
  cand_z0.push_back( z0 );
  cand_remarks.push_back( remarks );

 return;
}


//****************************************************************************
void L1PixelTrackFit::PrintHitPermutationCandidate(bool bSort)
//****************************************************************************
{

  if (bSort) SortHitPermutationCandidates();

  Table candInfo("Permutation | i0 | i1 | i2 | i3 | Fit Hits | chi2 | chi2/dof | fabs(1-chi2/dof) | rinv | phi0 | d0 | t | z0 | Remarks", "Text");
  int permutations = (int) cand_i0.size();
  
  // For-loop: All pixel hit permutations attempted for fitting tracks
  for(int i = 0; i < permutations; i++){

    int i0temp     = cand_i0.at(i);
    int i1temp     = cand_i1.at(i);
    int i2temp     = cand_i2.at(i);
    int i3temp     = cand_i3.at(i);
    int nHits      = cand_nHits.at(i);
    int nCandHits  = candidatePixHits.size();
    double chiSq   = cand_chiSq.at(i);
    double chiSqDof= cand_chiSqDof.at(i);
    double rinv    = cand_rinv.at(i);
    double phi     = cand_phi.at(i);
    double d0      = cand_d0.at(i);
    double t       = cand_t.at(i);
    double z0      = cand_z0.at(i);
    string remarks = cand_remarks.at(i);

    candInfo.AddRowColumn(i, auxTools.ToString( i+1 ) + " (" + auxTools.ToString(permutations) + ")" );
    candInfo.AddRowColumn(i, auxTools.ToString( i0temp ) );
    candInfo.AddRowColumn(i, auxTools.ToString( i1temp ) );
    candInfo.AddRowColumn(i, auxTools.ToString( i2temp ) );
    candInfo.AddRowColumn(i, auxTools.ToString( i3temp ) );   
    candInfo.AddRowColumn(i, auxTools.ToString( nHits ) + " (" + auxTools.ToString( nCandHits ) + ")" );
    candInfo.AddRowColumn(i, auxTools.ToString( chiSq ) );
    candInfo.AddRowColumn(i, auxTools.ToString( chiSqDof ) ); // Note: chi2=chi2/dof for 3 hits since (3*2)-5=1
    candInfo.AddRowColumn(i, auxTools.ToString( fabs(1.0-chiSqDof) ) ); // Note: chi2=chi2/dof for 3 hits since (3*2)-5=1
    candInfo.AddRowColumn(i, auxTools.ToString( rinv ) );
    candInfo.AddRowColumn(i, auxTools.ToString( phi ) );
    candInfo.AddRowColumn(i, auxTools.ToString( d0 ) );
    candInfo.AddRowColumn(i, auxTools.ToString( t ) );
    candInfo.AddRowColumn(i, auxTools.ToString( z0 ) );
    candInfo.AddRowColumn(i, remarks );
        
  }
  
  candInfo.Print();
  ClearCandInfo();
  
  return;
}


//****************************************************************************
void L1PixelTrackFit::SortHitPermutationCandidates(void)
//****************************************************************************
{

  // Do nothing if this function has already been called
  if (cand_sorted) return;
    
  // Build a vector consisting of the vector indices
  int permutations = (int) cand_i0.size(); 
  vector< int > cand_indices(permutations);
  for(int i = 0; i < permutations; i++){ cand_indices.at(i) = i; }

  // Sort indices according to chiSqDof (descending)
  // std::sort( std::begin(cand_indices), std::end(cand_indices), SortDescendingTarget( cand_chiSqDof ) );
  std::sort( std::begin(cand_indices), std::end(cand_indices), SortClosestToOneTarget( cand_chiSqDof ) );

  // Create temporary containers for sorted values (to be copied later to permanent containers)
  vector< int >    cand_i0_tmp;
  vector< int >    cand_i1_tmp;
  vector< int >    cand_i2_tmp;
  vector< int >    cand_i3_tmp;
  vector< int >    cand_nHits_tmp;
  vector< double > cand_chiSq_tmp;
  vector< double > cand_chiSqDof_tmp;
  vector< double > cand_rinv_tmp;
  vector< double > cand_phi_tmp;
  vector< double > cand_d0_tmp;
  vector< double > cand_t_tmp;
  vector< double > cand_z0_tmp;
  vector< string > cand_remarks_tmp;
  // For-loop: All pixel hit permutations
  for(int i = 0; i < permutations; i++){
    int j = cand_indices.at(i);
    cand_i0_tmp.push_back( cand_i0.at(j) );
    cand_i1_tmp.push_back( cand_i1.at(j) );
    cand_i2_tmp.push_back( cand_i2.at(j) );
    cand_i3_tmp.push_back( cand_i3.at(j) );
    cand_nHits_tmp.push_back( cand_nHits.at(j) );
    cand_chiSq_tmp.push_back( cand_chiSq.at(j) );
    cand_chiSqDof_tmp.push_back( cand_chiSqDof.at(j) );
    cand_rinv_tmp.push_back( cand_rinv.at(j) );
    cand_phi_tmp.push_back( cand_phi.at(j) );
    cand_d0_tmp.push_back( cand_d0.at(j) );
    cand_t_tmp.push_back( cand_t.at(j) );
    cand_z0_tmp.push_back( cand_z0.at(j) );
    cand_remarks_tmp.push_back( cand_remarks.at(j) );
  }

  // Now copied the customly sorted containers to their permanent containers
  cand_i0       = cand_i0_tmp;
  cand_i1       = cand_i1_tmp;
  cand_i2       = cand_i2_tmp;
  cand_i3       = cand_i3_tmp;
  cand_nHits    = cand_nHits_tmp;
  cand_chiSq    = cand_chiSq_tmp;
  cand_chiSqDof = cand_chiSqDof_tmp;
  cand_rinv     = cand_rinv_tmp;
  cand_phi      = cand_phi_tmp;
  cand_d0       = cand_d0_tmp;
  cand_t        = cand_t_tmp;
  cand_z0       = cand_z0_tmp;
  cand_remarks  = cand_remarks_tmp;
  cand_sorted   = true;
    
  return;
}


//****************************************************************************
void L1PixelTrackFit::ClearCandInfo(void)
//****************************************************************************
{

  cand_i0.clear();
  cand_i1.clear();
  cand_i2.clear();
  cand_i3.clear();
  cand_nHits.clear();
  cand_chiSq.clear();
  cand_chiSqDof.clear();
  cand_rinv.clear();
  cand_phi.clear(); 
  cand_d0.clear();
  cand_t.clear();
  cand_z0.clear();
  cand_remarks.clear();

  return;

}
#endif
