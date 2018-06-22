#ifndef TrackingParticle_cxx
#define TrackingParticle_cxx

// User
#include "../interface/TrackingParticle.h"

// ****************************************************************************
TrackingParticle::TrackingParticle()
// ****************************************************************************
{

  _InitVars();

}

// ****************************************************************************
void TrackingParticle::_InitVars(void)
// ****************************************************************************
{


  // Variable declaration
  theIndex = 0.0;
  theMomentum.SetXYZ(0,0,0);
  theD0propagated = 0.0;
  theZ0propagated = 0.0;
  theCharge = 0.0;
  theQ = "N/A";
  theDxy = 0.0;
  theD0produced = 0.0;
  theZ0produced = 0.0;
  theD0Sign = 0.0;
  theD0Phi = 0.0;
  theNMatch = 0;
  theNStubs = 0;
  theTTTrackPt     = 0.0;
  theTTTrackEta    = 0.0;
  theTTTrackPhi    = 0.0;
  theTTTrackZ0     = 0.0;
  theTTTrackD0     = 0.0;
  theTTTrackChi2   = 0.0;
  theTTTrackNStubs = 0;
  theEventId = 0;
  
  return;
}

// ****************************************************************************
TrackingParticle::~TrackingParticle()
// ****************************************************************************
{

}


// ****************************************************************************
TrackingParticle::TrackingParticle(unsigned short index,
				   TVector3 momentum,
				   float d0_propagated,
				   float z0_propagated,
				   float d0_produced,
				   float z0_produced,
				   float dxy,
				   int charge,
				   int pdgId,
				   int nMatch,
				   int nStubs,
				   float ttTrackPt,
				   float ttTrackEta,
				   float ttTrackPhi,
				   float ttTrackZ0,
				   float ttTrackD0,
				   float ttTrackChi2,
				   int ttTrackNstubs,
				   int eventId)
// ****************************************************************************
{

  _InitVars();
  theIndex        = index;
  theMomentum     = momentum;
  theD0propagated = d0_propagated;
  theZ0propagated = z0_propagated;
  theCharge       = charge;
  theQ            = _getQ();
  thePdgId        = pdgId;
  theNMatch       = nMatch;
  theNStubs       = nStubs;
  theDxy          = dxy; //_getDxy();
  theD0produced  = d0_produced; //_getD0();
  theZ0produced  = z0_produced;
  //theD0Sign       = _getD0Sign();
  //theD0Phi        = _getD0Phi();
  theNMatch       = nMatch;
  theTTTrackPt    = ttTrackPt;
  theTTTrackEta   = ttTrackEta;
  theTTTrackPhi   = ttTrackPhi;
  theTTTrackZ0    = ttTrackZ0;
  theTTTrackD0    = ttTrackD0;
  theTTTrackChi2  = ttTrackChi2;
  theTTTrackNStubs= ttTrackNstubs;
  theEventId      = eventId;

  if (0) PrintProperties();
}

/*
// ****************************************************************************
float TrackingParticle::_getD0Phi(void)
// ****************************************************************************
{

  float d0phi = atan2(thePOCA.Y(), thePOCA.X()); // range(-pi, +pi)
  return d0phi;
}


// ****************************************************************************
float TrackingParticle::_getD0Sign(void)
// ****************************************************************************
{

  // The track d0 sign has the same sign as the curvature (rho)
  // [and hence charge (Q) since Sgn(Q) = Sgn(rho)] if the z-axis is OUTSIDE
  // the circle in the x-y plane formed by the track helix projection.
  // In other words, if the track is going clockwise in the x-y plane it has a
  // positive d0 sign, otherwise negative. The latter is irrespective of charge.
  // One way to determine the track direction (clockwise or anticlockwise) is to
  // first take the cross product of the POCA vector with the instantaneous
  // momentum at the POCA and then project along the unit vector of the viewing direction.
  // If you are looking AT the x-y plane the we are looking along the negative z-axis.
  // Therefore we need to project along (0, 0, -1)
  // doublem tp_RInv = t->TP_RInv->at(iTP);

  // Q. How to determine if a vector is moving clockwise or anti-clockwise to another
  // A. The question only makes sense in 3-D since the notion only makes sense
  // once you have established a viewing direction. If V is the viewing direction vector
  // and A and B are the individual vectors, then what you look at is the sign of
  // (A x B) o V
  // where A x B is the vector cross product and "o" is the dot product.
  // If your direction vector is looking down on the Cartesian plane then V = (0, 0, -1) and
  // then a positive sign indicates a clockwise move. 
  // To show this, let A = (1, 0, 0) and let B = (0, 1, 0); that is, A is along the x-axis and
  // vector B is along the y-axis. Then A x B = (0, 0, 1) and it is pointing up along the Z axis.
  // The dot product (0, 0, 1) o (0, 0, -1) = -1 and, indeed, from the X axis to the Y axis is
  // counter-clockwise when viewed from the +Z side of the plane.

  // Check the dot product sign: vec{a} dot vec{b}  = ax*bx + ay*by + az*bz
  // a = thePOCA vector
  // b = theMomentum unit vector at the POCA
  // TVector3 cross_product = thePOCA.Cross( theMomentum.Unit() );
  ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag> cross_product = thePOCA.Cross( theMomentum.Unit() );

  // For looking down on the Cartesian plane the viewing direction vector (0, 0, -1) is (0, 0, -1)
  TVector3 v_z(0, 0, -1);
  float dot_product = cross_product.Dot(v_z);  
  
  // If dot product is positive then track is clockwise. Otherwise anticlockwise
  if (dot_product >= 0.0) return +1.0;
  else return -1.0;

}
*/

// ****************************************************************************
string TrackingParticle::_getQ(void)
// ****************************************************************************
{

  string theQ = "N/A";
  if (theCharge > 0) theQ = "+";
  else if (theCharge < 0) theQ = "-";
  else
    {
      theQ = "0";
      // std::cout << "E R R O R ! TrackingParticle::getQ(...) - Invalid value for theCharge \"" << theCharge << "\". EXIT" << std::endl;
      // exit(1);
    }

  return theQ;
}
  

// ****************************************************************************
void TrackingParticle::PrintProperties(bool bPrintTitleRow)
// ****************************************************************************
{
  
  Table info("Index | Pt | Eta | Phi | PdgId | Q | z0 propagated | d0 propagated | z0 produced | d0 produced | dxy | NMatch | NStubs | TTTrackPt | TTTrackEta | TTTrackPhi | TTTrackZ0 | TTTrackD0 | TTTrackChi2 | TTTrackNStubs | Event-Id", "Text");

  info.AddRowColumn(0, auxTools.ToString( theIndex) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Perp(), 3) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Eta() , 3) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Phi() , 3) );
  info.AddRowColumn(0, auxTools.ToString( getPdgId()) );
  info.AddRowColumn(0, getQ());
  info.AddRowColumn(0, auxTools.ToString( getZ0propagated(), 3)    );
  info.AddRowColumn(0, auxTools.ToString( getD0propagated(), 3)    );
  info.AddRowColumn(0, auxTools.ToString( getZ0produced(), 3)    );
  info.AddRowColumn(0, auxTools.ToString( getD0produced(), 3)    );
  info.AddRowColumn(0, auxTools.ToString( theDxy , 3)    );
  //info.AddRowColumn(0, auxTools.ToString( getD0Sign() )  );
  //info.AddRowColumn(0, auxTools.ToString( getD0Phi(), 3) );
  info.AddRowColumn(0, auxTools.ToString( theNMatch)       );
  info.AddRowColumn(0, auxTools.ToString( theNStubs)       );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackPt)    );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackEta)   );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackPhi)   );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackZ0)    );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackD0)    );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackChi2)  );
  info.AddRowColumn(0, auxTools.ToString( theTTTrackNStubs));
  info.AddRowColumn(0, auxTools.ToString( theEventId)      );
  info.Print(bPrintTitleRow);

  return;
}

// ****************************************************************************
TLorentzVector TrackingParticle::p4(float mass)
// ****************************************************************************
{  
  TLorentzVector p4;
  p4.SetPtEtaPhiM(getPt(), getEta(), getPhi(), mass);
  return p4;
}



#endif
