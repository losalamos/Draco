//----------------------------------*-C++-*----------------------------------//
// radphys_pt.cc
// Randy M. Roberts
// Tue Nov 10 09:18:00 1998
//---------------------------------------------------------------------------//
// @> Template instantiations for the radphyslibrary.
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "RadiationPhysics.t.hh"

namespace rtt_radphys
{

 typedef double T1;

 template void RadiationPhysics::getElectIonCoupling<T1>(const T1 &,
							 const T1 &,
							 const T1 &,
							 double,
							 T1 &) const;

 template void RadiationPhysics::getElectronConductionCoeff<T1>(const T1 &,
								const T1 &,
								const T1 &,
								double,
								T1 &) const;

 template void RadiationPhysics::getIonConductionCoeff<T1>(const T1 &,
							   const T1 &,
							   const T1 &,
							   double,
							   T1 &) const;

 template void RadiationPhysics::getPlanck(const double &TElectron,
					   double &planckian) const;

} // end namespace rtt_radphys

//---------------------------------------------------------------------------//
//                              end of radphys_pt.cc
//---------------------------------------------------------------------------//
