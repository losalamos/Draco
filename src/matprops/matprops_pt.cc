//----------------------------------*-C++-*----------------------------------//
// matprops_pt.cc
// Randy M. Roberts
// Thu May 28 14:56:12 1998
//---------------------------------------------------------------------------//
// @> Template instantiations for the matprops library.
//---------------------------------------------------------------------------//

#include "radphys/RadiationPhysics.t.cc"

typedef double T1;

template void XTM::RadiationPhysics::getElectIonCoupling<T1>(const T1 &,
							     const T1 &,
							     const T1 &,
							     double,
							     T1 &) const;

template void XTM::RadiationPhysics::getElectronConductionCoeff<T1>(const T1 &,
								    const T1 &,
								    const T1 &,
								    double,
								    T1 &) const;

template void XTM::RadiationPhysics::getIonConductionCoeff<T1>(const T1 &,
							       const T1 &,
							       const T1 &,
							       double,
							       T1 &) const;

//---------------------------------------------------------------------------//
//                              end of matprops_pt.cc
//---------------------------------------------------------------------------//
