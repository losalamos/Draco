//----------------------------------*-C++-*----------------------------------//
// testRadPhys_pt.cc
// Randy M. Roberts
// Thu May 14 15:00:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "testRadPhys_pt.hh"

#include <vector>
typedef std::vector<double> VD;

#include "../RadiationPhysics.t.hh"
using namespace XTM;

template 
void RadiationPhysics::getPlanck(const VD &TElectron, VD &planckian) const;

template 
void RadiationPhysics::getPlanckTemperatureDerivative(const VD &TElectron,
						      VD &planckian) const;

template
void RadiationPhysics::getElectIonCoupling(const double &, const double &,
					   const double &, double,
					   double &) const;

VD func()
{
    VD a(3, 42.);
    const VD b(3, 23.);

    Units units(Units::getAstroPhysUnits());	       

    VD c = units.ConvertDensity(a) / b;

    return c;
}

#if 0
template
void RadiationPhysics::getElectIonCoulombLog(const VD &, const VD &,
					     const VD &, double,
					     VD &) const;
#endif

template
void RadiationPhysics::getElectIonCoupling(const VD &, const VD &,
					   const VD &, double,
					   VD &) const;

//---------------------------------------------------------------------------//
//                              end of testRadPhys_pt.cc
//---------------------------------------------------------------------------//
