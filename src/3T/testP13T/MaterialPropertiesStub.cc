//----------------------------------*-C++-*----------------------------------//
// MaterialPropertiesStub.cc
// Randy M. Roberts
// Mon Mar 23 08:41:02 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/MaterialPropertiesStub.hh"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

#include "3T/Units.hh"

BEGIN_NS_XTM

template<class MT>
MaterialPropertiesStub<MT>::MaterialPropertiesStub(const Units &units_,
						   const SP<MT> &spMesh_)
    : units(units_), spMesh(spMesh_)
{
    // empty
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaTotal(const MaterialStateField &matState,
	      GT group, FT &sigmaTotal) const
{
    sigmaTotal = 2.0;
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaAbsorption(const MaterialStateField &matState,
		   GT group, FT &sigmaAbsorption) const
{
    sigmaAbsorption = 1.0;
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaEmission(const MaterialStateField &matState,
		 GT group, FT &sigmaEmission) const
{
    sigmaEmission = 1.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronIonCoupling(const MaterialStateField &matState,
		       FT &gamma) const
{
    gamma = 0.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronTemperature(const MaterialStateField &matState,
		       FT &TElectron) const
{
    TElectron = 2.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonTemperature(const MaterialStateField &matState,
		  FT &TIon) const
{
    TIon = 3.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronConductionCoeff(const MaterialStateField &matState,
			   FT &kappaElectron) const
{
    kappaElectron = 4.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonConductionCoeff(const MaterialStateField &matState,
		      FT &kappaIon) const
{
    kappaIon = 4.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronSpecificHeat(const MaterialStateField &matState,
			FT &CvElectron) const
{
    CvElectron = 1.0;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonSpecificHeat(const MaterialStateField &matState,
		   FT &CvIon) const
{
    CvIon = 2.0;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of MaterialPropertiesStub.cc
//---------------------------------------------------------------------------//
