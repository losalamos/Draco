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
MaterialPropertiesStub<MT>::
MaterialPropertiesStub(const Units &units_,
		       const SP<MT> &spMesh_,
		       double TElectron_,
		       double TIon_,
		       const std::vector<MatVals> &matVals_)
    : units(units_), spMesh(spMesh_), TElectron(TElectron_), TIon(TIon_),
      matVals(matVals_)
{
    // empty
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaTotal(const MaterialStateField &matState,
	      GT group, FT &sigmaTotal) const
{
    sigmaTotal = matVals[matState].sigmaTotal;
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaAbsorption(const MaterialStateField &matState,
		   GT group, FT &sigmaAbsorption) const
{
    sigmaAbsorption = matVals[matState].sigmaAbsorption;
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaEmission(const MaterialStateField &matState,
		 GT group, FT &sigmaEmission) const
{
    sigmaEmission = matVals[matState].sigmaEmission;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronIonCoupling(const MaterialStateField &matState,
		       FT &gamma) const
{
    gamma = matVals[matState].electronIonCoupling;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronTemperature(const MaterialStateField &matState,
		       FT &TElectron_) const
{
    TElectron_ = TElectron;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonTemperature(const MaterialStateField &matState,
		  FT &TIon_) const
{
    TIon_ = TIon;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronConductionCoeff(const MaterialStateField &matState,
			   FT &kappaElectron) const
{
    kappaElectron = matVals[matState].electronConductionCoeff;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonConductionCoeff(const MaterialStateField &matState,
		      FT &kappaIon) const
{
    kappaIon = matVals[matState].ionConductionCoeff;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronSpecificHeat(const MaterialStateField &matState,
			FT &CvElectron) const
{
    CvElectron = matVals[matState].electronSpecificHeat;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonSpecificHeat(const MaterialStateField &matState,
		   FT &CvIon) const
{
    CvIon = matVals[matState].ionSpecificHeat;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of MaterialPropertiesStub.cc
//---------------------------------------------------------------------------//
