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
getSigmaTotal(const MaterialStateField<FT> &matState,
	      GT group, FT &sigmaTotal) const
{
    sigmaTotal = matVals[matState.val].sigmaTotal;
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaAbsorption(const MaterialStateField<FT> &matState,
		   GT group, FT &sigmaAbsorption) const
{
    sigmaAbsorption = matVals[matState.val].sigmaAbsorption;
}

template<class MT>
template<class FT, class GT>
void MaterialPropertiesStub<MT>::
getSigmaEmission(const MaterialStateField<FT> &matState,
		 GT group, FT &sigmaEmission) const
{
    sigmaEmission = matVals[matState.val].sigmaEmission;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronIonCoupling(const MaterialStateField<FT> &matState,
		       FT &gamma) const
{
    gamma = matVals[matState.val].electronIonCoupling;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronTemperature(const MaterialStateField<FT> &matState,
		       FT &TElectron_) const
{
    TElectron_ = TElectron;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonTemperature(const MaterialStateField<FT> &matState,
		  FT &TIon_) const
{
    TIon_ = TIon;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronConductionCoeff(const MaterialStateField<FT> &matState,
			   FT &kappaElectron) const
{
    kappaElectron = matVals[matState.val].electronConductionCoeff;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonConductionCoeff(const MaterialStateField<FT> &matState,
		      FT &kappaIon) const
{
    kappaIon = matVals[matState.val].ionConductionCoeff;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getElectronSpecificHeat(const MaterialStateField<FT> &matState,
			FT &CvElectron) const
{
    CvElectron = matVals[matState.val].electronSpecificHeat;
}

template<class MT>
template<class FT>
void MaterialPropertiesStub<MT>::
getIonSpecificHeat(const MaterialStateField<FT> &matState,
		   FT &CvIon) const
{
    CvIon = matVals[matState.val].ionSpecificHeat;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of MaterialPropertiesStub.cc
//---------------------------------------------------------------------------//
