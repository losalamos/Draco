#include "3T/testP13T/testP13T.hh"

#include "3T/testP13T/MeshTypeStub.hh"
#include "3T/testP13T/MaterialPropertiesStub.hh"
#include "3T/testP13T/DiffusionSolverStub.hh"
#include "3T/P13T.hh"
#include "3T/RadiationPhysics.hh"
#include "3T/P13TOptions.hh"
#include "3T/Units.hh"

#include "3T/testP13T/MaterialPropertiesStub.cc"
#include "3T/testP13T/DiffusionSolverStub.cc"
#include "3T/P13T.cc"
#include "3T/RadiationPhysics.cc"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

typedef MeshTypeStub MT;
typedef MaterialPropertiesStub<MT> MP;
typedef DiffusionSolverStub<MT> DS;

template class MaterialPropertiesStub<MT>;
template class DiffusionSolverStub<MT>;
template class P13T<MT,MP,DS>;
template void RadiationPhysics::getPlanck(const MT::ccsf &TElectron,
					  MT::ccsf &planckian) const;
template void RadiationPhysics::
    getPlanckTemperatureDerivative(const MT::ccsf &TElectron,
				   MT::ccsf &dplanckdT) const;

typedef MeshTypeStub::ccsf FT;
typedef int GT;

template void MaterialPropertiesStub<MT>::
getSigmaTotal(const MaterialStateField &matState,
	      GT group, FT &sigmaTotal) const;

template void MaterialPropertiesStub<MT>::
getSigmaAbsorption(const MaterialStateField &matState,
		   GT group, FT &sigmaAbsorption) const;

template void MaterialPropertiesStub<MT>::
getSigmaEmission(const MaterialStateField &matState,
		 GT group, FT &sigmaEmission) const;

template void MaterialPropertiesStub<MT>::
getElectronIonCoupling(const MaterialStateField &matState,
		       FT &gamma) const;

template void MaterialPropertiesStub<MT>::
getElectronTemperature(const MaterialStateField &matState,
		       FT &TElectron) const;

template void MaterialPropertiesStub<MT>::
getIonTemperature(const MaterialStateField &matState,
		  FT &TIon) const;

template void MaterialPropertiesStub<MT>::
getElectronConductionCoeff(const MaterialStateField &matState,
			   FT &kappaElectron) const;

template void MaterialPropertiesStub<MT>::
getIonConductionCoeff(const MaterialStateField &matState,
		      FT &kappaIon) const;

template void MaterialPropertiesStub<MT>::
getElectronSpecificHeat(const MaterialStateField &matState,
			FT &CvElectron) const;

template void MaterialPropertiesStub<MT>::
getIonSpecificHeat(const MaterialStateField &matState,
		   FT &CvIon) const;

END_NS_XTM  // namespace XTM
