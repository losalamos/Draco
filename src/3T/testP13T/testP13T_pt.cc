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

typedef MeshTypeStub::ccsf ccsf;
typedef MeshTypeStub::fcdsf fcdsf;
typedef int GT;

template void MaterialPropertiesStub<MT>::
getSigmaTotal(const MaterialStateField<fcdsf> &matState,
	      GT group, fcdsf &sigmaTotal) const;

template void MaterialPropertiesStub<MT>::
getSigmaAbsorption(const MaterialStateField<ccsf> &matState,
		   GT group, ccsf &sigmaAbsorption) const;

template void MaterialPropertiesStub<MT>::
getSigmaEmission(const MaterialStateField<ccsf> &matState,
		 GT group, ccsf &sigmaEmission) const;

template void MaterialPropertiesStub<MT>::
getElectronIonCoupling(const MaterialStateField<ccsf> &matState,
		       ccsf &gamma) const;

template void MaterialPropertiesStub<MT>::
getElectronTemperature(const MaterialStateField<ccsf> &matState,
		       ccsf &TElectron) const;

template void MaterialPropertiesStub<MT>::
getIonTemperature(const MaterialStateField<ccsf> &matState,
		  ccsf &TIon) const;

template void MaterialPropertiesStub<MT>::
getElectronConductionCoeff(const MaterialStateField<fcdsf> &matState,
			   fcdsf &kappaElectron) const;

template void MaterialPropertiesStub<MT>::
getIonConductionCoeff(const MaterialStateField<fcdsf> &matState,
		      fcdsf &kappaIon) const;

template void MaterialPropertiesStub<MT>::
getElectronSpecificHeat(const MaterialStateField<ccsf> &matState,
			ccsf &CvElectron) const;

template void MaterialPropertiesStub<MT>::
getIonSpecificHeat(const MaterialStateField<ccsf> &matState,
		   ccsf &CvIon) const;

END_NS_XTM  // namespace XTM
