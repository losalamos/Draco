#include "3T/testP13T/testP13T.hh"

#include "3T/testP13T/MeshTypeStub.hh"
#include "matprops/InterpedMaterialProps.hh"
#include "3T/testP13T/DiffusionSolverStub.hh"
#include "3T/P13T.hh"
#include "radphys/RadiationPhysics.hh"
#include "3T/P13TOptions.hh"
#include "units/Units.hh"

#include "3T/testP13T/DiffusionSolverStub.cc"
#include "3T/P13T.cc"
#include "radphys/RadiationPhysics.cc"
#include "matprops/InterpedMaterialProps.cc"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

typedef MeshTypeStub MT;
typedef InterpedMaterialProps MP;
typedef DiffusionSolverStub<MT> DS;

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

typedef MeshTypeStub::scalar T1;
typedef MeshTypeStub::scalar T2;

template
MP::MaterialStateField<T1>
MP::getMaterialState<T1, T2>(const T1 &, const T1 &, const T1 &,
			     const T2 &) const;

END_NS_XTM  // namespace XTM
