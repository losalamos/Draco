//----------------------------------*-C++-*----------------------------------//
// P13T.cc
// Randy M. Roberts
// Wed Mar 11 11:17:54 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/P13T.hh"
#include "3T/RadiationPhysics.hh"
#include "3T/PhysicalConstants.hh"
#include <cmath>

#include <iostream>
using std::cerr;
using std::endl;

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
// CREATORS

//---------------------------------------------------------------------------//
// P13T:
//     Constructor
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>::P13T(const P13TOptions &options_,
		       const SP<MP> &spProp_,
		       const SP<DS> &spDiffSolver_)
    : options(options_), spProp(spProp_), spDiffSolver(spDiffSolver_)
{
    // empty
}

//---------------------------------------------------------------------------//
// P13T:
//     Copy constructor
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>::P13T(const P13T<MT, MP, DS> &rhs)
    : options(rhs.options), spProp(rhs.spProp), spDiffSolver(rhs.spDiffSolver)
{
    // empty
}

//---------------------------------------------------------------------------//
// ~P13T:
//     Destructor
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>::~P13T()
{
    // empty
}



// MANIPULATORS

//---------------------------------------------------------------------------//
// operator=:
//   Assignment operator
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>& P13T<MT, MP, DS>::operator=(const P13T &rhs)
{
    options = rhs.options;
    spProp = rhs.spProp;
    spDiffSolver = rhs.spDiffSolver;
    return *this;
}


//---------------------------------------------------------------------------//
// setMaterialProperties:
//     Set the contained material property to the new value.
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::
setMaterialProperties(const SP<MaterialProperties> &spProp_)
{
    spProp = spProp_;
}

//---------------------------------------------------------------------------//
// setOptions:
//     Set the contained options object to the new value.
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::setOptions(const P13TOptions options_)
{
    options = options_;
}

//---------------------------------------------------------------------------//
// setDiffSolver:
//     Set the diffusion solver to be used during the solves to the new
//     value.
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::setDiffSolver(const SP<DiffusionSolver> &spDiffSolver_)
{
    spDiffSolver = spDiffSolver_;
}


// ACCESSORS

//------------------------------------------------------------------------//
// print:
//     Print itself (for debug mostly)
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
std::ostream &P13T<MT, MP, DS>::print(std::ostream &os) const
{
    os << "(P13T::this: " << (void *)this
       << " spProp: " << *spProp
       << " spDiffSolver: " << *spDiffSolver
       << ")";
    return os;
}

//---------------------------------------------------------------------------//
// initializeRadiationState:
//     Initialize the radiation field to Planckian
//     based on material electron temperatures.
//---------------------------------------------------------------------------//
    
template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::
initializeRadiationState(const MaterialStateField &matState,
			 RadiationStateField &resultsStateField) const
{
    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    ccsf TElectron(spMesh);
    spProp->getElectronTemperature(matState, TElectron);
    
    // Set the radiation physics to the given units.

    const RadiationPhysics radPhys(spProp->getUnits());

    // The radiation field is set to the 4pi*planckian with no current flow.

    radPhys.getPlanck(TElectron, resultsStateField.phi);

    using XTM::PhysicalConstants::pi;
    
    resultsStateField.phi *= 4.0*pi;
    resultsStateField.F = 0.0;
}

//---------------------------------------------------------------------------//
// solve:
//     Solve for the new radiation field, the electron/ion energy depositions,
//     and the momentom deposition.
//
//     The P13TOptions object (P13T state variable "options")
//     determines whether this solve is with or without the
//     electron/ion conduction equations.
//---------------------------------------------------------------------------//
    
template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::solve(double dt,
			     const MaterialStateField &matState,
			     const RadiationStateField &prevStateField,
			     const ccsf QRad,
			     const ccsf QElectron,
			     const ccsf QIon,
			     const bsbf boundary,
			     RadiationStateField &resultsStateField,
			     ccsf &electronEnergyDeposition,
			     ccsf &ionEnergyDeposition,
			     ncvf &momentumDeposition,
			     ccsf &Tnp1Electron,
			     ccsf &Tnp1Ion) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);

    // This is a one-group problem

    int groupNo = 1;
    int numGroups = 1;

    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Let's save the temperatures at time t^n.
    
    ccsf TnElectron(spMesh);
    ccsf TnIon(spMesh);
    spProp->getElectronTemperature(matState, TnElectron);
    spProp->getIonTemperature(matState, TnIon);

    // If there is conduction we will need to differentiate
    // between the temperatures at time t^n and t^n+1/2.
    
    ccsf Tnp12Electron = TnElectron;
    ccsf Tnp12Ion = TnIon;

    // If they want ion conduction, we do it here.
    
    if (options.wantIonConduction())
    {
	// Calculate the conduction coefficient
	
	fcdsf kappaIon(spMesh);
	spProp->getIonConductionCoeff(matState, kappaIon);
	kappaIon *= dt;

	// Calculate the removal coefficient.
	
	ccsf removalCoeff(spMesh);
	spProp->getIonSpecificHeat(matState, removalCoeff);

	// Calculate the source term.
	ccsf source = removalCoeff * TnIon;

	// Do the diffusion solve for the temperature at n+1/2.
	
	spDiffSolver->solve(kappaIon, removalCoeff, source,
			    boundary, Tnp12Ion);
    }

    
    // If they want electron conduction, we do it here.
    
    if (options.wantElectronConduction())
    {
	// Calculate the conduction coefficient
	
	fcdsf kappaElectron(spMesh);
	spProp->getElectronConductionCoeff(matState, kappaElectron);
	kappaElectron *= dt;
	
	// Calculate the removal coefficient.
	
	ccsf removalCoeff(spMesh);
	spProp->getElectronSpecificHeat(matState, removalCoeff);

	// Calculate the source term.
	ccsf source = removalCoeff * TnElectron;

	// Do the diffusion solve for the temperature at n+1/2.
	
	spDiffSolver->solve(kappaElectron, removalCoeff, source,
			    boundary, Tnp12Electron);
    }

    // Calculate the new radiation state due to the radiation P1,
    // electron and ion equations ***without*** the conduction equations.
    
    calcNewRadState(dt, groupNo, matState, prevStateField,
		    QRad, QElectron, QIon,
		    Tnp12Electron, Tnp12Ion,
		    boundary,
		    resultsStateField);

    // Calculate the delta electron temperature from the new radiation
    // state, (Te^n+1 - Te^n+1/2).
    
    ccsf deltaTElectron(spMesh);
    calcDeltaTElectron(dt, numGroups, matState, prevStateField, QElectron, QIon,
		       Tnp12Electron, Tnp12Ion,
		       resultsStateField, deltaTElectron);
    
    // Calculate the delta ion temperature from the delta electron
    // temperature, (Ti^n+1 - Ti^n+1/2).
    
    ccsf deltaTIon(spMesh);
    calcDeltaTIon(dt, matState, prevStateField, QIon, Tnp12Electron, Tnp12Ion,
		  deltaTElectron, deltaTIon);

    // Now include the contribution from conduction.
    
    if (options.wantIonConduction())
    {
	deltaTIon += Tnp12Ion - TnIon;
    }

    if (options.wantElectronConduction())
    {
	deltaTElectron += Tnp12Electron - TnElectron;
    }

    // Calculate electron and ion energy deposition

    ccsf Cv(spMesh);

    spProp->getElectronSpecificHeat(matState, Cv);
    electronEnergyDeposition = Cv * deltaTElectron;

    Tnp1Electron = TnElectron + deltaTElectron;

    spProp->getIonSpecificHeat(matState, Cv);
    ionEnergyDeposition = Cv * deltaTIon;

    Tnp1Ion = TnIon + deltaTIon;
    
    // Calculate the momentum deposition.
    // For now we just zero it out.

    momentumDeposition = 0.0;
}

//---------------------------------------------------------------------------//
// clacNewRadState:
//     calculate the new radiation state using the previous state,
//     material properties, and sources.
//     This solves the coupled radiation, electron, and ion equations
//     ***without*** the conduction equations.
//---------------------------------------------------------------------------//


template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::
calcNewRadState(double dt,
		int groupNo,
		const MaterialStateField &matState,
		const RadiationStateField &prevStateField,
		const ccsf QRad,
		const ccsf QElectron,
		const ccsf QIon,
		const ccsf TElectron,
		const ccsf TIon,
		const bsbf boundary,
		RadiationStateField &resultsStateField) const
{
    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Calculate the coefficients needed by the diffusion solver.

    fcdsf D(spMesh);
    DiscFluxField Fprime(spMesh);
    ccsf sigmaAbsBar(spMesh);
    ccsf QRadBar(spMesh);
    
    calcP1Coeffs(dt, groupNo, matState, prevStateField, QRad, QElectron, QIon,
		 TElectron, TIon, D, Fprime, sigmaAbsBar, QRadBar);

    // Set up aliases for the long names.
    
    ccsf &phi = resultsStateField.phi;
    FluxField &F = resultsStateField.F;

    // Call the diffusion solver to solve for phi and F.
    // Since phi and F are aliased to the resultsStateField members,
    // this is all we have to do.
    
    spDiffSolver->solve(D, sigmaAbsBar, QRadBar, Fprime, boundary, phi, F);

    // ENSURE that F is indeed continuous, if we can!
}
    
//---------------------------------------------------------------------------//
// calcP1Coeffs:
//     Calculate the coefficients, e.g. diffusion and removal, and
//     source terms for solving the P1 equation.
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::
calcP1Coeffs(double dt,
	     int groupNo,
	     const MaterialStateField &matState,
	     const RadiationStateField &prevStateField,
	     const ccsf &QRad,
	     const ccsf &QElectron,
	     const ccsf &QIon,
	     const ccsf &TElectron,
	     const ccsf &TIon,
	     fcdsf &D,
	     DiscFluxField &Fprime,
	     ccsf &sigmaAbsBar,
	     ccsf &QRadBar) const
{
    // Set the radiation physics to the given units.
    
    const RadiationPhysics radPhys(spProp->getUnits());

    // set up some needed scalars, like tau

    double c = radPhys.getLightSpeed();
    double tau = 1.0/(c*dt);

    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Ask the material properties for sigma total.
    // It is the material properties responsibility to do
    // any averaging of temperatures, etc. to achieve the correct
    // resulting sigmaTotal.

    fcdsf sigmaTotal(spMesh);
    spProp->getSigmaTotal(matState, groupNo, sigmaTotal);

    //
    // We can now calculate the results
    //

    ccsf sigmaAbs(spMesh);
    spProp->getSigmaAbsorption(matState, groupNo, sigmaAbs);
    
    // Calculate the diffusion constant.
    
    D = (1.0/3.0) / (sigmaTotal + tau);

    // We need nu and QElecStar, we get CvStar for free.

    // Get sigma emission
    ccsf sigmaEmission(spMesh);
    spProp->getSigmaEmission(matState, groupNo, sigmaEmission);

    ccsf QElecStar(spMesh);
    ccsf CvStar(spMesh);
    ccsf nu(spMesh);
    
    calcStarredFields(dt, groupNo, matState, radPhys,
		      QElectron, QIon, TElectron, TIon,
		      sigmaEmission, QElecStar, CvStar, nu);

    // Calculate modified sigma absorption

    sigmaAbsBar = (1.0 - nu) * sigmaAbs + tau;

    // Calculated modified radiation source

    // We need the Planckian.

    ccsf planck(spMesh);
    radPhys.getPlanck(TElectron, planck);

    // We need our Planckian multiplied by 4pi

    using XTM::PhysicalConstants::pi;
    planck *= 4.0*pi;

    QRadBar = tau*prevStateField.phi + (1.0 - nu)*sigmaEmission*planck
	+ nu*QElecStar;

    // Calculate the "telegraph" term to the P1 equation.

    Fprime = tau*prevStateField.F / (sigmaTotal + tau);

}

//------------------------------------------------------------------------//
// calcStarredFields:
//    Calculate Qe*, Cv*, but not nu.
//    These are needed to calculate other coefficients
//    and delta temperatures.
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::calcStarredFields(double dt,
					 int groupNo,
					 const MaterialStateField &matState,
					 const RadiationPhysics &radPhys,
					 const ccsf &QElectron,
					 const ccsf &QIon,
					 const ccsf &TElectron,
					 const ccsf &TIon,
					 const ccsf &sigmaEmission,
					 ccsf &QElecStar,
					 ccsf &CvStar,
					 ccsf &nu) const
{

    // Calculate Qe* and Cv*.
    // We will then calculate nu ourself.
    
    calcStarredFields(dt, groupNo, matState, radPhys,
		      QElectron, QIon, TElectron, TIon,
		      sigmaEmission, QElecStar, CvStar);

    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Calculate the Planckian's temperature derivative.

    ccsf dPlanckdT(spMesh);
    radPhys.getPlanckTemperatureDerivative(TElectron, dPlanckdT);

    // We need our Planckian multiplied by 4pi

    using XTM::PhysicalConstants::pi;
    dPlanckdT *= 4.0*pi;
    
    // Calculate the "nu" used in the 3T modification of sigmaAbs
    
    nu = dt * sigmaEmission * dPlanckdT / (CvStar + dt * dPlanckdT);

}

//------------------------------------------------------------------------//
// calcStarredFields:
//    Calculate Qe*, Cv*, but not nu.
//    These are needed to calculate other coefficients
//    and delta temperatures.
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::calcStarredFields(double dt,
					 int groupNo,
					 const MaterialStateField &matState,
					 const RadiationPhysics &radPhys,
					 const ccsf &QElectron,
					 const ccsf &QIon,
					 const ccsf &TElectron,
					 const ccsf &TIon,
					 const ccsf &sigmaEmission,
					 ccsf &QElecStar,
					 ccsf &CvStar) const
{
    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Get the electron and ion heat capacities.
    
    ccsf CvElec(spMesh);
    spProp->getElectronSpecificHeat(matState, CvElec);

    ccsf CvIon(spMesh);
    spProp->getIonSpecificHeat(matState, CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    spProp->getElectronIonCoupling(matState, gamma);

    // tmpCoeff is a common term to two calculations.
    // Let's just do it once.
    
    const ccsf tmpCoeff = (gamma*dt) / (CvIon + gamma*dt);

    // CvStar is one of the results, as well as intermediate.
    
    CvStar = CvElec + CvIon * tmpCoeff;
    
    // Calculate QElecStar (Qe*).

    QElecStar = QElectron + (CvIon/dt * (TIon - TElectron) + QIon) * tmpCoeff;

}

//------------------------------------------------------------------------//
// calcDeltaTElectron:
//    Calculate the difference between T electron from timestep
//    n+1 to timestep n+1/2
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::calcDeltaTElectron(double dt,
					  int numGroups, 
					  const MaterialStateField &matState, 
					  const RadiationStateField &prevStateField,
					  const ccsf &QElectron, 
					  const ccsf &QIon,
					  const ccsf &TElectron,
					  const ccsf &TIon,
					  const RadiationStateField &resultsStateField, 
					  ccsf &deltaTelectron) const
{
    // Only one group.

    Require(numGroups == 1);
    int groupNo = 1;
    
    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Set the radiation physics to the given units.
    
    const RadiationPhysics radPhys(spProp->getUnits());

    ccsf sigmaEmission(spMesh);
    spProp->getSigmaEmission(matState, groupNo, sigmaEmission);

    // Calculate QElecStar and CvStar.

    ccsf QElecStar(spMesh);
    ccsf CvStar(spMesh);
    
    calcStarredFields(dt, groupNo, matState, radPhys,
		      QElectron, QIon, TElectron, TIon,
		      sigmaEmission, QElecStar, CvStar);

    ccsf sigmaAbs(spMesh);
    spProp->getSigmaAbsorption(matState, groupNo, sigmaAbs);

    ccsf planck(spMesh);
    radPhys.getPlanck(TElectron, planck);
    
    ccsf dPlanckdT(spMesh);
    radPhys.getPlanckTemperatureDerivative(TElectron, dPlanckdT);

    // The planckian and derivative must be multiplies by 4pi
    
    using XTM::PhysicalConstants::pi;
    planck *= 4.0*pi;
    dPlanckdT *= 4.0*pi;
    
    // Get shorthand for phi^n+1
    const ccsf &phi_np1 = resultsStateField.phi;

    // calculate delta T electron
    
    deltaTelectron = dt * (sigmaAbs*phi_np1 - sigmaEmission*planck + QElecStar)
	/ (CvStar + dt*sigmaEmission*dPlanckdT);
}

//------------------------------------------------------------------------//
// calcDeltaTIon:
//    Calculate the difference between T ion from timestep
//    n+1 to timestep n+1/2
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::calcDeltaTIon(double dt,
				     const MaterialStateField &matState, 
				     const RadiationStateField &prevStateField, 
				     const ccsf &QIon,
				     const ccsf &TElectron,
				     const ccsf &TIon,
				     const ccsf &deltaTelectron,
				     ccsf &deltaTIon) const
{
    // We need a smart pointer to a mesh

    const SP<MeshType> spMesh = spDiffSolver->getMesh();

    // Get the ion heat capacity.
    
    ccsf CvIon(spMesh);
    spProp->getIonSpecificHeat(matState, CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    spProp->getElectronIonCoupling(matState, gamma);

    deltaTIon = dt * (gamma*(TElectron - TIon + deltaTelectron) + QIon)
	/ (CvIon + dt*gamma);

}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of P13T.cc
//---------------------------------------------------------------------------//
