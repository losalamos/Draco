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

// CREATORS

//---------------------------------------------------------------------------//
// P13T:
//     Constructor
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>::P13T(const P13TOptions &options_,
		       const SP<MaterialProperties> &spProp_,
		       const SP<DiffusionSolver> &spDiffSolver_)
    : options(options_), spProp(spProp_), spDiffSolver(spDiffSolver)
{
    // empty
}

//---------------------------------------------------------------------------//
// P13T:
//     Copy constructor
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>::P13T(const P13T &rhs)
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
P13T& P13T<MT, MP, DS>::operator=(const P13T &rhs)
{
    options = rhs.options;
    spProp = rhs.spProp;
    spDiffSolver = rhs.spDiffSolver;
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
    
    radPhysics.getPlanck(TElectron, resultsStateField.phi);

    using PhysicalConstants::pi;
    
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
			     ncvf &momentumDeposition) const
{
    // Require dt > 0, etc.
    //
    // REQUIRE(dt > 0.0)

   
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
    
    ccsf Tnp12Electron(spMesh) = TnElectron;
    ccsf Tnp12Ion(spMesh) = TnIon;

    // If they want ion conduction, we do it here.
    
    if (options.wantIonConduction())
    {
	// Calculate the conduction coefficient
	
	fcdsf kappaIon(spMesh);
	spProp->getIonConductionCoeff(matState, kappaIon);

	// Calculate the removal coefficient.
	
	ccsf removalCoeff(spMesh);
	spProp->getIonSpecificHeat(matState, removalCoeff);
	removalCoeff /= dt;

	// Calculate the source term.
	ccsf source = removalCoeff * TnIon;

	// Do the diffusion solve for the temperature at n+1/2.
	
	spDiffSolver->solve(kappaIon, removalCoeff, source,
			    Tnp12Ion);
    }

    
    // If they want electron conduction, we do it here.
    
    if (options.wantElectronConduction())
    {
	// Calculate the conduction coefficient
	
	fcdsf kappaElectron(spMesh);
	spProp->getElectronConductionCoeff(matState, kappaElectron);

	// Calculate the removal coefficient.
	
	ccsf removalCoeff(spMesh);
	spProp->getElectronSpecificHeat(matState, removalCoeff);
	removalCoeff /= dt;

	// Calculate the source term.
	ccsf source = removalCoeff * TnElectron;

	// Do the diffusion solve for the temperature at n+1/2.
	
	spDiffSolver->solve(kappaElectron, removalCoeff, source,
			    Tnp12Electron);
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
    calcDeltaElectron(dt, numGroups, matState, prevStateField, QElectron, QIon,
		      resultsStateField, deltaTelectron);
    
    // Calculate the delta ion temperature from the delta electron
    // temperature, (Ti^n+1 - Ti^n+1/2).
    
    ccsf deltaTIon(spMesh);
    calcDeltaTIon(dt, matState, prevStateField, QIon, deltaTElectron,
		  deltaTIon);

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

    spProp->getIonSpecificHeat(matState, Cv);
    ionEnergyDeposition = Cv * deltaTIon;
    
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
    
    calcP1Coefs(dt, groupNo, matState, prevStateField, QRad, QElec, QIon,
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
	     const ccsf QRad,
	     const ccsf QElectron,
	     const ccsf QIon,
	     const ccsf TElectron,
	     const ccsf TIon,
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

    // Get sigma emission
    ccsf sigmaEmission;
    spProp->getSigmaEmission(matState, groupNo, sigmaEmission);

    // Get the electron and ion heat capacities.
    
    ccsf CvElec(spMesh);
    spProp->getElectronSpecificHeat(matState, CvElec);

    ccsf CvIon(spMesh);
    spProp->getIonSpecificHeat(matState, CvIon);

    // Calculate Planckian and its temperature derivative

    ccsf planck(spMesh);
    ccsf dPlanckdT(spMesh);
    radPhysics.getPlanck(TElectron, planck);
    radPhysics.getPlanckTemperatureDerivative(TElectron, dPlanckdT);

    // We need our Planckian multiplied by 4pi

    using PhysicalConstants::pi;
    
    planck *= 4.0*pi;
    dPlanckdT *= 4.0*pi;
    
    // We need gamma, the electron-ion coupling coefficient,
    // and CvStar (Cv*) in order to calculate nu.

    ccsf gamma(spMesh);
    spProp->getElectronIonCoupling(matState, gamma);

    // tmpCoeff is a common term to two calculations.
    // Let's just do it once.
    
    const ccsf tmpCoeff = (gamma*dt) / (CvIon + gamma*dt);

    const ccsf CvStar = CvElec + CvIon * tmpCoeff;
    
    // Calculate the "nu" used in the 3T modification of sigmaAbs
    
    const ccsf nu = dt * sigmaEmission * dPlanckdT / (CvStar + dt * dPlanckdT);

    // Calculate QElecStar (Qe*)

    const ccsf QElecStar = QElec +
	(CvIon/dt * (TIon - TElectron) + QIon) * tmpCoeff;

    //
    // We can now calculate the results
    //

    ccsf sigmaAbs(spMesh);
    spProp->getSigmaAbsorption(matState, groupNo, sigmaAbs);
    
    // Calculate the diffusion constant.
    
    D = (1.0/3.0) / (sigmaTotal + tau);

    // Calculate modified sigma absorption
    
    sigmaAbsBar = (1.0 - nu) * sigmaAbs + tau;

    // Calculated modified radiation source

    QRadBar = tau*prevStateField.phi + (1.0 - nu)*sigmaEmission*planck
	+ nu*QElecStar;

    // Calculate the "telegraph" term to the P1 equation.

    Fprime = tau*prevStateField.F / (sigmaTotal + tau);
    
}
 
//---------------------------------------------------------------------------//
//                              end of P13T.cc
//---------------------------------------------------------------------------//
