//----------------------------------*-C++-*----------------------------------//
// P13T.cc
// Randy M. Roberts
// Wed Mar 11 11:17:54 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/P13T.hh"
#include "radphys/RadiationPhysics.hh"
#include "units/PhysicalConstants.hh"
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
		       const SP<MeshType> &spMesh_)
    : options(options_), spMesh(spMesh_)
{
    // empty
}

//---------------------------------------------------------------------------//
// P13T:
//     Copy constructor
//---------------------------------------------------------------------------//

template<class MT, class MP, class DS>
P13T<MT, MP, DS>::P13T(const P13T<MT, MP, DS> &rhs)
    : options(rhs.options), spMesh(rhs.spMesh)
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
    spMesh = rhs.spMesh;
    return *this;
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

// ACCESSORS

//------------------------------------------------------------------------//
// print:
//     Print itself (for debug mostly)
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
std::ostream &P13T<MT, MP, DS>::print(std::ostream &os) const
{
    os << "(P13T::this: " << (void *)this
       << " spMesh: " << *spMesh
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
initializeRadiationState(const CCMaterialStateField &matStateCC,
			 RadiationStateField &resultsStateField) const
{
    ccsf TElectron(spMesh);
    matStateCC.getElectronTemperature(TElectron);
    
    // Set the radiation physics to the given units.

    const RadiationPhysics radPhys(matStateCC.getUnits());

    // The radiation field is set to the 4pi*planckian with no current flow.

    getBhat(radPhys, TElectron, resultsStateField.phi);

    resultsStateField.F = 0.0;
}

//---------------------------------------------------------------------------//
// solveElectConduction:
//     Solve for the energy deposition and new temperature due to  
//     the conduction equation split.
//---------------------------------------------------------------------------//
    
template<class MT, class MP, class DS>
void
P13T<MT, MP, DS>::solveElectConduction(double dt,
				       const CCMaterialStateField &matStateCC,
				       const FCMaterialStateField &matStateFC,
				       DiffusionSolver &solver,
				       ccsf &electronEnergyDeposition,
				       ccsf &Tnp1Electron) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);
    Require(matStateCC.getUnits() == matStateFC.getUnits());
    
    // This is a one-group problem

    int groupNo = 1;
    int numGroups = 1;

    // Let's save the temperatures at time t^n.
    
    ccsf TnElectron(spMesh);

    matStateCC.getElectronTemperature(TnElectron);

    // Calculate the conduction coefficient
	
    fcdsf kappaElectron(spMesh);
    matStateFC.getElectronConductionCoeff(kappaElectron);
    kappaElectron *= dt;
	
    // Calculate the removal coefficient.
	
    ccsf removalCoeff(spMesh);
    matStateCC.getElectronSpecificHeat(removalCoeff);

    // Calculate the source term.
    ccsf source = removalCoeff * TnElectron;

    // Do the diffusion solve for the temperature at n+1.
	
    solver.solve(kappaElectron, removalCoeff, source, Tnp1Electron);

    // deltaTElectron for radiation will
    // be Te^n+1 - Te^n

    ccsf deltaTElectron = Tnp1Electron - TnElectron;

    // Calculate electron energy deposition

    ccsf Cv(spMesh);

    matStateCC.getElectronSpecificHeat(Cv);
    electronEnergyDeposition = Cv * deltaTElectron;
}

//---------------------------------------------------------------------------//
// solveIonConduction:
//     Solve for the energy deposition and new temperature due to  
//     the conduction equation split.
//---------------------------------------------------------------------------//
    
template<class MT, class MP, class DS>
void
P13T<MT, MP, DS>::solveIonConduction(double dt,
				     const CCMaterialStateField &matStateCC,
				     const FCMaterialStateField &matStateFC,
				     DiffusionSolver &solver,
				     ccsf &ionEnergyDeposition,
				     ccsf &Tnp1Ion) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);
    Require(matStateCC.getUnits() == matStateFC.getUnits());
    
    // This is a one-group problem

    int groupNo = 1;
    int numGroups = 1;

    // Let's save the temperatures at time t^n.
    
    ccsf TnIon(spMesh);

    matStateCC.getIonTemperature(TnIon);

    // Calculate the conduction coefficient
	
    fcdsf kappaIon(spMesh);
    matStateFC.getIonConductionCoeff(kappaIon);
    kappaIon *= dt;
	
    // Calculate the removal coefficient.
	
    ccsf removalCoeff(spMesh);
    matStateCC.getIonSpecificHeat(removalCoeff);

    // Calculate the source term.
    ccsf source = removalCoeff * TnIon;

    // Do the diffusion solve for the temperature at n+1.
	
    solver.solve(kappaIon, removalCoeff, source, Tnp1Ion);

    // deltaTIon for radiation will
    // be Ti^n+1 - Ti^n

    ccsf deltaTIon = Tnp1Ion - TnIon;

    // Calculate ion energy deposition

    ccsf Cv(spMesh);

    matStateCC.getIonSpecificHeat(Cv);
    ionEnergyDeposition = Cv * deltaTIon;
}

//---------------------------------------------------------------------------//
// solve3T:
//     Solve for the new radiation field, the electron/ion energy depositions,
//     and the momentom deposition.
//---------------------------------------------------------------------------//
    
template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::solve3T(double dt,
			       const CCMaterialStateField &matStateCC,
			       const FCMaterialStateField &matStateFC,
			       const RadiationStateField &prevStateField,
			       const ccsf QRad,
			       const ccsf QElectron,
			       const ccsf QIon,
			       const bsbf boundary,
			       DiffusionSolver &solver,
			       RadiationStateField &resultsStateField,
			       ccsf &electronEnergyDeposition,
			       ccsf &ionEnergyDeposition,
#if 0
			       ncvf &momentumDeposition,
#endif
			       ccsf &Tnp1Electron,
			       ccsf &Tnp1Ion) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);
    Require(matStateCC.getUnits() == matStateFC.getUnits());
    
    // This is a one-group problem

    int groupNo = 1;
    int numGroups = 1;

    // Let's save the temperatures at time t^n.
    
    ccsf TnElectron(spMesh);
    ccsf TnIon(spMesh);
    matStateCC.getElectronTemperature(TnElectron);
    matStateCC.getIonTemperature(TnIon);

    // Calculate the new radiation state due to the radiation P1,
    // electron and ion equations ***without*** the conduction equations.
    
    calcNewRadState(dt, groupNo, matStateCC, matStateFC, prevStateField,
		    QRad, QElectron, QIon, TnElectron, TnIon,
		    boundary, solver, resultsStateField);

    // Calculate the delta electron temperature from the new radiation
    // state, (Te^n+1 - Te^n).
    
    ccsf deltaTElectron(spMesh);
    calcDeltaTElectron(dt, numGroups, matStateCC, prevStateField,
		       QElectron, QIon,
		       TnElectron, TnIon,
		       resultsStateField, deltaTElectron);
    
    // Calculate the delta ion temperature from the delta electron
    // temperature, (Ti^n+1 - Ti^n).
    
    ccsf deltaTIon(spMesh);
    calcDeltaTIon(dt, matStateCC, prevStateField, QIon,
		  TnElectron, TnIon,
		  deltaTElectron, deltaTIon);

    // Calculate electron and ion energy deposition

    ccsf Cv(spMesh);

    matStateCC.getElectronSpecificHeat(Cv);
    electronEnergyDeposition = Cv * deltaTElectron;

    Tnp1Electron = TnElectron + deltaTElectron;

    matStateCC.getIonSpecificHeat(Cv);
    ionEnergyDeposition = Cv * deltaTIon;

    Tnp1Ion = TnIon + deltaTIon;
    
#if 0
    // Calculate the momentum deposition.
    momentumDeposition = 0.0;
#endif
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
		const CCMaterialStateField &matStateCC,
		const FCMaterialStateField &matStateFC,
		const RadiationStateField &prevStateField,
		const ccsf QRad,
		const ccsf QElectron,
		const ccsf QIon,
		const ccsf TElectron,
		const ccsf TIon,
		const bsbf boundary,
		DiffusionSolver &solver,
		RadiationStateField &resultsStateField) const
{
    // Calculate the coefficients needed by the diffusion solver.

    fcdsf D(spMesh);
    DiscFluxField Fprime(spMesh);
    ccsf sigmaAbsBar(spMesh);
    ccsf QRadBar(spMesh);
    
    calcP1Coeffs(dt, groupNo, matStateCC, matStateFC, prevStateField,
		 QRad, QElectron, QIon,
		 TElectron, TIon, D, Fprime, sigmaAbsBar, QRadBar);

    // Set up aliases for the long names.
    
    ccsf &phi = resultsStateField.phi;
    FluxField &F = resultsStateField.F;

    // Call the diffusion solver to solve for phi and F.
    // Since phi and F are aliased to the resultsStateField members,
    // this is all we have to do.
    
    solver.solve(D, sigmaAbsBar, QRadBar, Fprime, boundary, phi, F);

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
	     const CCMaterialStateField &matStateCC,
	     const FCMaterialStateField &matStateFC,
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
    
    const RadiationPhysics radPhys(matStateCC.getUnits());

    // set up some needed scalars, like tau

    double c = radPhys.getLightSpeed();
    double tau = 1.0/(c*dt);

    // Ask the material properties for sigma total.
    // It is the material properties' responsibility to do
    // any averaging of temperatures, etc. to achieve the correct
    // resulting sigmaTotal.

    fcdsf sigmaTotal(spMesh);
    matStateFC.getSigmaTotal(groupNo, sigmaTotal);

    //
    // We can now calculate the results
    //

    ccsf sigmaAbs(spMesh);
    matStateCC.getSigmaAbsorption(groupNo, sigmaAbs);
    
    // Calculate the diffusion constant.
    
    D = (1.0/3.0) / (sigmaTotal + tau);

    // We need nu and QElecStar, we get CvStar for free.

    // Get sigma emission
    ccsf sigmaEmission(spMesh);
    matStateCC.getSigmaEmission(groupNo, sigmaEmission);

    ccsf QElecStar(spMesh);
    ccsf CvStar(spMesh);
    ccsf nu(spMesh);
    
    calcStarredFields(dt, groupNo, matStateCC, radPhys,
		      QElectron, QIon, TElectron, TIon,
		      sigmaEmission, QElecStar, CvStar, nu);

    // Calculate modified sigma absorption

    sigmaAbsBar = (1.0 - nu) * sigmaAbs + tau;

    // Calculated modified radiation source

    // We need the Bhat.

    ccsf Bhat(spMesh);
    getBhat(radPhys, TElectron, Bhat);
    
    QRadBar = tau*prevStateField.phi + (1.0 - nu)*sigmaEmission*Bhat
	+ nu*QElecStar + QRad;

    // Calculate the "telegraph" term to the P1 equation.

    Fprime = tau*prevStateField.F / (sigmaTotal + tau);

}

//------------------------------------------------------------------------//
// calcStarredFields:
//    Calculate Qe*, Cv*, and nu.
//    These are needed to calculate other coefficients
//    and delta temperatures.
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::calcStarredFields(double dt,
					 int groupNo,
					 const CCMaterialStateField &matStateCC,
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
    
    calcStarredFields(dt, groupNo, matStateCC, radPhys,
		      QElectron, QIon, TElectron, TIon,
		      sigmaEmission, QElecStar, CvStar);

    // Calculate the 4pi*Planckian's temperature derivative.

    ccsf dBhatdT(spMesh);
    getdBhatdT(radPhys, TElectron, dBhatdT);
    
    // Calculate the "nu" used in the 3T modification of sigmaAbs
    
    nu = dt * sigmaEmission * dBhatdT / (CvStar + dt * dBhatdT);

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
					 const CCMaterialStateField &matStateCC,
					 const RadiationPhysics &radPhys,
					 const ccsf &QElectron,
					 const ccsf &QIon,
					 const ccsf &TElectron,
					 const ccsf &TIon,
					 const ccsf &sigmaEmission,
					 ccsf &QElecStar,
					 ccsf &CvStar) const
{
    // Get the electron and ion heat capacities.
    
    ccsf CvElec(spMesh);
    matStateCC.getElectronSpecificHeat(CvElec);

    ccsf CvIon(spMesh);
    matStateCC.getIonSpecificHeat(CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    matStateCC.getElectronIonCoupling(gamma);

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
					  const CCMaterialStateField &matStateCC, 
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
    
    // Set the radiation physics to the given units.
    
    const RadiationPhysics radPhys(matStateCC.getUnits());

    ccsf sigmaEmission(spMesh);
    matStateCC.getSigmaEmission(groupNo, sigmaEmission);

    // Calculate QElecStar and CvStar.

    ccsf QElecStar(spMesh);
    ccsf CvStar(spMesh);
    
    calcStarredFields(dt, groupNo, matStateCC, radPhys,
		      QElectron, QIon, TElectron, TIon,
		      sigmaEmission, QElecStar, CvStar);

    ccsf sigmaAbs(spMesh);
    matStateCC.getSigmaAbsorption(groupNo, sigmaAbs);

    // Get the 4pi*planckian and its temperature derivative
    
    ccsf Bhat(spMesh);
    ccsf dBhatdT(spMesh);

    getBhat(radPhys, TElectron, Bhat);
    getdBhatdT(radPhys, TElectron, dBhatdT);

    // Get shorthand for phi^n+1
    const ccsf &phi_np1 = resultsStateField.phi;

    // calculate delta T electron
    
    deltaTelectron = dt * (sigmaAbs*phi_np1 - sigmaEmission*Bhat + QElecStar)
	/ (CvStar + dt*sigmaEmission*dBhatdT);
}

//------------------------------------------------------------------------//
// calcDeltaTIon:
//    Calculate the difference between T ion from timestep
//    n+1 to timestep n+1/2
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::calcDeltaTIon(double dt,
				     const CCMaterialStateField &matStateCC,
				     const RadiationStateField &prevStateField, 
				     const ccsf &QIon,
				     const ccsf &TElectron,
				     const ccsf &TIon,
				     const ccsf &deltaTelectron,
				     ccsf &deltaTIon) const
{
    // Get the ion heat capacity.
    
    ccsf CvIon(spMesh);
    matStateCC.getIonSpecificHeat(CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    matStateCC.getElectronIonCoupling(gamma);

    deltaTIon = dt * (gamma*(TElectron - TIon + deltaTelectron) + QIon)
	/ (CvIon + dt*gamma);

}

//------------------------------------------------------------------------//
// getBhat:
//    get the 4pi*planckian
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::getBhat(const RadiationPhysics &radPhys,
			       const ccsf &TElectron, ccsf &Bhat) const
{
    radPhys.getPlanck(TElectron, Bhat);

    // We need our Planckian multiplied by 4pi

    using XTM::PhysicalConstants::pi;
    Bhat *= 4.0*pi;
}

//------------------------------------------------------------------------//
// getdBhatdT:
//    get the 4pi*dPlanckiandT
//------------------------------------------------------------------------//

template<class MT, class MP, class DS>
void P13T<MT, MP, DS>::getdBhatdT(const RadiationPhysics &radPhys,
				  const ccsf &TElectron, ccsf &dBhatdT) const
{
    radPhys.getPlanckTemperatureDerivative(TElectron, dBhatdT);

    // We need our Planckian multiplied by 4pi

    using XTM::PhysicalConstants::pi;
    dBhatdT *= 4.0*pi;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of P13T.cc
//---------------------------------------------------------------------------//
