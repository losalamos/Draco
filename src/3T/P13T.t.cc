//----------------------------------*-C++-*----------------------------------//
// P13T.cc
// Randy M. Roberts
// Wed Mar 11 11:17:54 1998
//---------------------------------------------------------------------------//
// @> P1 3T Radiation Solver.
//---------------------------------------------------------------------------//

#include "3T/P13T.hh"
#include "radphys/RadiationPhysics.hh"
#include "units/PhysicalConstants.hh"
#include "timestep/ts_manager.hh"
#include "timestep/field_ts_advisor.hh"
#include "matprops/MaterialProperties.hh"
#include <cmath>

#include <iostream>
using std::cerr;
using std::endl;

using XTM::P13T;

// CREATORS

//---------------------------------------------------------------------------//
// P13T:
//     Constructor
//---------------------------------------------------------------------------//

template<class DS>
P13T<DS>::P13T(const P13TOptions &options_, const dsxx::SP<MeshType> &spMesh_)
    : options(options_), spMesh(spMesh_)
{
    // empty
}

//---------------------------------------------------------------------------//
// P13T:
//     Constructor
//---------------------------------------------------------------------------//

template<class DS>
P13T<DS>::P13T(const P13TOptions &options_, const dsxx::SP<MeshType> &spMesh_,
	       dsxx::SP<ts_manager> &spTsManager_)
    : options(options_), spMesh(spMesh_), spTsManager(spTsManager_)
{
    // Set up timestep advisors, and add them to the manager.

    Require(spTsManager);

    // We need to add a way to use non-default values on the construction
    // of these advisors.

    // We will make sure that they are not activated until used.
	
    spRadTsAdvisor = new field_ts_advisor("P13T Radiation Intensity");
    spRadTsAdvisor->deactivate();
	
    spElecTsAdvisor = new field_ts_advisor("P13T Electron Temperature");
    spElecTsAdvisor->deactivate();
	
    spIonTsAdvisor = new field_ts_advisor("P13T Ion Temperature");
    spIonTsAdvisor->deactivate();
	
    spElecCondTsAdvisor = new field_ts_advisor("P13T Elec Conduction");
    spElecCondTsAdvisor->deactivate();
	
    spIonCondTsAdvisor = new field_ts_advisor("P13T Ion Conduction");
    spIonCondTsAdvisor->deactivate();

    spTsManager->add_advisor(spRadTsAdvisor);
    spTsManager->add_advisor(spElecTsAdvisor);
    spTsManager->add_advisor(spIonTsAdvisor);
    spTsManager->add_advisor(spElecCondTsAdvisor);
    spTsManager->add_advisor(spIonCondTsAdvisor);
}

//---------------------------------------------------------------------------//
// ~P13T:
//     Destructor
//---------------------------------------------------------------------------//

template<class DS>
P13T<DS>::~P13T()
{
    // If we have a timestep manager, we must remove the registered
    // timestep advisors.
    
    if (spTsManager)
    {
	spTsManager->remove_advisor(spRadTsAdvisor);
	spTsManager->remove_advisor(spElecTsAdvisor);
	spTsManager->remove_advisor(spIonTsAdvisor);
	spTsManager->remove_advisor(spElecCondTsAdvisor);
	spTsManager->remove_advisor(spIonCondTsAdvisor);
    }
}



// MANIPULATORS

//---------------------------------------------------------------------------//
// setOptions:
//     Set the contained options object to the new value.
//---------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::setOptions(const P13TOptions options_)
{
    options = options_;
}

// ACCESSORS

//------------------------------------------------------------------------//
// print:
//     Print itself (for debug mostly)
//------------------------------------------------------------------------//

template<class DS>
std::ostream &P13T<DS>::print(std::ostream &os) const
{
    os << "(P13T::this: " << (void *)this
       << ")";
    return os;
}

//---------------------------------------------------------------------------//
// initializeRadiationState:
//     Initialize the radiation field to Planckian
//     based on material electron temperatures.
//---------------------------------------------------------------------------//
    
template<class DS>
void P13T<DS>::
initializeRadiationState(const MaterialProperties &matprops,
			 RadiationStateField &resultsStateField) const
{
    ccsf TElectron(spMesh);
    matprops.getElectronTemperature(TElectron);
    
    // Set the radiation physics to the given units.

    const RadiationPhysics radPhys(matprops.getUnits());

    // The radiation field is set to the 4pi*planckian with no current flow.

    getBhat(resultsStateField.phi, radPhys, TElectron);

    resultsStateField.F = 0.0;
}

//---------------------------------------------------------------------------//
// solveElectConduction:
//     Solve for the energy deposition and new temperature due to  
//     the conduction equation split.
//---------------------------------------------------------------------------//
    
template<class DS>
void P13T<DS>::solveElectConduction(ccsf &electronEnergyDeposition,
				    ccsf &Tnp1Electron,
				    DiffusionSolver &solver,
				    double dt,
				    const MaterialProperties &matprops,
				    const bssf &alpha,
				    const bssf &beta,
				    const bssf &bSrc) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);

    // Let's save the temperatures at time t^n.
    
    ccsf TnElectron(spMesh);

    matprops.getElectronTemperature(TnElectron);

    // Calculate the conduction coefficient
	
    fcdsf kappaElectron(spMesh);
    matprops.getElectronConductionCoeff(kappaElectron);
    kappaElectron = kappaElectron * dt;
	
    // Calculate the removal coefficient.
	
    ccsf removalCoeff(spMesh);
    matprops.getElectronSpecificHeat(removalCoeff);

    // Calculate the source term.
    ccsf source(spMesh);
    source = removalCoeff * TnElectron;

    // Do the diffusion solve for the temperature at n+1.

    FluxField Flux(spMesh);
    const DiscFluxField Fprime(spMesh);
    solver.solve(Tnp1Electron, Flux, kappaElectron, removalCoeff, source,
		 Fprime, alpha, beta, bSrc);

    // deltaTElectron for radiation will
    // be Te^n+1 - Te^n

    ccsf deltaTElectron(spMesh);
    deltaTElectron = Tnp1Electron - TnElectron;

    // Calculate electron energy deposition

    ccsf Cv(spMesh);

    matprops.getElectronSpecificHeat(Cv);
    electronEnergyDeposition = Cv * deltaTElectron;

    // Let's activate and update the electron conduction timestep advisor.
    
    if (spElecCondTsAdvisor)
    {
	Assert(spTsManager);
	spElecCondTsAdvisor->activate();
	spElecCondTsAdvisor->update_tstep(*spTsManager, TnElectron,
					  Tnp1Electron);
    }
}

//---------------------------------------------------------------------------//
// solveIonConduction:
//     Solve for the energy deposition and new temperature due to  
//     the conduction equation split.
//---------------------------------------------------------------------------//
    
template<class DS>
void P13T<DS>::solveIonConduction(ccsf &ionEnergyDeposition,
				  ccsf &Tnp1Ion,
				  DiffusionSolver &solver,
				  double dt,
				  const MaterialProperties &matprops,
				  const bssf &alpha,
				  const bssf &beta,
				  const bssf &bSrc) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);
    
    // Let's save the temperatures at time t^n.
    
    ccsf TnIon(spMesh);

    matprops.getIonTemperature(TnIon);

    // Calculate the conduction coefficient
	
    fcdsf kappaIon(spMesh);
    matprops.getIonConductionCoeff(kappaIon);
    kappaIon = kappaIon * dt;
	
    // Calculate the removal coefficient.
	
    ccsf removalCoeff(spMesh);
    matprops.getIonSpecificHeat(removalCoeff);

    // Calculate the source term.
    ccsf source(spMesh);
    source = removalCoeff * TnIon;

    // Do the diffusion solve for the temperature at n+1.
	
    FluxField Flux(spMesh);
    const DiscFluxField Fprime(spMesh);
    solver.solve(Tnp1Ion, Flux, kappaIon, removalCoeff, source,
		 Fprime, alpha, beta, bSrc);

    // deltaTIon for radiation will
    // be Ti^n+1 - Ti^n

    ccsf deltaTIon(spMesh);
    deltaTIon = Tnp1Ion - TnIon;

    // Calculate ion energy deposition

    ccsf Cv(spMesh);

    matprops.getIonSpecificHeat(Cv);
    ionEnergyDeposition = Cv * deltaTIon;

    // Let's activate and update the ion conduction timestep advisor.
    
    if (spIonCondTsAdvisor)
    {
	Assert(spTsManager);
	spIonCondTsAdvisor->activate();
	spIonCondTsAdvisor->update_tstep(*spTsManager, TnIon, Tnp1Ion);
    }
}

//---------------------------------------------------------------------------//
// solve3T:
//     Solve for the new radiation field, the electron/ion energy depositions,
//     and the momentom deposition.
//---------------------------------------------------------------------------//
    
template<class DS>
void P13T<DS>::solve3T(RadiationStateField &resultsStateField,
		       ccsf &QEEM,
		       ccsf &REEM,
		       ccsf &electronEnergyDeposition,
		       ccsf &ionEnergyDeposition,
		       ncvsf &momentumDeposition,
		       ccsf &Tnp1Electron,
		       ccsf &Tnp1Ion,
		       DiffusionSolver &solver,
		       double dt,
		       const MaterialProperties &matprops,
		       const ncvsf &velocity,
		       const RadiationStateField &prevStateField,
		       const ccsf &QRad,
		       const ccsf &QElectron,
		       const ccsf &QIon,
		       const bssf &alpha,
		       const bssf &beta,
		       const bssf &bSrc) const
{
    // Require dt > 0, etc.

    Require(dt > 0.0);
    
    // This is a one-group problem

    int groupNo = 1;
    int numGroups = 1;

    // Let's save the temperatures at time t^n.
    
    ccsf TnElectron(spMesh);
    ccsf TnIon(spMesh);
    matprops.getElectronTemperature(TnElectron);
    matprops.getIonTemperature(TnIon);

    // Calculate the new radiation state due to the radiation P1,
    // electron and ion equations ***without*** the conduction equations.
    
    calcNewRadState(resultsStateField, QEEM, REEM, solver,
		    dt, groupNo, matprops, velocity, prevStateField,
		    QRad, QElectron, QIon, TnElectron, TnIon, alpha, beta,
		    bSrc);

    // Calculate the delta electron temperature from the new radiation
    // state, (Te^n+1 - Te^n).
    
    ccsf deltaTElectron(spMesh);
    calcDeltaTElectron(deltaTElectron, dt, numGroups, matprops,
		       prevStateField,
		       QElectron, QIon,
		       TnElectron, TnIon,
		       resultsStateField);
    
    // Calculate the delta ion temperature from the delta electron
    // temperature, (Ti^n+1 - Ti^n).
    
    ccsf deltaTIon(spMesh);
    calcDeltaTIon(deltaTIon, dt, matprops, prevStateField, QIon,
		  TnElectron, TnIon,
		  deltaTElectron);

    // Calculate electron and ion energy deposition

    ccsf Cv(spMesh);

    matprops.getElectronSpecificHeat(Cv);
    electronEnergyDeposition = Cv * deltaTElectron;

    Tnp1Electron = TnElectron + deltaTElectron;

    matprops.getIonSpecificHeat(Cv);
    ionEnergyDeposition = Cv * deltaTIon;

    Tnp1Ion = TnIon + deltaTIon;
    
#if 0
    // Calculate the momentum deposition.
    momentumDeposition = 0.0;
#endif

    // Update and activate the timestep advisors.
    
    if (spTsManager)
    {
	Assert(spRadTsAdvisor);
	Assert(spElecTsAdvisor);
	Assert(spIonTsAdvisor);
	
	spRadTsAdvisor->activate();
	spRadTsAdvisor->update_tstep(*spTsManager, prevStateField.phi,
				     resultsStateField.phi);
	
	spElecTsAdvisor->activate();
	spElecTsAdvisor->update_tstep(*spTsManager, TnElectron, Tnp1Electron);
	
	spIonTsAdvisor->activate();
	spIonTsAdvisor->update_tstep(*spTsManager, TnIon, Tnp1Ion);
    }
}

//---------------------------------------------------------------------------//
// clacNewRadState:
//     calculate the new radiation state using the previous state,
//     material properties, and sources.
//     This solves the coupled radiation, electron, and ion equations
//     ***without*** the conduction equations.
//---------------------------------------------------------------------------//


template<class DS>
void P13T<DS>::calcNewRadState(RadiationStateField &resultsStateField,
			       ccsf &QEEM,
			       ccsf &REEM,
			       DiffusionSolver &solver,
			       double dt,
			       int groupNo,
			       const MaterialProperties &matprops,
			       const ncvsf &velocity,
			       const RadiationStateField &prevStateField,
			       const ccsf &QRad,
			       const ccsf &QElectron,
			       const ccsf &QIon,
			       const ccsf &TElectron,
			       const ccsf &TIon,
			       const bssf &alpha,
			       const bssf &beta,
			       const bssf &bSrc) const
{
    // Calculate the coefficients needed by the diffusion solver.

    fcdsf D(spMesh);
    DiscFluxField Fprime(spMesh);
    ccsf sigmaAbsBar(spMesh);
    ccsf QRadBar(spMesh);
    
    calcP1Coeffs(D, Fprime, sigmaAbsBar, QEEM, REEM, QRadBar,
		 dt, groupNo, matprops, velocity, prevStateField,
		 QRad, QElectron, QIon,
		 TElectron, TIon);

    cerr << "sigmaAbsBar: " << *sigmaAbsBar.begin() << endl;
    cerr << "D: " << *D.begin() << endl;
    cerr << "QRadBar: " << *QRadBar.begin() << endl;
    
    // Set up aliases for the long names.
    
    ccsf &phi = resultsStateField.phi;
    FluxField &F = resultsStateField.F;

    // Call the diffusion solver to solve for phi and F.
    // Since phi and F are aliased to the resultsStateField members,
    // this is all we have to do.
    
    solver.solve(phi, F, D, sigmaAbsBar, QRadBar, Fprime, alpha, beta,
		 bSrc);

    cerr << "QRadBar/sigmaAbsBar: " <<
	(*QRadBar.begin())/(*sigmaAbsBar.begin()) << endl;
    
    cerr << "phi: " << *phi.begin() << endl;

    // ENSURE that F is indeed continuous, if we can!
}
    
//---------------------------------------------------------------------------//
// calcP1Coeffs:
//     Calculate the coefficients, e.g. diffusion and removal, and
//     source terms for solving the P1 equation.
//---------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::calcP1Coeffs(fcdsf &D,
			    DiscFluxField &Fprime,
			    ccsf &sigmaAbsBar,
			    ccsf &QEEM,
			    ccsf &REEM,
			    ccsf &QRadBar,
			    double dt,
			    int groupNo,
			    const MaterialProperties &matprops,
			    const ncvsf &velocity,
			    const RadiationStateField &prevStateField,
			    const ccsf &QRad,
			    const ccsf &QElectron,
			    const ccsf &QIon,
			    const ccsf &TElectron,
			    const ccsf &TIon) const
{
    // Set the radiation physics to the given units.
    
    const RadiationPhysics radPhys(matprops.getUnits());

    // set up some needed scalars, like tau

    double c = radPhys.getLightSpeed();
    double tau = 1.0/(c*dt);

    // If this is diffusion, then tauP1 = 0.0.
    // If this is full P1 then tauP1 = tau.
    
    double tauP1 = tau * options.getP1TauMultiplier();

    // Ask the material properties for sigma total.
    // It is the material properties' responsibility to do
    // any averaging of temperatures, etc. to achieve the correct
    // resulting sigmaTotal.

    fcdsf sigmaTotal(spMesh);
    matprops.getSigmaTotal(groupNo, sigmaTotal);

    //
    // We can now calculate the results
    //

    ccsf sigmaAbs(spMesh);
    matprops.getSigmaAbsorption(groupNo, sigmaAbs);
    
    // Calculate the diffusion constant.

    D = (1.0/3.0) / (sigmaTotal + tauP1);

    // We need nu and QElecStar, we get CvStar for free.

    // Get sigma emission
    ccsf sigmaEmission(spMesh);
    matprops.getSigmaEmission(groupNo, sigmaEmission);

    ccsf QElecStar(spMesh);
    ccsf CvStar(spMesh);
    ccsf nu(spMesh);

    if (options.getIsCoupledMaterial())
    {
	calcStarredFields(QElecStar, CvStar, nu,
			  dt, groupNo, matprops, radPhys,
			  QElectron, QIon, TElectron, TIon,
			  sigmaEmission);

	cerr << "QElecStar: " << *QElecStar.begin() << endl;
	cerr << "CvStar: " << *CvStar.begin() << endl;
	cerr << "nu: " << *nu.begin() << endl;
    
	// Calculate modified sigma absorption

	sigmaAbsBar = (1.0 - nu) * sigmaAbs + tau;

	// Calculate the emmissive removal coefficient
	
	REEM = -1.0 * nu * sigmaAbs;
	
	// We need the Bhat.

	ccsf Bhat(spMesh);
	getBhat(Bhat, radPhys, TElectron);

	// Calculate the emmissive source term.
	
	QEEM = (1.0 - nu)*sigmaEmission*Bhat + nu*QElecStar;

	// Calculated modified radiation source

	QRadBar = tau*prevStateField.phi + QEEM + QRad;
    }
    else
    {
	sigmaAbsBar = sigmaAbs + tau;
	REEM = 0.0;
	QEEM = 0.0;
	QRadBar = tau*prevStateField.phi + QRad;
    }

    // Calculate the "telegraph" term to the P1 equation.

    Fprime = tauP1*prevStateField.F / (sigmaTotal + tauP1);

}

//------------------------------------------------------------------------//
// calcStarredFields:
//    Calculate Qe*, Cv*, and nu.
//    These are needed to calculate other coefficients
//    and delta temperatures.
//------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::calcStarredFields(ccsf &QElecStar,
				 ccsf &CvStar,
				 ccsf &nu,
				 double dt,
				 int groupNo,
				 const MaterialProperties &matprops,
				 const RadiationPhysics &radPhys,
				 const ccsf &QElectron,
				 const ccsf &QIon,
				 const ccsf &TElectron,
				 const ccsf &TIon,
				 const ccsf &sigmaEmission) const
{
    // Calculate Qe* and Cv*.
    // We will then calculate nu ourself.
    
    calcStarredFields(QElecStar, CvStar,
		      dt, groupNo, matprops, radPhys,
		      QElectron, QIon, TElectron, TIon);

    // Calculate the 4pi*Planckian's temperature derivative.

    ccsf dBhatdT(spMesh);
    getdBhatdT(dBhatdT, radPhys, TElectron);
    
    // Calculate the "nu" used in the 3T modification of sigmaAbs
    
    nu = dt * sigmaEmission * dBhatdT / (CvStar + dt * sigmaEmission * dBhatdT);

}

//------------------------------------------------------------------------//
// calcStarredFields:
//    Calculate Qe*, Cv*, but not nu.
//    These are needed to calculate other coefficients
//    and delta temperatures.
//------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::calcStarredFields(ccsf &QElecStar,
				 ccsf &CvStar,
				 double dt,
				 int groupNo,
				 const MaterialProperties &matprops,
				 const RadiationPhysics &radPhys,
				 const ccsf &QElectron,
				 const ccsf &QIon,
				 const ccsf &TElectron,
				 const ccsf &TIon) const
{
    // Get the electron and ion heat capacities.
    
    ccsf CvElec(spMesh);
    matprops.getElectronSpecificHeat(CvElec);

    ccsf CvIon(spMesh);
    matprops.getIonSpecificHeat(CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    matprops.getElectronIonCoupling(gamma);

    // tmpCoeff is a common term to two calculations.
    // Let's just do it once.
    
    ccsf tmpCoeff(spMesh);
    // tmpCoeff = (gamma*dt) / (CvIon + gamma*dt);
    tmpCoeff = gamma / (CvIon/dt + gamma);

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

template<class DS>
void P13T<DS>::calcDeltaTElectron(ccsf &deltaTelectron,
				  double dt,
				  int numGroups, 
				  const MaterialProperties &matprops, 
				  const RadiationStateField &prevStateField,
				  const ccsf &QElectron, 
				  const ccsf &QIon,
				  const ccsf &TElectron,
				  const ccsf &TIon,
				  const RadiationStateField &resultsStateField) const
{
    // Only one group.

    Require(numGroups == 1);
    int groupNo = 1;
    
    // Set the radiation physics to the given units.
    
    const RadiationPhysics radPhys(matprops.getUnits());

    // Calculate QElecStar and CvStar.

    ccsf QElecStar(spMesh);
    ccsf CvStar(spMesh);
    
    calcStarredFields(QElecStar, CvStar, dt, groupNo, matprops, radPhys,
		      QElectron, QIon, TElectron, TIon);

    if (!options.getIsCoupledMaterial())
    {
	// calculate delta T electron
    
	cerr << "QElecStar: " << *QElecStar.begin() << endl;
	cerr << "CvStar: " << *CvStar.begin() << endl;
	deltaTelectron = dt * QElecStar / CvStar;
	return;
    }

    ccsf sigmaAbs(spMesh);
    matprops.getSigmaAbsorption(groupNo, sigmaAbs);

    // Get the 4pi*planckian and its temperature derivative
    
    ccsf Bhat(spMesh);
    ccsf dBhatdT(spMesh);

    getBhat(Bhat, radPhys, TElectron);
    getdBhatdT(dBhatdT, radPhys, TElectron);

    // Get shorthand for phi^n+1
    const ccsf &phi_np1 = resultsStateField.phi;

    ccsf sigmaEmission(spMesh);
    matprops.getSigmaEmission(groupNo, sigmaEmission);

    // calculate delta T electron
    
    deltaTelectron = dt * (sigmaAbs*phi_np1 - sigmaEmission*Bhat + QElecStar)
	/ (CvStar + dt*sigmaEmission*dBhatdT);
}

//------------------------------------------------------------------------//
// calcDeltaTIon:
//    Calculate the difference between T ion from timestep
//    n+1 to timestep n+1/2
//------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::calcDeltaTIon(ccsf &deltaTIon,
			     double dt,
			     const MaterialProperties &matprops,
			     const RadiationStateField &prevStateField, 
			     const ccsf &QIon,
			     const ccsf &TElectron,
			     const ccsf &TIon,
			     const ccsf &deltaTelectron) const
{
    // Get the ion heat capacity.
    
    ccsf CvIon(spMesh);
    matprops.getIonSpecificHeat(CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    matprops.getElectronIonCoupling(gamma);

    deltaTIon = dt * (gamma*(TElectron - TIon + deltaTelectron) + QIon)
	/ (CvIon + dt*gamma);

}

//------------------------------------------------------------------------//
// getBhat:
//    get the 4pi*planckian
//------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::getBhat(ccsf &Bhat, const RadiationPhysics &radPhys,
		       const ccsf &TElectron) const
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

template<class DS>
void P13T<DS>::getdBhatdT(ccsf &dBhatdT,
			  const RadiationPhysics &radPhys,
			  const ccsf &TElectron) const
{
    radPhys.getPlanckTemperatureDerivative(TElectron, dBhatdT);

    // We need our Planckian multiplied by 4pi

    using XTM::PhysicalConstants::pi;
    dBhatdT *= 4.0*pi;
}

//---------------------------------------------------------------------------//
//                              end of P13T.cc
//---------------------------------------------------------------------------//
