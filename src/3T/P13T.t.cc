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

namespace rtt_3T
{
 
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

    // Construct a P1Coeffs object to get derived quantities.
    
    P1Coeffs p1coeffs(*this, dt, groupNo, options, solver, matprops, velocity,
		      prevStateField, QRad, QElectron, QIon, TnElectron, TnIon);

    QEEM = p1coeffs.QEEM();
    REEM = p1coeffs.REEM();
    
    calcNewRadState(resultsStateField, solver, p1coeffs, alpha, beta, bSrc);

    // Calculate the delta electron temperature from the new radiation
    // state, (Te^n+1 - Te^n).
    
    ccsf deltaTElectron(spMesh);
    calcDeltaTElectron(deltaTElectron, dt, numGroups, p1coeffs,
		       resultsStateField);
    
    // Calculate the delta ion temperature from the delta electron
    // temperature, (Ti^n+1 - Ti^n).
    
    ccsf deltaTIon(spMesh);
    calcDeltaTIon(deltaTIon, dt, p1coeffs, deltaTElectron);

    // Calculate electron and ion energy deposition

    ccsf Cv(spMesh);

    matprops.getElectronSpecificHeat(Cv);
    electronEnergyDeposition = Cv * deltaTElectron;

    Tnp1Electron = TnElectron + deltaTElectron;

    matprops.getIonSpecificHeat(Cv);
    ionEnergyDeposition = Cv * deltaTIon;

    Tnp1Ion = TnIon + deltaTIon;
    
    // Calculate the momentum deposition.

    calcMomentumDeposition(momentumDeposition, resultsStateField,
                           solver, matprops, groupNo);

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
// calcNewRadState:
//     calculate the new radiation state using the previous state,
//     material properties, and sources.
//     This solves the coupled radiation, electron, and ion equations
//     ***without*** the conduction equations.
//---------------------------------------------------------------------------//


template<class DS>
void P13T<DS>::calcNewRadState(RadiationStateField &resultsStateField,
			       DiffusionSolver &solver,
			       const P1Coeffs &p1coeffs,
			       const bssf &alpha,
			       const bssf &beta,
			       const bssf &bSrc) const
{
    // Calculate the coefficients needed by the diffusion solver.

    cerr << "sigmaAbsBar: " << *p1coeffs.sigmaAbsBar().begin() << endl;
    cerr << "D: " << *p1coeffs.D().begin() << endl;
    cerr << "QRadBar: " << *p1coeffs.QRadBar().begin() << endl;
    
    // Set up aliases for the long names.
    
    ccsf &phi = resultsStateField.phi;
    FluxField &F = resultsStateField.F;

    // Call the diffusion solver to solve for phi and F.
    // Since phi and F are aliased to the resultsStateField members,
    // this is all we have to do.
    
    solver.solve(phi, F, p1coeffs.D(), p1coeffs.sigmaAbsBar(),
		 p1coeffs.QRadBar(), p1coeffs.Fprime(), alpha, beta,
		 bSrc);

    cerr << "QRadBar/sigmaAbsBar: " <<
	(*p1coeffs.QRadBar().begin())/(*p1coeffs.sigmaAbsBar().begin()) << endl;
    
    cerr << "phi: " << *phi.begin() << endl;

    // ENSURE that F is indeed continuous, if we can!
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
				  const P1Coeffs &p1coeffs,
				  const RadiationStateField &radStateField)
    const
{
    // Only one group.

    Require(numGroups == 1);
    int groupNo = 1;
    
    // Set the radiation physics to the given units.

    const MaterialProperties &matprops = p1coeffs.getMatprops();
    const RadiationPhysics radPhys(matprops.getUnits());

    if (!options.getIsCoupledMaterial())
    {
	// calculate delta T electron
    
	cerr << "QElecStar: " << *p1coeffs.QElecStar().begin() << endl;
	cerr << "CvStar: " << *p1coeffs.CvStar().begin() << endl;
	deltaTelectron = dt * p1coeffs.QElecStar() / p1coeffs.CvStar();
	return;
    }

    ccsf sigmaAbs(spMesh);
    matprops.getSigmaAbsorption(groupNo, sigmaAbs);

    // Get the 4pi*planckian and its temperature derivative
    
    ccsf Bhat(spMesh);
    ccsf dBhatdT(spMesh);

    getBhat(Bhat, radPhys, p1coeffs.getTElectron());
    getdBhatdT(dBhatdT, radPhys, p1coeffs.getTElectron());

    // Get shorthand for phi^n+1
    const ccsf &phi_np1 = radStateField.phi;

    ccsf sigmaEmission(spMesh);
    matprops.getSigmaEmission(groupNo, sigmaEmission);

    // calculate delta T electron
    
    deltaTelectron = dt * (sigmaAbs*phi_np1 - sigmaEmission*Bhat
			   + p1coeffs.QElecStar())
	/ (p1coeffs.CvStar() + dt*sigmaEmission*dBhatdT);
}

//------------------------------------------------------------------------//
// calcDeltaTIon:
//    Calculate the difference between T ion from timestep
//    n+1 to timestep n+1/2
//------------------------------------------------------------------------//

template<class DS>
void P13T<DS>::calcDeltaTIon(ccsf &deltaTIon,
			     double dt,
			     const P1Coeffs &p1coeffs,
			     const ccsf &deltaTelectron) const
{
    // Get the ion heat capacity.
    
    ccsf CvIon(spMesh);
    const MaterialProperties &matprops = p1coeffs.getMatprops();
    matprops.getIonSpecificHeat(CvIon);

    // We need gamma, the electron-ion coupling coefficient.

    ccsf gamma(spMesh);
    matprops.getElectronIonCoupling(gamma);

    deltaTIon = dt * (gamma*(p1coeffs.getTElectron() -
			     p1coeffs.getTIon() + deltaTelectron)
		      + p1coeffs.getQIon())
	/ (CvIon + dt*gamma);

}

//-----------------------------------------------------------------------//
// calcMomentumDeposition:
//    Calculate the momentum deposition from the radiation
//    to the material
//-----------------------------------------------------------------------//

template<class DS>
void P13T<DS>::calcMomentumDeposition
(MomentumField &momentumDeposition,
 const RadiationStateField &resultsStateField,
 const DiffusionSolver &solver,
 const MaterialProperties &matprops, const int groupNo) const
{
    // Set the radiation physics to the given units.

    const RadiationPhysics radPhys(matprops.getUnits());

    // obtain the speed of light

    double c = radPhys.getLightSpeed();

    // Obtain unit vectors
    DiscMomentumField e1Field( spMesh ), e2Field( spMesh ), e3Field( spMesh );
    DiscMomentumField::value_type e1, e2, e3;
    e1(0) = 1.;
    e1(1) = 0.;
    e1(2) = 0.;
    e2(0) = 0.;
    e2(1) = 1.;
    e2(2) = 0.;
    e3(0) = 0.;
    e3(1) = 0.;
    e3(2) = 1.;
    e1Field = e1;
    e2Field = e2;
    e3Field = e3;

    // determine the vertex to node volume ratios

    DiscKineticEnergyField vertex_volumes(spMesh);
    spMesh->get_vertex_volumes(vertex_volumes);
    ncsf node_volumes(spMesh);
    spMesh->get_node_volumes(node_volumes);
    DiscKineticEnergyField vc_node_volumes(spMesh);
    MT::gather ( vc_node_volumes, node_volumes, MT::OpAssign() );
    DiscKineticEnergyField vc_volume_ratios(spMesh);
    vc_volume_ratios = vertex_volumes/vc_node_volumes;

    // calculate momentum deposition

    fcdsf sigmaTotal(spMesh);
    matprops.getSigmaTotal(groupNo, sigmaTotal);
    DiscFluxField sigmaF(spMesh);
    sigmaF = sigmaTotal*resultsStateField.F;
    DiscKineticEnergyField sigmaFe1(spMesh),
                           sigmaFe2(spMesh),
                           sigmaFe3(spMesh);
    solver.dotProduct(sigmaFe1, sigmaF, e1Field);
    solver.dotProduct(sigmaFe2, sigmaF, e2Field);
    solver.dotProduct(sigmaFe3, sigmaF, e3Field);
    // add in the other flux term here

    sigmaFe1 *= vc_volume_ratios;
    sigmaFe2 *= vc_volume_ratios;
    sigmaFe3 *= vc_volume_ratios;
    sigmaFe1 /= c;
    sigmaFe2 /= c;
    sigmaFe3 /= c;

    ncsf momentum1(spMesh), momentum2(spMesh), momentum3(spMesh);
    MT::scatter ( momentum1, sigmaFe1, MT::OpAddAssign() );
    MT::scatter ( momentum2, sigmaFe2, MT::OpAddAssign() );
    MT::scatter ( momentum3, sigmaFe3, MT::OpAddAssign() );
    ncsf::iterator iter1 = momentum1.begin(),
                   iter2 = momentum2.begin(),
                   iter3 = momentum3.begin();
    for (MomentumField::iterator mom_iter = momentumDeposition.begin();
         mom_iter != momentumDeposition.end(); mom_iter++)
    {
        (*mom_iter)(0) = *iter1++;
        (*mom_iter)(1) = *iter2++;
        (*mom_iter)(2) = *iter3++;
    }

#ifdef SDP
    DiscMomentumField vc_momentum(spMesh);
    DiscKineticEnergyField::iterator iter1 = sigmaFe1.begin(),
                                     iter2 = sigmaFe2.begin(),
                                     iter3 = sigmaFe3.begin();
    for (DiscMomentumField::iterator mom_iter = vc_momentum.begin();
         mom_iter != vc_momentum.end(); mom_iter++)
    {
        (*mom_iter)(0) = *iter1++;
        (*mom_iter)(1) = *iter2++;
        (*mom_iter)(2) = *iter3++;
    }
    MT::scatter ( momentumDeposition, vc_momentum, MT::OpAddAssign() );
#endif
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

} // end namespace rtt_3T

#include "3T/P1Coeffs.t.cc"

//---------------------------------------------------------------------------//
//                              end of P13T.cc
//---------------------------------------------------------------------------//
