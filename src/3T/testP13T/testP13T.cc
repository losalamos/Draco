//----------------------------------*-C++-*----------------------------------//
// testP13T.cc
// Randy M. Roberts
// Fri Mar 20 12:04:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testP13T.hh"

#include "3T/testP13T/MeshTypeStub.hh"
#include "3T/testP13T/MaterialPropertiesStub.hh"
#include "3T/testP13T/DiffusionSolverStub.hh"
#include "3T/P13T.hh"
#include "3T/P13TOptions.hh"
#include "3T/Units.hh"
#include "3T/RadiationPhysics.hh"

#include <new>
#include <iostream>
using std::cerr;
using std::endl;

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
testP13T::testP13T()
{
    SP<MT> spmesh = new MT;
    
    SP<DS> spdiffSolver = new DS(spmesh);
    
    Units units = Units::getAstroPhysUnits();

    double TElectron = 2.0;
    double TIon = 4.0;

#if 0
    std::vector<MP::MatVals> matVals(3);
    matVals[0] = MP::MatVals(/*sigmaTotal*/ 2.0,
			     /* sigmaAbsorption */ 1.0,
			     /* sigmaEmission */ 1.0,
			     /* electronIonCoupling */ 1.0e10,
			     /* electronConductionCoeff */ 1.0,
			     /* ionConductionCoeff */ 2.0,
			     /* electronSpecificHeat */ 1.0,
			     /* ionSpecificHeat */ 2.0);
    matVals[1] = MP::MatVals(/*sigmaTotal*/ 2.0,
			     /* sigmaAbsorption */ 1.0,
			     /* sigmaEmission */ 1.0,
			     /* electronIonCoupling */ 1.0,
			     /* electronConductionCoeff */ 1.0,
			     /* ionConductionCoeff */ 2.0,
			     /* electronSpecificHeat */ 1.0,
			     /* ionSpecificHeat */ 2.0);
    matVals[2] = MP::MatVals(/*sigmaTotal*/ 2.0,
			     /* sigmaAbsorption */ 0.0,
			     /* sigmaEmission */ 0.0,
			     /* electronIonCoupling */ 0.0,
			     /* electronConductionCoeff */ 1.0,
			     /* ionConductionCoeff */ 2.0,
			     /* electronSpecificHeat */ 1.0,
			     /* ionSpecificHeat */ 2.0);
    
    SP<MP> spmatprop = new MP(units, spmesh, TElectron, TIon, matVals);
#endif
    
    P13TOptions options(true, true);

    spP13T = new P13TStub(options, spdiffSolver);
}

void testP13T::solve() const
{
    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::bsbf bsbf;
    typedef MT::ncvf ncvf;

    SP<MT> spmesh = spP13T->getMesh();
    SP<MP> spProp = spP13T->getMaterialProperties();

    MP::MaterialStateField<ccsf> matStateCC(spmesh);

    matStateCC = 0;

    // Here is where you would do anything fancy to the cell-temperatures,
    // like averaging, etc., in order to obtain the face-temperatures.
    
    MP::MaterialStateField<fcdsf> matStateFC = matStateCC;
    
    P13T<MT,MP,DS>::RadiationStateField radState(spmesh);

    spP13T->initializeRadiationState(matStateCC, radState);

    ccsf TElec(spmesh);
    ccsf TIon(spmesh);

    spProp->getElectronTemperature(matStateCC, TElec);
    spProp->getIonTemperature(matStateCC, TIon);

    cerr << "TElectron: " << TElec << endl;
    cerr << "TIon: " << TIon << endl;
    cerr << "radState.phi: " << radState.phi
	 << " radState.F: " << radState.F << endl;

    double dt = 1.0e-2;

    ccsf QRad(spmesh);
    ccsf QElectron(spmesh);
    ccsf QIon(spmesh);
    bsbf boundary(spmesh);
    P13T<MT,MP,DS>::RadiationStateField newRadState(spmesh);
    ccsf electEnergyDep(spmesh);
    ccsf ionEnergyDep(spmesh);
    ncvf momDep(spmesh);
    
    QRad = 2.5;
    QElectron = 10.0;
    QIon = 4.25;
    boundary = 0.0;

    spP13T->solve(dt, matStateCC, matStateFC, radState, QRad, QElectron, QIon,
		  boundary, newRadState, electEnergyDep, ionEnergyDep,
		  momDep, TElec, TIon);

    const RadiationPhysics radPhys(spP13T->getMaterialProperties()->getUnits());
    double c = radPhys.getLightSpeed();
    
    ccsf deltaRadEnergy = (newRadState.phi - radState.phi) / c;

    cerr << "deltaRadEnergy: " << deltaRadEnergy
	 << "\trate: " << deltaRadEnergy/dt << endl;
    cerr << "electEnergyDep: " << electEnergyDep
	 << "\trate: " << electEnergyDep/dt << endl;
    cerr << "  ionEnergyDep: " << ionEnergyDep
	 << "\trate: " << ionEnergyDep/dt << endl;
    cerr << "Energy rate: "
	 << (deltaRadEnergy + electEnergyDep + ionEnergyDep) / dt << endl;
    cerr << "newRadState.phi: " << newRadState.phi << endl;
    cerr << "TElectron: " << TElec << endl;
    cerr << "TIon: " << TIon << endl;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of testP13T.cc
//---------------------------------------------------------------------------//
