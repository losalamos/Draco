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
    cerr << "In testP13T()" << endl;

    SP<MT> spmesh = new MT;
    
    SP<DS> spdiffSolver = new DS(spmesh);
    
    Units units = Units::getAstroPhysUnits();
    
    SP<MP> spmatprop = new MP(units, spmesh);

    P13TOptions options(true, true);

    spP13T = new P13TStub(options, spmatprop, spdiffSolver);

    cerr << "spP13T: " << *spP13T << endl;
}

void testP13T::solve() const
{
    cerr << "in testP13T::solve()" << endl;
    
    SP<MT> spmesh = spP13T->getMesh();
    
    MP::MaterialStateField matState(spmesh);
    matState = 1;
    
    P13T<MT,MP,DS>::RadiationStateField radState(spmesh);

    cerr << "before spP13T->initializeRadiationState" << endl;
    
    spP13T->initializeRadiationState(matState, radState);

    cerr << "after spP13T->initializeRadiationState" << endl;
    cerr << "radState.phi: " << radState.phi
	 << " radState.F: " << radState.F << endl;

    double dt = 1.0;

    typedef MT::ccsf ccsf;
    typedef MT::bsbf bsbf;
    typedef MT::ncvf ncvf;

    ccsf QRad(spmesh);
    ccsf QElectron(spmesh);
    ccsf QIon(spmesh);
    bsbf boundary(spmesh);
    P13T<MT,MP,DS>::RadiationStateField newRadState(spmesh);
    ccsf electEnergyDep(spmesh);
    ccsf ionEnergyDep(spmesh);
    ncvf momDep(spmesh);
    
    QRad = 0.0;
    QElectron = 0.0;
    QIon = 0.0;
    boundary = 0.0;

    cerr << "before spP13T->solve()" << endl;
    
    ccsf TElec(spmesh);
    ccsf TIon(spmesh);

    spP13T->solve(dt, matState, radState, QRad, QElectron, QIon,
		  boundary, newRadState, electEnergyDep, ionEnergyDep,
		  momDep, TElec, TIon);
    cerr << "after spP13T->solve()" << endl;

    cerr << "electEnergyDep: " << electEnergyDep << endl;
    cerr << "ionEnergyDep: " << ionEnergyDep << endl;
    cerr << "newRadState.phi: " << newRadState.phi << endl;
    cerr << "TElectron: " << TElec << endl;
    cerr << "TIon: " << TIon << endl;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of testP13T.cc
//---------------------------------------------------------------------------//
