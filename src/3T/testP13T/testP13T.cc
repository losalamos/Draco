//----------------------------------*-C++-*----------------------------------//
// testP13T.cc
// Randy M. Roberts
// Fri Mar 20 12:04:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testP13T.hh"

#include "matprops/FifiMatPropsReader.hh"
#include "3T/testP13T/MeshTypeStub.hh"
#include "3T/testP13T/DiffusionSolverStub.hh"
#include "3T/P13T.hh"
#include "3T/P13TOptions.hh"
#include "units/Units.hh"
#include "radphys/RadiationPhysics.hh"

#include <new>
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

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

    P13TOptions options(false, false);

    spP13T = new P13T(options, spdiffSolver);
}

void testP13T::solve() const
{
    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::bsbf bsbf;
    typedef MT::ncvf ncvf;

    SP<MT> spmesh = spP13T->getMesh();
    
    Units units = Units::getAstroPhysUnits();

    std::ifstream ifs("testP13T.opac.inp");

    typedef FifiMatPropsReader::MaterialDefinition MatDef;
    vector<MatDef> matdefs;
    matdefs.push_back(MatDef("BeO, 1 Group", 1, 25.0));
    
    FifiMatPropsReader reader(matdefs, units, ifs);

    vector<int> matIds(matdefs.size());
    for (int i=0; i<matdefs.size(); i++)
	matIds[i] = matdefs[i].matid;
    
    MP matProp(matIds, reader);

    double TElect0 = 2.0;
    double TIon0 = 4.0;
    double density = 1.0;
    int matid = 1;
    
    MP::MaterialStateField<ccsf> matStateCC
	= matProp.getMaterialState(ccsf(spmesh, density),
				   ccsf(spmesh, TElect0),
				   ccsf(spmesh, TIon0), ccsf(spmesh, matid));

    // Here is where you would do anything fancy to the cell-temperatures,
    // like averaging, etc., in order to obtain the face-temperatures.
    
    MP::MaterialStateField<fcdsf> matStateFC = matStateCC;
    
    P13T::RadiationStateField radState(spmesh);

    spP13T->initializeRadiationState(matStateCC, radState);

    ccsf TElec(spmesh);
    ccsf TIon(spmesh);
    ccsf CvElec(spmesh);
    ccsf CvIon(spmesh);
    ccsf sigAbs(spmesh);
    ccsf alpha(spmesh);
    ccsf kappaElec(spmesh);
    ccsf kappaIon(spmesh);
    
    matStateCC.getElectronTemperature(TElec);
    matStateCC.getIonTemperature(TIon);
    matStateCC.getElectronSpecificHeat(CvElec);
    matStateCC.getIonSpecificHeat(CvIon);
    matStateCC.getSigmaAbsorption(1,sigAbs);
    matStateCC.getElectronIonCoupling(alpha);
    matStateCC.getElectronConductionCoeff(kappaElec);
    matStateCC.getIonConductionCoeff(kappaIon);

    cout << "TElectron: " << TElec << endl;
    cout << "TIon: " << TIon << endl;
    cout << "Cv Electron: " << CvElec << endl;
    cout << "Cv Ion: " << CvIon << endl;
    cout << "Absorption: " << sigAbs << endl;
    cout << "alpha: " << alpha << endl;
    cout << "kappa Electron: " << kappaElec << endl;
    cout << "kappa Ion: " << kappaIon << endl;
    cout << "radState.phi: " << radState.phi
	 << " radState.F: " << radState.F << endl;

    double dt = 1.0e-2;

    ccsf QRad(spmesh);
    ccsf QElectron(spmesh);
    ccsf QIon(spmesh);
    bsbf boundary(spmesh);
    P13T::RadiationStateField newRadState(spmesh);
    ccsf electEnergyDep(spmesh);
    ccsf ionEnergyDep(spmesh);
    ncvf momDep(spmesh);
    
    QRad = 2.5;
    QElectron = 10.0;
    QIon = 4.25;
    boundary = 0.0;

    spP13T->solve(dt, matStateCC, matStateFC, radState, QRad, QElectron, QIon,
		  boundary, newRadState, electEnergyDep, ionEnergyDep,
		  /* momDep, */ TElec, TIon);

    const RadiationPhysics radPhys(matProp.getUnits());
    double c = radPhys.getLightSpeed();
    
    ccsf deltaRadEnergy = (newRadState.phi - radState.phi) / c;

    cout << "deltaRadEnergy: " << deltaRadEnergy
	 << "\trate: " << deltaRadEnergy/dt << endl;
    cout << "electEnergyDep: " << electEnergyDep
	 << "\trate: " << electEnergyDep/dt << endl;
    cout << "  ionEnergyDep: " << ionEnergyDep
	 << "\trate: " << ionEnergyDep/dt << endl;
    cout << "Energy rate: "
	 << (deltaRadEnergy + electEnergyDep + ionEnergyDep) / dt << endl;
    cout << "newRadState.phi: " << newRadState.phi << endl;
    cout << "TElectron: " << TElec << endl;
    cout << "TIon: " << TIon << endl;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of testP13T.cc
//---------------------------------------------------------------------------//
