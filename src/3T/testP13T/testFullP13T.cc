//----------------------------------*-C++-*----------------------------------//
// testFullP13T.cc
// Randy M. Roberts
// Sat May 30 15:09:58 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testFullP13T.hh"

#include "matprops/FifiMatPropsReader.hh"
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

using namespace XTM;

std::ostream &operator<<(std::ostream &os, const testFullP13T::MT::ccsf rhs)
{
    typedef testFullP13T::MT::ccsf FT;
    
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
    {
	os << *it << " ";
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const testFullP13T::MT::fcdsf rhs)
{
    typedef testFullP13T::MT::fcdsf FT;
    
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
    {
	os << *it << " ";
    }
    return os;
}

testFullP13T::testFullP13T()
{
    Mesh_DB mdb;
    SP<MT> spmesh = new MT(mdb);
    
    P13TOptions options;

    spP13T = new P13T<MT,MP,DS>(options, spmesh);
}

void testFullP13T::run() const
{
    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::bsbf bsbf;

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

    ccsf TElect0(spmesh);
    ccsf TIon0(spmesh);
    ccsf density(spmesh);
    ccsf matid(spmesh);

    TElect0 = 2.0;
    TIon0 = 4.0;
    density = 1.0;
    matid = 1;

    MP::MaterialStateField<ccsf> matStateCC
	= matProp.getMaterialState(density, TElect0, TIon0,  matid);

    fcdsf TElectFC(spmesh);
    fcdsf TIonFC(spmesh);
    fcdsf densityFC(spmesh);
    fcdsf matidFC(spmesh);

    TElectFC = TElect0;
    TIonFC = TIon0;
    densityFC = density;
    matidFC = matid;
    
    MP::MaterialStateField<fcdsf> matStateFC
	= matProp.getMaterialState(densityFC, TElectFC, TIonFC,  matidFC);

    P13T<MT,MP,DS>::RadiationStateField radState(spmesh);

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
    P13T<MT,MP,DS>::RadiationStateField newRadState(spmesh);
    ccsf electEnergyDep(spmesh);
    ccsf ionEnergyDep(spmesh);
    
    QRad = 2.5;
    QElectron = 10.0;
    QIon = 4.25;
    // boundary = 0.0;

    Diffusion_DB diffdb;
    pcg_DB pcg_db("pcg");
    DS diffSolver(diffdb, spmesh, pcg_db);
    
    spP13T->solve3T(dt, matStateCC, matStateFC, radState, QRad, QElectron, QIon,
		    boundary, diffSolver, newRadState,
		    electEnergyDep, ionEnergyDep, TElec, TIon);

    const RadiationPhysics radPhys(matProp.getUnits());
    double c = radPhys.getLightSpeed();
    
    ccsf deltaRadEnergy(spmesh);
    deltaRadEnergy = (newRadState.phi - radState.phi) / c;

    ccsf rate(spmesh);

    rate = deltaRadEnergy/dt;
    cout << "deltaRadEnergy: " << deltaRadEnergy
	 << "\trate: " << rate << endl;

    rate = electEnergyDep/dt;
    cout << "electEnergyDep: " << electEnergyDep
	 << "\trate: " << rate << endl;

    rate = ionEnergyDep/dt;
    cout << "  ionEnergyDep: " << ionEnergyDep
	 << "\trate: " << rate << endl;

    rate = (deltaRadEnergy + electEnergyDep + ionEnergyDep) / dt;
    cout << "Energy rate: " << rate << endl;
    
    cout << "newRadState.phi: " << newRadState.phi << endl;
    cout << "TElectron: " << TElec << endl;
    cout << "TIon: " << TIon << endl;
}

//---------------------------------------------------------------------------//
//                              end of testFullP13T.cc
//---------------------------------------------------------------------------//
