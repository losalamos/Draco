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
#include "nml/Group.hh"

#include <new>
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
#include <vector>
using std::vector;
#include <string>
using std::string;

using namespace XTM;

namespace
{
    
double sum(const testFullP13T::MT::ccsf rhs)
{
    typedef testFullP13T::MT::ccsf FT;

    double results = 0.0;
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	results += *it;
    return results;
}

}

std::ostream &operator<<(std::ostream &os, const testFullP13T::MT::ccsf rhs)
{
    typedef testFullP13T::MT::ccsf FT;

    std::ios_base::fmtflags fmtflags = os.flags();

    os << std::setw(16) << std::scientific << std::setprecision(6);
    
    int iline = 0;
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
    {
	os << *it << " ";
	if (++iline % 6 == 0)
	    os << endl;
    }
    if (iline % 6 != 0)
	os << endl;

    // restore the original flags.
    
    os.flags(fmtflags);

    return os;
}

std::ostream &operator<<(std::ostream &os, const testFullP13T::MT::fcdsf rhs)
{
    typedef testFullP13T::MT::fcdsf FT;
    
    std::ios_base::fmtflags fmtflags = os.flags();

    os << std::setw(16) << std::scientific << std::setprecision(6);
    
    int iline = 0;
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
    {
	os << *it << " ";
	if (++iline % 6 == 0)
	    os << endl;
    }
    if (iline % 6 != 0)
	os << endl;

    // restore the original flags.
    
    os.flags(fmtflags);

    return os;
}

testFullP13T::testFullP13T(const string &infile)
    : pcg_db("pcg")
{
    NML_Group g("testFullP13T");

    pdb.setup_namelist(g);
    
    diffdb.setup_namelist(g);

    Mesh_DB mdb;
    mdb.setup_namelist(g);
    
    pcg_db.setup_namelist(g);

    cout << "Reading input from " << infile << ".\n";

    g.readgroup(infile.c_str());
    g.writegroup("testFullP13T.out");

    SP<MT> spmesh = new MT(mdb);
    
    P13TOptions options(pdb.P1TauMultiplier, pdb.IsCoupledMaterial);

    spP13T = new P13T<MT,MP,DS>(options, spmesh);
}

void testFullP13T::run() const
{
    typedef MT::ccsf ccsf;
    typedef MT::ccif ccif;
    typedef MT::fcdsf fcdsf;
    typedef MT::fcdif fcdif;
    typedef MT::bssf bssf;

    SP<MT> spmesh = spP13T->getMesh();

    
    Units units;

    switch (pdb.units)
    {
    case AstroPhysical:
	units = Units::getAstroPhysUnits();
	break;
    case SI:
	/* nothing to do. */
	break;
    default:
	throw std::runtime_error("Unknown units specification.");
    }

    std::ifstream ifs(pdb.opacityFile);

    typedef FifiMatPropsReader::MaterialDefinition MatDef;
    vector<MatDef> matdefs;
    matdefs.push_back(MatDef(std::string(pdb.materialName),
			     pdb.materialId, pdb.abar));
    
    FifiMatPropsReader reader(matdefs, units, ifs);

    vector<int> matIds(matdefs.size());
    for (int i=0; i<matdefs.size(); i++)
	matIds[i] = matdefs[i].matid;
    
    MP matProp(matIds, reader);

    ccsf TElect0(spmesh);
    ccsf TIon0(spmesh);
    ccsf density(spmesh);
    ccif matid(spmesh);
    // ccsf matid(spmesh);

    TElect0 = pdb.Te;
    TIon0 = pdb.Ti;
    density = pdb.rho;
    matid = pdb.materialId;

    MP::MaterialStateField<ccsf> matStateCC
	= matProp.getMaterialState(density, TElect0, TIon0,  matid);

    fcdsf TElectFC(spmesh);
    fcdsf TIonFC(spmesh);
    fcdsf densityFC(spmesh);
    fcdif matidFC(spmesh);

    ccsf  oneCC(spmesh);
    oneCC = 1.0;
    
    fcdsf twoFC(spmesh);
    MT::scatter<MT::AddOp>(twoFC, oneCC);

    cout << "twoFC: " << twoFC << endl;
    
    MT::scatter<MT::AddOp>(TElectFC, TElect0);
    TElectFC /= twoFC;
    
    cout << "TElectFC: " << TElectFC << endl;

    MT::scatter<MT::AddOp>(TIonFC, TIon0);
    TIonFC /= twoFC;

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

    double dt = pdb.dt;

    ccsf QRad(spmesh);
    ccsf QElectron(spmesh);
    ccsf QIon(spmesh);
    bssf boundary(spmesh);
    P13T<MT,MP,DS>::RadiationStateField newRadState(spmesh);
    ccsf electEnergyDep(spmesh);
    ccsf ionEnergyDep(spmesh);
    
    QRad = pdb.Qr;
    QElectron = pdb.Qe;
    QIon = pdb.Qi;

    boundary.face(0) = pdb.src_left;
    boundary.face(1) = pdb.src_right;
    boundary.face(2) = pdb.src_front;
    boundary.face(3) = pdb.src_back;
    boundary.face(4) = pdb.src_bottom;
    boundary.face(5) = pdb.src_top;

    cerr << "Made it before diffSolver ctor" << endl;
    
    DS diffSolver(diffdb, spmesh, pcg_db);
    
    cerr << "Made it after diffSolver ctor" << endl;
    
    cerr << "Made it before solve3T" << endl;
    
    spP13T->solve3T(dt, matStateCC, matStateFC, radState, QRad, QElectron, QIon,
		    boundary, diffSolver, newRadState,
		    electEnergyDep, ionEnergyDep, TElec, TIon);

    cerr << "Made it after solve3T" << endl;
    
    cout << "newRadState.phi: " << newRadState.phi << endl;
    cout << "newRadState.F: " << newRadState.F << endl;
    cout << "TElectron: " << TElec << endl;
    cout << "TIon: " << TIon << endl;

    const RadiationPhysics radPhys(matProp.getUnits());
    double c = radPhys.getLightSpeed();
    
    ccsf deltaRadEnergy(spmesh);
    deltaRadEnergy = (newRadState.phi - radState.phi) / c;

    double dx = spmesh->get_dx();
    double dy = spmesh->get_dy();
    double dz = spmesh->get_dz();
    double nx = spmesh->get_ncx();
    double ny = spmesh->get_ncy();
    double nz = spmesh->get_ncz();
    
    double volpcell = dx*dy*dz;
    double volume = volpcell*nx*ny*nz;
    
    cout << "volume: " << volume << endl;
    
    double delta;
    ccsf temp(spmesh);

    delta = sum(deltaRadEnergy)*volpcell;
    cout << "deltaRadEnergy: " << delta
	 << "\trate: " << delta/dt << endl;
    double radrate = delta/dt;

    delta = sum(electEnergyDep)*volpcell;
    cout << "electEnergyDep: " << delta
	 << "\trate: " << delta/dt << endl;

    temp = (TElec - TElect0)*CvElec;
    delta = sum(temp)*volpcell;
    
    // cout << "electEnergyDep (recalc): " << delta
    //      << "\trate: " << delta/dt << endl;

    delta = sum(ionEnergyDep)*volpcell;
    cout << "  ionEnergyDep: " << delta
	 << "\trate: " << delta/dt << endl;

    delta = (sum(deltaRadEnergy) + sum(electEnergyDep) + sum(ionEnergyDep))
           * volpcell;
    cout << "Energy rate: " << delta/dt << endl;

    double qrate = (pdb.Qr + pdb.Qe + pdb.Qi)*volume;

    cout << "Volume src: " << qrate << endl;

    double brate;
    brate  = (pdb.src_left + pdb.src_right) * nx*dx * ny*dy;
    brate += (pdb.src_front + pdb.src_back) * nz*dz * nx*dx;
    brate += (pdb.src_bottom + pdb.src_top) * ny*dy * nz*dz;

    cout << "Boundary src: " << brate << endl;

    double totsrc = brate+qrate;
    cout << "Total src: " << totsrc << endl;

    temp = newRadState.phi*sigAbs;
    double absorption = sum(temp) * volpcell;
    cout << "absorption: " << absorption << endl;

    ccsf planck(spmesh);
    radPhys.getPlanck(TElec, planck);

    temp = sigAbs*planck;
    double emission = sum(temp) * 4.0*PhysicalConstants::pi *
	volpcell;
    cout << "emission: " << emission << endl;

    cout << "emission-absorption: " << emission - absorption << endl;

    double leakage = 0.0;

    for (int k=0; k<nx; k++)
    {
	for (int l=0; l<ny; l++)
	{
	    if (diffdb.alpha_bottom != 0)
	    {
		double phi_bottom =
		    (pdb.src_bottom -
		     diffdb.beta_bottom*newRadState.F(k,l,0,4)) /
		    diffdb.alpha_bottom;
		leakage += phi_bottom/4.0 + newRadState.F(k,l,0,4)/2.0;
	    }

	    if (diffdb.alpha_top != 0)
	    {
		double phi_top = (pdb.src_top -
				  diffdb.beta_top*newRadState.F(k,l,nz-1,5)) /
		    diffdb.alpha_top;
		leakage += phi_top/4.0 + newRadState.F(k,l,nz-1,5)/2.0;
	    }
	}
    }
    leakage *= dx*dy;

    cout << "leakage (x-y): " << leakage << endl;
    
    double balance = radrate + absorption + leakage - totsrc;
    cout << "radrate + absorption + leakage: " <<
	radrate + absorption + leakage << endl;
    cout << "totsrc: " << totsrc << endl;
    cout << "Rad only Balance: " << balance << endl;
    double relbal = balance / radrate;
    cout << "Relative Balance: " << relbal << endl;

    std::ofstream ofs("testFullP13T.dat");
    for (int m=0; m<nz; m++)
	ofs << m << '\t' << newRadState.phi(0,0,m)
	    << '\t' << newRadState.F(0,0,m,4)
	    << '\t' << newRadState.F(0,0,m,5)
	    << endl;
}

//---------------------------------------------------------------------------//
//                              end of testFullP13T.cc
//---------------------------------------------------------------------------//
