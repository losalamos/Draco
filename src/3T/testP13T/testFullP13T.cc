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
#include "timestep/ts_manager.hh"
#include "timestep/fixed_ts_advisor.hh"
#include "timestep/ratio_ts_advisor.hh"

#include <functional>
#include <new>
#include <strstream>
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
#include <vector>
using std::vector;
#include <string>
using std::string;

using dsxx::SP;

using namespace XTM;

namespace
{

 int nx;
 int ny;
 int nz;

 double sum(const testFullP13T::MT::ccsf &rhs)
 {
     typedef testFullP13T::MT::ccsf FT;

     double results = 0.0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 results += *it;
     return results;
 }

 template<class FT>
     typename FT::value_type min(const FT &rhs)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (*it < results)
	     results = *it;
     return results;
 }

 template<class FT>
     typename FT::value_type max(const FT &rhs)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (*it > results)
	     results = *it;
     return results;
 }

 template<class FT, class OP>
     typename FT::value_type min(const FT &rhs, OP &op)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (op(*it) < results)
	     results = *it;
     return results;
 }

 template<class FT, class OP>
     typename FT::value_type max(const FT &rhs, OP &op)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (op(*it) > results)
	     results = *it;
     return results;
 }

 class GmvDump
 {
   private:
     
     struct VertId
     {
	 int nx;
	 int ny;
	 int nz;
	 VertId(int nx_, int ny_, int nz_) : nx(nx_+1), ny(ny_+1), nz(nz_+1) { }
	 int operator()(int i, int j, int k)
	 {
	     return k*nx*ny + j*nx + i + 1;
	 }
     };

     std::ostream &os;
     const SP<testFullP13T::MT> &spMesh;

     int nx;
     int ny;
     int nz;

     VertId vid;

     bool variablePrinted;

     int cycle;
     double time;
     
   public:
     
     GmvDump(std::ostream &os_, const SP<testFullP13T::MT> &spMesh_,
	     int cycle_, double time_)
	 : os(os_), spMesh(spMesh_), nx(spMesh->get_ncx()),
	   ny(spMesh->get_ncy()), nz(spMesh->get_ncz()), vid(nx, ny, nz),
	   variablePrinted(false), cycle(cycle_), time(time_)
     {
	 std::ios_base::fmtflags fmtflags = os.flags();

	 os << std::setw(16) << std::scientific << std::setprecision(6);

	 os << "gmvinput ascii" << endl;

	 double dx = spMesh->get_dx();
	 double dy = spMesh->get_dy();
	 double dz = spMesh->get_dz();

	 int nnodes = (nx+1)*(ny+1)*(nz+1);
	 os << "nodes " << nnodes << endl;

	 const int nitemsPerLine = 5;
	 int nitems = 0;

	 for (int k=0; k<nz+1; k++)
	 {
	     for (int j=0; j<ny+1; j++)
	     {
		 for (int i=0; i<nx+1; i++)
		 {
		     os << i*dx;
		     if (++nitems % nitemsPerLine)
			 os << " ";
		     else
			 os << endl;
		 }
	     }
	 }
	 if (nitems % nitemsPerLine)
	     os << endl;

	 nitems = 0;
	 for (int k=0; k<nz+1; k++)
	 {
	     for (int j=0; j<ny+1; j++)
	     {
		 for (int i=0; i<nx+1; i++)
		 {
		     os << j*dy;
		     if (++nitems % nitemsPerLine)
			 os << " ";
		     else
			 os << endl;
		 }
	     }
	 }
	 if (nitems % nitemsPerLine)
	     os << endl;

	 nitems = 0;
	 for (int k=0; k<nz+1; k++)
	 {
	     for (int j=0; j<ny+1; j++)
	     {
		 for (int i=0; i<nx+1; i++)
		 {
		     os << k*dz;
		     if (++nitems % nitemsPerLine)
			 os << " ";
		     else
			 os << endl;
		 }
	     }
	 }
	 if (nitems % nitemsPerLine)
	     os << endl;

	 os << "cells " << nx*ny*nz << endl;

	 for (int k=0; k<nz; k++)
	 {
	     for (int j=0; j<ny; j++)
	     {
		 for (int i=0; i<nx; i++)
		 {
		     os << "hex 8" << endl;
		     os << vid(i,   j,   k+1) << " ";
		     os << vid(i+1, j,   k+1) << " ";
		     os << vid(i+1, j+1, k+1) << " ";
		     os << vid(i,   j+1, k+1) << " ";
		     os << vid(i,   j,   k  ) << " ";
		     os << vid(i+1, j,   k  ) << " ";
		     os << vid(i+1, j+1, k  ) << " ";
		     os << vid(i,   j+1, k  ) << " ";
		     os << endl;
		 }
	     }
	 }

	 // restore the original flags.
	 
	 os.flags(fmtflags);
     }
 
     ~GmvDump()
     {
	 if (variablePrinted)
	     os << "endvars" << endl;
	 os << "probtime " << time << endl;
	 os << "cycleno " << cycle << endl;
	 os << "endgmv" << endl;
     }

     void dump(const testFullP13T::MT::ccsf &var, const string &name)
     {
	 std::ios_base::fmtflags fmtflags = os.flags();

	 os << std::setw(16) << std::scientific << std::setprecision(6);
	 
	 if (!variablePrinted)
	 {
	     os << "variable" << endl;
	     variablePrinted = true;
	 }
	 
	 os << name << " 0" << endl;

	 const int nitemsPerLine = 5;
	 int nitems = 0;
	 for (int k=0; k<nz; k++)
	 {
	     for (int j=0; j<ny; j++)
	     {
		 for (int i=0; i<nx; i++)
		 {
		     os << var(i,j,k);
		     if (++nitems % nitemsPerLine)
			 os << " ";
		     else
			 os << endl;
		 }
	     }
	 }

	 if (nitems % nitemsPerLine)
	     os << endl;

	 // restore the original flags.
	 
	 os.flags(fmtflags);
     }
 };

}

std::ostream &operator<<(std::ostream &os, const testFullP13T::MT::ccsf &rhs)
{
    typedef testFullP13T::MT::ccsf FT;

    std::ios_base::fmtflags fmtflags = os.flags();

    os << std::setw(16) << std::scientific << std::setprecision(6);
    
#if 0
    int iline = 0;
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
    {
	os << *it << " ";
	if (++iline % 6 == 0)
	    os << endl;
    }
    if (iline % 6 != 0)
	os << endl;
#else
    os << endl;
    int icell = 0;
    for (int k=0; k<nz; k++)
	for (int j=0; j<ny; j++)
	    for (int i=0; i<nx; i++)
	    {
		os << icell++ << ":";
		os << " " << rhs(i,j,k);
		os << endl;
	    }
#endif    

    // restore the original flags.
    
    os.flags(fmtflags);

    return os;
}

std::ostream &operator<<(std::ostream &os, const testFullP13T::MT::fcdsf &rhs)
{
    typedef testFullP13T::MT::fcdsf FT;
    
    std::ios_base::fmtflags fmtflags = os.flags();

    os << std::setw(16) << std::scientific << std::setprecision(6);

#if 0
    int iline = 0;
    for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
    {
	os << *it << " ";
	if (++iline % 6 == 0)
	    os << endl;
    }
    if (iline % 6 != 0)
	os << endl;
#else
    os << endl;

    int icell = 0;
    for (int k=0; k<nz; k++)
	for (int j=0; j<ny; j++)
	    for (int i=0; i<nx; i++)
	    {
		os << icell++ << ":";
		for (int f=0; f<6; f++)
		    os << " " << rhs(i,j,k,f);
		os << endl;
	    }
#endif

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

    spMesh = new MT(mdb);
    
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

    getMatProp();

    P13TOptions options(pdb.P1TauMultiplier, pdb.IsCoupledMaterial);

    spTsManager = new rtt_timestep::ts_manager();
	
    spP13T = new P13T(options, spMesh, spTsManager);

    nx = spMesh->get_ncx();
    ny = spMesh->get_ncy();
    nz = spMesh->get_ncz();
}

testFullP13T::~testFullP13T()
{
    // empty
}

void testFullP13T::getMatProp()
{
    getMatProp(spMatProp);
}

void testFullP13T::getMatProp(SP<MarshakMaterialProps> &spMatProp_) const
{
    spMatProp_ = new MarshakMaterialProps(units, pdb.kappa0, pdb.abar,
					  pdb.kappaPower);
}

void testFullP13T::getMatProp(SP<InterpedMaterialProps> &spMatProp_) const
{
    std::ifstream ifs(pdb.opacityFile);

    using rtt_matprops::FifiMatPropsReader;

    typedef FifiMatPropsReader::MaterialDefinition MatDef;
    vector<MatDef> matdefs;
    matdefs.push_back(MatDef(std::string(pdb.materialName),
			     pdb.materialId, pdb.abar));
    
    FifiMatPropsReader reader(matdefs, units, ifs);

    vector<int> matIds(matdefs.size());
    for (int i=0; i<matdefs.size(); i++)
	matIds[i] = matdefs[i].matid;
    
    spMatProp_ = new InterpedMaterialProps(matIds, reader);
}

testFullP13T::MatStateCC testFullP13T::getMatStateCC(const ccsf &TElect,
						     const ccsf &TIon,
						     const ccsf &density,
						     const ccif &matid) const
{
    return spMatProp->getMaterialState(density, TElect, TIon,  matid);
}

testFullP13T::MatStateFC testFullP13T::getMatStateFC(const ccsf &TElect,
						     const ccsf &TIon,
						     const ccsf &density,
						     const ccif &matid) const
{
    fcdsf TElectFC(spMesh);
    fcdsf TIonFC(spMesh);
    fcdsf densityFC(spMesh);
    fcdif matidFC(spMesh);

    ccsf  oneCC(spMesh);
    oneCC = 1.0;
    
    fcdsf twoFC(spMesh);
    MT::scatter(twoFC, oneCC, MT::OpAddAssign());

    MT::scatter(TElectFC, TElect, MT::OpAddAssign());
    TElectFC /= twoFC;

    MT::scatter(TIonFC, TIon, MT::OpAddAssign());
    TIonFC /= twoFC;

    MT::gather(densityFC, density, MT::OpAssign());
    MT::gather(matidFC, matid, MT::OpAssign());
    
    if (pdb.verbose)
    {
	cout << "In testFullP13T::getMatStateFC" << endl;
	cout << endl;
	cout << "twoFC: " << twoFC << endl;
	cout << "TElectFC: " << TElectFC << endl;
	cout << "TElect: " << TElect << endl;
	cout << endl;
    }

    return spMatProp->getMaterialState(densityFC, TElectFC, TIonFC, matidFC);
}

void testFullP13T::gmvDump(const RadiationStateField &radState,
			   const ccsf &TElec, const ccsf &TIon,
			   int cycle, double time) const
{
    std::ostrstream oss;
    oss << "testFullP13T.gmvout."
	<< std::setw(3) << std::setfill('0') << cycle
	<< std::ends;
    std::ofstream ofs(oss.str());

    using PhysicalConstants::pi;
    
    const RadiationPhysics radphys(units);
    const double a = radphys.getRadConstant();
    const double c = radphys.getLightSpeed();
    
    ccsf TRad(spMesh);

    ccsf::iterator trit = TRad.begin();
    for (ccsf::const_iterator pit = radState.phi.begin();
	 pit != radState.phi.end(); pit++)
    {
	*trit++ = std::pow((*pit) / (a*c), 0.25);
    }

    GmvDump gmv(ofs, spMesh, cycle, time);
    gmv.dump(radState.phi, "phi");
    gmv.dump(TRad, "TRad");
    gmv.dump(TElec, "TElec");
    gmv.dump(TIon, "TIon");
}
	
void testFullP13T::run() const
{
    ccsf TElect0(spMesh);
    ccsf TIon0(spMesh);
    ccsf density(spMesh);
    ccif matid(spMesh);
    // ccsf matid(spMesh);

    TElect0 = pdb.Te;
    TIon0 = pdb.Ti;
    density = pdb.rho;
    matid = pdb.materialId;

    MatStateCC matStateCC = getMatStateCC(TElect0, TIon0, density, matid);
    MatStateFC matStateFC = getMatStateFC(TElect0, TIon0, density, matid);
   
    RadiationStateField radState(spMesh);

    spP13T->initializeRadiationState(matStateCC, radState);

    ccsf QRad(spMesh);
    ccsf QElectron(spMesh);
    ccsf QIon(spMesh);
    bssf boundary(spMesh);
    ccsf electEnergyDep(spMesh);
    ccsf ionEnergyDep(spMesh);

    if (pdb.Qloc < 0)
    {
	QRad = pdb.Qr;
	QElectron = pdb.Qe;
	QIon = pdb.Qi;
    }
    else
    {
	QRad(pdb.Qloc) = pdb.Qr;
	QElectron(pdb.Qloc) = pdb.Qe;
	QIon(pdb.Qloc) = pdb.Qi;
    }

    setBoundary(boundary);
    
    cerr << "Made it before diffSolver ctor" << endl;
    
    spDiffSolver = new DS(diffdb, spMesh, pcg_db);
    
    cerr << "Made it after diffSolver ctor" << endl;

    using dsxx::SP;
    using rtt_timestep::ts_advisor;
    using rtt_timestep::fixed_ts_advisor;
    using rtt_timestep::ratio_ts_advisor;
    
    // Set up a min timestep

    SP<fixed_ts_advisor> spTsMin =
	new fixed_ts_advisor("Minimum",
			     ts_advisor::min, 
			     ts_advisor::small());
    spTsManager->add_advisor(spTsMin);
    spTsMin->set_fixed_value(pdb.dtMin);

    // Set up a max timestep

    SP<fixed_ts_advisor> spTsMax =
	new fixed_ts_advisor("Maximum",
			     ts_advisor::max, 
			     ts_advisor::small());
    spTsManager->add_advisor(spTsMax);
    spTsMax->set_fixed_value(pdb.dtMax);

    // Set up a lower limit on the timestep rate of change

    SP<ratio_ts_advisor> spTsLowRateChange = 
	new ratio_ts_advisor("Rate of Change Lower Limit",
			     ts_advisor::min, 0.8);
    spTsManager->add_advisor(spTsLowRateChange);

    // Set up an upper limit on the time-step rate of change

    SP<ratio_ts_advisor>  spTsHighRateChange =
	new ratio_ts_advisor("Rate of Change Upper Limit");
    spTsManager->add_advisor(spTsHighRateChange);

    const int ncycles = pdb.nsteps;
    double dt = pdb.dt;
    double time = 0.0;
    
    for (int cycle = 1; cycle <= ncycles; cycle++)
    {
	timestep(time, dt, cycle, matStateCC, matStateFC, radState,
		 electEnergyDep, ionEnergyDep,
		 QRad, QElectron, QIon, boundary);
    }
}

void testFullP13T::timestep(double &time, double &dt, int &cycle,
			    MatStateCC &matStateCC,
			    MatStateFC &matStateFC,
			    RadiationStateField &radState,
			    ccsf &electEnergyDep, ccsf &ionEnergyDep,
			    const ccsf &QRad, const ccsf &QElectron,
			    const ccsf &QIon, const bssf &boundary) const
{
    // end of cycle time
    
    time += dt;

    // advance the timestep manager
	
    spTsManager->set_cycle_data(dt, cycle, time);
    
    ccsf TElec(spMesh);
    ccsf TIon(spMesh);
    ccsf CvElec(spMesh);
    ccsf CvIon(spMesh);
    ccsf sigAbs(spMesh);
    ccsf alpha(spMesh);
    ccsf kappaElec(spMesh);
    ccsf kappaIon(spMesh);
    
    matStateCC.getElectronTemperature(TElec);
    matStateCC.getIonTemperature(TIon);
    matStateCC.getElectronSpecificHeat(CvElec);
    matStateCC.getIonSpecificHeat(CvIon);
    matStateCC.getSigmaAbsorption(1,sigAbs);
    matStateCC.getElectronIonCoupling(alpha);
    matStateCC.getElectronConductionCoeff(kappaElec);
    matStateCC.getIonConductionCoeff(kappaIon);

    if (pdb.verbose)
    {
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
    }
    else
    {
	double (*pabs)(double) = std::abs;
	std::pointer_to_unary_function<double,double> opAbs =
	    std::ptr_fun(pabs);
	    
	cout << "TElec (min), (max): " << min(TElec, opAbs) << " " <<
	    max(TElec, opAbs) << endl;
	cout << "TIon (min), (max): " << min(TIon, opAbs) << " " <<
	    max(TIon, opAbs) << endl;
	cout << "CvElec (min), (max): " << min(CvElec, opAbs) << " " <<
	    max(CvElec, opAbs) << endl;
	cout << "CvIon (min), (max): " << min(CvIon, opAbs) << " " <<
	    max(CvIon, opAbs) << endl;
	cout << "sigAbs (min), (max): " << min(sigAbs, opAbs) << " " <<
	    max(sigAbs, opAbs) << endl;
	cout << "alpha (min), (max): " << min(alpha, opAbs) << " " <<
	    max(alpha, opAbs) << endl;
	cout << "kappaElec (min), (max): "
	     << min(kappaElec, opAbs) << " "
	     << max(kappaElec, opAbs) << endl;
	cout << "kappaIon (min), (max): " << min(kappaIon, opAbs)
	     << " " << max(kappaIon, opAbs) << endl;
	cout << "radState.phi (min), (max): " << min(radState.phi, opAbs)
	     << " " << max(radState.phi, opAbs) << endl;
	cout << "radState.F (min), (max): " << min(radState.F, opAbs)
	     << " " << max(radState.F, opAbs) << endl;
    }
    
    cerr << "Made it before solve3T" << endl;
    
    RadiationStateField newRadState(spMesh);
    ccsf QEEM(spMesh);
    ccsf REEM(spMesh);
	
    spP13T->solve3T(dt, matStateCC, matStateFC, radState, QRad,
		    QElectron, QIon,
		    boundary, *spDiffSolver, newRadState, QEEM, REEM,
		    electEnergyDep, ionEnergyDep, TElec, TIon);

    cerr << "Made it after solve3T" << endl;

    gmvDump(newRadState, TElec, TIon, cycle, time);

    ccsf density(spMesh);
    matStateCC.getDensity(density);

    ccif matid(spMesh);
    matStateCC.getMatId(matid);
    
    MatStateCC newMatStateCC = getMatStateCC(TElec, TIon, density, matid);
    MatStateFC newMatStateFC = getMatStateFC(TElec, TIon, density, matid);
	
    postProcess(radState, newRadState, matStateCC, newMatStateCC,
		electEnergyDep, ionEnergyDep, QRad, QElectron, QIon,
		QEEM, REEM, dt);

    dt = spTsManager->compute_new_timestep();
    spTsManager->print_summary();

    matStateCC = newMatStateCC;
    matStateFC = newMatStateFC;
    radState = newRadState;
}

void testFullP13T::setBoundary(bssf &boundary) const
{
    for (int j=0; j<ny; j++)
    {
	for (int k=0; k<nz; k++)
	{
	    boundary(0   , j, k, 0) = pdb.src_left;
	    boundary(nx-1, j, k, 1) = pdb.src_right;
	}
    }

    for (int i=0; i<nx; i++)
    {
	for (int k=0; k<nz; k++)
	{
	    boundary(i, 0   , k, 2) = pdb.src_front;
	    boundary(i, ny-1, k, 3) = pdb.src_back;
	}
    }
    
    for (int i=0; i<nx; i++)
    {
	for (int j=0; j<ny; j++)
	{
	    boundary(i, j, 0   , 4) = pdb.src_bottom;
	    boundary(i, j, nz-1, 5) = pdb.src_top;
	}
    }
}    

void testFullP13T::postProcess(const RadiationStateField &radState,
			       const RadiationStateField &newRadState,
			       const MatStateCC &matStateCC,
			       const MatStateCC &newMatStateCC,
			       const ccsf &electEnergyDepCC,
			       const ccsf &ionEnergyDepCC,
			       const ccsf &QRad,
			       const ccsf &QElectron,
			       const ccsf &QIon,
			       const ccsf &QEEM,
			       const ccsf &REEM,
			       double dt) const
{
    std::ios_base::fmtflags oldOptions = cout.flags(std::ios_base::scientific);
    
    ccsf TElec0(spMesh);
    ccsf TIon0(spMesh);
    matStateCC.getElectronTemperature(TElec0);
    matStateCC.getIonTemperature(TIon0);

    ccsf TElec(spMesh);
    ccsf TIon(spMesh);
    newMatStateCC.getElectronTemperature(TElec);
    newMatStateCC.getIonTemperature(TIon);

    if (pdb.verbose)
    {
	cout << "newRadState.phi: " << newRadState.phi << endl;
	cout << "newRadState.F: " << newRadState.F << endl;
	cout << "TElectron: " << TElec << endl;
	cout << "TIon: " << TIon << endl;
    }
    else
    {
	double (*pabs)(double) = std::abs;
	std::pointer_to_unary_function<double,double> opAbs =
	    std::ptr_fun(pabs);

	cout << "TElec0 (min), (max): " << min(TElec0, opAbs) << " " <<
	    max(TElec0, opAbs) << endl;
	cout << "TElec (min), (max): " << min(TElec, opAbs) << " " <<
	    max(TElec, opAbs) << endl;
	cout << "TIon0 (min), (max): " << min(TIon0, opAbs) << " " <<
	    max(TIon0, opAbs) << endl;
	cout << "TIon (min), (max): " << min(TIon, opAbs) << " " <<
	    max(TIon, opAbs) << endl;
	cout << "radState.phi (min), (max): "
	     << min(radState.phi, opAbs) << " " <<
	    max(radState.phi, opAbs) << endl;
	cout << "newRadState.phi (min), (max): "
	     << min(newRadState.phi, opAbs) << " "
	     << max(newRadState.phi, opAbs) << endl;
    }

    const double dx = spMesh->get_dx();
    const double dy = spMesh->get_dy();
    const double dz = spMesh->get_dz();
    const double nx = spMesh->get_ncx();
    const double ny = spMesh->get_ncy();
    const double nz = spMesh->get_ncz();
    
    double volpcell = dx*dy*dz;
    double volume = volpcell*nx*ny*nz;

    cout << endl;
    cout << "volume: " << volume << endl;
    
    const RadiationPhysics radPhys(spMatProp->getUnits());
    const double c = radPhys.getLightSpeed();
    
    ccsf deltaRadEnergyCC(spMesh);
    deltaRadEnergyCC = (newRadState.phi - radState.phi) / c;

    const double deltaRadEnergy = sum(deltaRadEnergyCC)*volpcell;
    cout << "deltaRadEnergy: " << deltaRadEnergy
	 << "\trate: " << deltaRadEnergy/dt << endl;

    const double electEnergyDep = sum(electEnergyDepCC)*volpcell;
    cout << "electEnergyDep: " << electEnergyDep
	 << "\trate: " << electEnergyDep/dt << endl;

    ccsf CvElec(spMesh);
    ccsf CvIon(spMesh);
    matStateCC.getElectronSpecificHeat(CvElec);
    matStateCC.getIonSpecificHeat(CvIon);
    
    // ccsf temp(spMesh);
    // temp = (TElec - TElec0)*CvElec;
    // delta = sum(temp)*volpcell;
    // cout << "electEnergyDep (recalc): " << delta
    //      << "\trate: " << delta/dt << endl;

    const double ionEnergyDep = sum(ionEnergyDepCC)*volpcell;
    cout << "  ionEnergyDep: " << ionEnergyDep
	 << "\trate: " << ionEnergyDep/dt << endl;

    const double energyDep = deltaRadEnergy + electEnergyDep + ionEnergyDep;

    cout << "energyDep: " << energyDep
	 << "\trate: " << energyDep/dt << endl;

    const double inhomosrc = sum(QRad) * volpcell;
    const double qsrc = inhomosrc + (sum(QElectron) + sum(QIon))*volpcell;

    cout << "Volume src: " << qsrc << endl;

    double bndsrc;
    bndsrc  = (pdb.src_left + pdb.src_right) * ny*dy * nz*dz ;
    bndsrc += (pdb.src_front + pdb.src_back) * nz*dz * nx*dx;
    bndsrc += (pdb.src_bottom + pdb.src_top) * nx*dx * ny*dy;

    cout << "Boundary src: " << bndsrc << endl;

    const double externsrc = bndsrc + qsrc;
    cout << "Total external src: " << externsrc << endl;

    const double timedepsrc = sum(radState.phi) * volpcell / (c*dt);
    cout << "time dep src: " << timedepsrc << endl;
    
    const double timedeprem = sum(newRadState.phi) * volpcell / (c*dt);
    cout << "time dep removal: " << timedeprem << endl;
    
    ccsf sigAbs(spMesh);
    matStateCC.getSigmaAbsorption(1,sigAbs);

    ccsf temp(spMesh);
    temp = newRadState.phi*sigAbs;
    const double absorption = sum(temp) * volpcell;
    cout << "absorption: " << absorption << endl;

    temp = newRadState.phi*REEM;
    const double emisrem = sum(temp) * volpcell;
    cout << "emissive removal: " << emisrem << endl;
    
    const double emission = sum(QEEM) *	volpcell;
    cout << "emission: " << emission << endl;

    // cout << "emission-absorption: " << emission - absorption << endl;

    const double leakage = calcLeakage(newRadState);
    
    cout << "leakage (x-y): " << leakage << endl;

    const double totsrc = inhomosrc + bndsrc + timedepsrc + emission;
    cout << "total radiation sources: " << totsrc << endl;

    const double totrem = leakage + absorption + emisrem + timedeprem;
    cout << "total radiation removal: " << totrem << endl;
    
    const double balance = totsrc - totrem;
    cout << "Rad only Balance: " << balance << endl;

    const double relbal = balance / totsrc;
    cout << "Relative Balance: " << relbal << endl;

#if 0
    std::ofstream ofs("testFullP13T.dat");
    for (int m=0; m<nz; m++)
	ofs << m << '\t' << newRadState.phi(0,0,m)
	    << '\t' << newRadState.F(0,0,m,4)
	    << '\t' << newRadState.F(0,0,m,5)
	    << endl;
#endif

    cout << endl;
    cout.flags(oldOptions);
}

double testFullP13T::calcLeakage(const RadiationStateField &radstate) const
{
    const double dx = spMesh->get_dx();
    const double dy = spMesh->get_dy();
    const double dz = spMesh->get_dz();
    const double nx = spMesh->get_ncx();
    const double ny = spMesh->get_ncy();
    const double nz = spMesh->get_ncz();

    double leakage = 0.0;

    for (int k=0; k<nx; k++)
    {
	for (int l=0; l<ny; l++)
	{
	    if (diffdb.alpha_bottom != 0)
	    {
		double alpha_bottom = diffdb.alpha_bottom;
		double beta_bottom = diffdb.beta_bottom;
		double src_bottom = pdb.src_bottom;

		double F_bottom = radstate.F(k,l,0,4);
		
		// calculate phi along bottom face
		double phi_bottom = (src_bottom - beta_bottom*F_bottom) /
		    alpha_bottom;
		
		leakage += phi_bottom/4.0 + F_bottom/2.0;
		// leakage += alpha_bottom*phi_bottom - beta_bottom*F_bottom;
	    }

	    if (diffdb.alpha_top != 0)
	    {
		double alpha_top = diffdb.alpha_top;
		double beta_top = diffdb.beta_top;
		double src_top = pdb.src_top;

		double F_top = radstate.F(k,l,nz-1,5);
		
		// calculate phi along top face
		double phi_top = (src_top - beta_top*F_top) / alpha_top;

		leakage += phi_top/4.0 + F_top/2.0;
		// leakage += alpha_top*phi_top - beta_top*F_top;
	    }
	}
    }
    leakage *= dx*dy;

    return leakage;
}

//---------------------------------------------------------------------------//
//                              end of testFullP13T.cc
//---------------------------------------------------------------------------//
