//----------------------------------*-C++-*----------------------------------//
// testFullP13T.cc
// Randy M. Roberts
// Sat May 30 15:09:58 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testFullP13T.hh"

#include "matprops/FifiMatPropsReader.hh"
#include "matprops/TempMapper.hh"
#include "3T/P13TOptions.hh"
#include "units/Units.hh"
#include "radphys/RadiationPhysics.hh"
#include "nml/Group.hh"
#include "timestep/ts_manager.hh"
#include "timestep/fixed_ts_advisor.hh"
#include "timestep/ratio_ts_advisor.hh"
#include "matprops/TempMapper.hh"
#include "3T/testP13T/GmvDump.hh"

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

 enum Faces { LEFT=0, RIGHT=1, FRONT=2, BACK=3, BOTTOM=4, TOP=5 };
 
 int nx;
 int ny;
 int nz;


 template<class FT1, class FT>
 void cpCell2Mats(FT1 &matValues, const FT &cellValue, const int nmat)
 {
     int ncell = cellValue.size();
     FT1::iterator mvit = matValues.begin();
     FT::const_iterator cit = cellValue.begin();
     for (int icell = 0; icell < ncell; icell++, cit++, mvit++)
     {
	 FT1::value_type tmp(nmat);
	 for (FT1::value_type::iterator tmpit = tmp.begin(); 
	      tmpit != tmp.end(); tmpit++)
	 {
	     *tmpit = *cit;
	 }
	 *mvit = tmp;
     }
 }

 template<class FT1>
 void assign2Mats(FT1 &matValues, const int ncell, const int nmat, 
		  const int value) 
 {
     FT::const_iterator cit = cellValue.begin();
     FT1::iterator mvit = matValues.begin();
     for (int icell = 0; icell < ncell; icell++, cit++, mvit++)
     {
	 FT1::value_type tmp(nmat);
	 for (FT1::value_type::iterator tmpit = tmp.begin(); 
	      tmpit != tmp.end(); tmpit++)
	 {
	     *tmpit = value;
	 }
	 *mvit = tmp;
     }
 }


 template<class UMCMP>
 bool isContinuous(const typename testFullP13T<UMCMP>::MT::fcdsf &rhs)
 {
     for (int i=0; i<nx; i++)
	 for (int j=0; j<ny; j++)
	     for (int k=0; k<nz; k++)
	     {
		 if (i > 0 && (rhs(i,j,k,LEFT) != -rhs(i-1,j,k,RIGHT)))
		     return false;
		 if (j > 0 && (rhs(i,j,k,FRONT) != -rhs(i,j-1,k,BACK)))
		     return false;
		 if (k > 0 && (rhs(i,j,k,BOTTOM) != -rhs(i,j,k-1,TOP)))
		     return false;
	     }
     return true;
 }
 
 template<class UMCMP>
 double sum(const typename testFullP13T<UMCMP>::MT::ccsf &rhs)
 {
     typedef testFullP13T<UMCMP>::MT::ccsf FT;

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

 std::ostream &operator<<(std::ostream &os, 
			  const Mesh_XYZ::ccsf &rhs)
 {
     typedef Mesh_XYZ::ccsf FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
    
#if 0
     int iline = 0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
     {
	 os << std::setw(16) << *it << " ";
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
		 os << std::setw(5) << icell++ << ":";
		 os << " " << std::setw(16) << rhs(i,j,k);
		 os << endl;
	     }
#endif    

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 std::ostream &operator<<(std::ostream &os, 
			  const Mesh_XYZ::cctf<std::vector<double> > &rhs)
 {
     typedef Mesh_XYZ::cctf<std::vector<double> > FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
    
#if 0
     int iline = 0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
     {
	 os << std::setw(16) << *it << " ";
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
		 os << std::setw(5) << icell++ << ":";
		 int nmat = rhs(i,j,k).size();
		 for (int imat = 0; imat < nmat; imat++)
		 {
		     os << " " << std::setw(16) << rhs(i,j,k)[imat];
		 }
		 os << endl;
	     }
#endif    

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }


 std::ostream &operator<<(std::ostream &os,
	     const Mesh_XYZ::fcdsf &rhs)
 {
     typedef Mesh_XYZ::fcdsf FT;
    
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);

#if 0
     int iline = 0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
     {
	 os << std::setw(16) << *it << " ";
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
		 os << std::setw(5) << icell++ << ":";
		 for (int f=0; f<6; f++)
		     os << " " << std::setw(16) << rhs(i,j,k,f);
		 os << endl;
	     }
#endif

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 std::ostream &operator<<(std::ostream &os,
			  const Mesh_XYZ::fcdtf<std::vector<double> > &rhs)
 {
     typedef Mesh_XYZ::fcdtf<std::vector<double> > FT;
    
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);

#if 0
     int iline = 0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
     {
	 os << std::setw(16) << *it << " ";
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
		 for (int f=0; f<6; f++)
		 {
		     os << std::setw(5) << icell++ << ":" << f << ":";
		     int nmat = rhs(i,j,k,f).size();
		     for (int imat=0; imat<nmat; imat++)
			 os << " " << std::setw(16) << rhs(i,j,k,f)[imat];
		     os << endl;
		 }
	     }
#endif

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 using rtt_matprops::MultiMatCellMatProps;

 template<class UMCMP>
 void getMatProp(dsxx::SP<UMCMP> &spUMatProp,
		 dsxx::SP<MultiMatCellMatProps<UMCMP> > &spMatProp,
		 testFullP13T<UMCMP> &tester,
		 const testFullP13T_DB &tdb)
 {
     Assert(0);
 }

 using rtt_matprops::MarshakMaterialProps;

 template<>
 void
 getMatProp(dsxx::SP<MarshakMaterialProps> &spUMatProp,
	    dsxx::SP<MultiMatCellMatProps<MarshakMaterialProps> > &spMatProp,
	    testFullP13T<MarshakMaterialProps> &tester,
	    const testFullP13T_DB &tdb) 
 {
     spUMatProp = new MarshakMaterialProps(tester.getUnits(), tdb.kappa0,
					   tdb.abar, tdb.kappaPower);

     spMatProp = new MultiMatCellMatProps<MarshakMaterialProps>(spUMatProp);
 }

 using rtt_matprops::InterpedMaterialProps;
    
 template<>
 void
 getMatProp(dsxx::SP<InterpedMaterialProps> &spUMatProp,
	    dsxx::SP<MultiMatCellMatProps<InterpedMaterialProps> > &spMatProp,
	    testFullP13T<InterpedMaterialProps> &tester,
	    const testFullP13T_DB &tdb)
 {
     std::ifstream ifs(tdb.opacityFile);

     using rtt_matprops::FifiMatPropsReader;

     typedef FifiMatPropsReader::MaterialDefinition MatDef;
     vector<MatDef> matdefs;
     matdefs.push_back(MatDef(std::string(tdb.materialName),
			      tdb.materialId, tdb.abar));
    
     FifiMatPropsReader reader(matdefs, tester.getUnits(), ifs);

     vector<int> matIds(matdefs.size());
     for (int i=0; i<matdefs.size(); i++)
	 matIds[i] = matdefs[i].matid;

     spUMatProp = new InterpedMaterialProps(matIds, reader);

     spMatProp = new MultiMatCellMatProps<InterpedMaterialProps>(spUMatProp);
 }

 void setTempFromFile(Mesh_XYZ::ccsf &Temp, const std::string &filename, double floor)
 {
     std::ifstream ifs(filename.c_str());

     ifs.ignore(10000, '\n');
     ifs.ignore(10000, '\n');

     for (int k=0; k<nz; k++)
     {
	 double r, T1, T2, T3;

	 ifs >> r >> T1 >> T2 >> T3;
	 if (ifs.eof())
	     throw std::runtime_error("Premature EOF reading temperatures.");

	 if (T1 < floor)
	     T1 = floor;
	 
	 for (int i=0; i<nx; i++)
	     for (int j=0; j<ny; j++)
		 Temp(i, j, k) = T1;
     }     
 }
 
} // end unnamed namespace

template<class UMCMP>
testFullP13T<UMCMP>::testFullP13T(const testFullP13T_DB &tdb_,
				  const Diffusion_DB &diffdb_,
				  const Mesh_DB &mdb,
				  const pcg_DB &pcg_db_)
    : diffdb(diffdb_), pcg_db(pcg_db_), tdb(tdb_)
{
    spMesh = new MT(mdb);
    
    switch (tdb.units)
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

    P13TOptions options(tdb.P1TauMultiplier, tdb.IsCoupledMaterial);

    spTsManager = new rtt_timestep::ts_manager();

    // Set up a informational advisor to
    // contain the current time-step for reference.
    // Activating this controller can also be used to
    // freeze the time-step at the current value.

    spTsCurrent =
	new rtt_timestep::fixed_ts_advisor("Current Time-Step",
					   rtt_timestep::ts_advisor::req,
					   tdb.dt, false);
    spTsManager->add_advisor(spTsCurrent);

	
    spP13T = new P13T(options, spMesh, spTsManager);

    nx = spMesh->get_ncx();
    ny = spMesh->get_ncy();
    nz = spMesh->get_ncz();

    // gamma is the high weight when computing the weighted mean of
    // cc temperatures to construct fc temperatures.
    
    double gamma = tdb.faceGamma;
    spTempMapper = new rtt_matprops::TempMapper<MT>(spMesh, gamma);
}

template<class UMCMP>
testFullP13T<UMCMP>::~testFullP13T()
{
    // empty
}

template<class UMCMP>
void testFullP13T<UMCMP>::getMatProp()
{
    ::getMatProp(spUMatProp, spMatProp, *this, tdb);
}

template<class UMCMP>
testFullP13T<UMCMP>::MatStateCC
testFullP13T<UMCMP>::getMatStateCC(const ccvsf &TElect,
				   const ccvsf &TIon,
				   const ccvsf &density,
				   const ccvsf &VolFrac,
				   const ccvif &matid) const
{
    return spMatProp->getMaterialState<ccsf, ccvsf, ccvif>
	(density, TElect, TIon, VolFrac, matid);
}

template<class UMCMP>
testFullP13T<UMCMP>::MatStateFC 
testFullP13T<UMCMP>::getMatStateFC(const MatStateCC &msfcc) const
{
    ccsf avgTemp(spMesh);
    ccvsf matTemps(spMesh);
    msfcc.getElectronTemperature(avgTemp);
    msfcc.getElectronTempByMat(matTemps);

    fcdvsf TElectFC(spMesh);
    spTempMapper->tempCC2FC(TElectFC, matTemps, avgTemp, MT::OpAssign());

    msfcc.getIonTemperature(avgTemp);
    msfcc.getIonTempByMat(matTemps);

    fcdvsf TIonFC(spMesh);
    spTempMapper->tempCC2FC(TIonFC, matTemps, avgTemp, MT::OpAssign());

    ccvsf density(spMesh);
    ccvsf volFrac(spMesh);
    ccvif matid(spMesh);

    msfcc.getDensity(density);
    msfcc.getVolumeFraction(volFrac);
    msfcc.getMatId(matid);

    fcdvsf densityFC(spMesh);
    fcdvsf volFracFC(spMesh);
    fcdvif matidFC(spMesh);

    MT::gather(densityFC, density, MT::OpAssign());
    MT::gather(volFracFC, volFrac, MT::OpAssign());
    MT::gather(matidFC, matid, MT::OpAssign());
    
    if (tdb.Te_bottom > 0.)
    {
	for (int i=0; i<nx; i++)
	    for (int j=0; j<ny; j++)
	    {
		using std::fill;
		fill(TElectFC(i, j, 0, BOTTOM).begin(),
		     TElectFC(i, j, 0, BOTTOM).end(), tdb.Te_bottom);
		fill(TIonFC(i, j, 0, BOTTOM).begin(),
		     TIonFC(i, j, 0, BOTTOM).end(), tdb.Te_bottom);
	    }
    }

    if (tdb.verbose)
    {
	cout << "In testFullP13T::getMatStateFC" << endl;
	cout << endl;
	cout << "TElectFC: " << TElectFC << endl;
	cout << endl;
    }

    return spMatProp->getMaterialState<fcdsf, fcdvsf, fcdvif>(
	densityFC, TElectFC, TIonFC,
	volFracFC, matidFC);
}

template<class UMCMP>
void testFullP13T<UMCMP>::gmvDump(const RadiationStateField &radState,
				  const ccsf &TElec, const ccsf &TIon,
				  int dumpno, int cycle, double time) const
{
    std::ostrstream oss;
    oss << "testFullP13T.gmvout."
	<< std::setw(5) << std::setfill('0') << dumpno
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

    rtt_3T_testP13T::GmvDump<MT> gmv(ofs, spMesh, cycle, time);
    gmv.dump(radState.phi, "phi");
    gmv.dump(TRad, "TRad");
    gmv.dump(TElec, "TElec");
    gmv.dump(TIon, "TIon");
}
	
template<class UMCMP>
void testFullP13T<UMCMP>::run() const
{
    ccsf TElect0_in(spMesh);
    ccsf TIon0_in(spMesh);
    ccsf density_in(spMesh);
    ccsf VolFrac_in(spMesh);
    ccif matid_in(spMesh);

    std::string TFile = static_cast<const char*>(tdb.TFile);

    if (TFile != "")
    {
	setTempFromFile(TElect0_in, TFile, tdb.Te);
	TIon0_in = TElect0_in;
    }
    else
    {
	TElect0_in = tdb.Te;
	TIon0_in = tdb.Ti;
    }
    
    density_in = tdb.rho;
    VolFrac_in = 1.;
    matid_in = tdb.materialId;

    ccvsf TElect0(spMesh);
    ccvsf TIon0(spMesh);
    ccvsf density(spMesh);
    ccvsf VolFrac(spMesh);
    ccvif matid(spMesh); 

    int nmat_prob = 1;
    cpCell2Mats<ccvsf,ccsf>(TElect0, TElect0_in, nmat_prob);
    cpCell2Mats<ccvsf,ccsf>(TIon0, TIon0_in, nmat_prob);
    cpCell2Mats<ccvsf,ccsf>(density, density_in, nmat_prob);
    cpCell2Mats<ccvsf,ccsf>(VolFrac, VolFrac_in, nmat_prob);
    cpCell2Mats<ccvif,ccif>(matid, matid_in, nmat_prob);

    MatStateCC matStateCC = getMatStateCC(TElect0, TIon0, density, 
					  VolFrac, matid);
    MatStateFC matStateFC = getMatStateFC(matStateCC);
   
    RadiationStateField radState(spMesh);

    spP13T->initializeRadiationState(matStateCC, radState);

    ccsf QRad(spMesh);
    ccsf QElectron(spMesh);
    ccsf QIon(spMesh);
    bssf alpha(spMesh);
    bssf beta(spMesh);
    bssf bSrc(spMesh);
    ccsf electEnergyDep(spMesh);
    ccsf ionEnergyDep(spMesh);
    ncvsf momentumDeposition(spMesh);

    if (tdb.Qloc < 0)
    {
	QRad = tdb.Qr;
	QElectron = tdb.Qe;
	QIon = tdb.Qi;
    }
    else
    {
	QRad(tdb.Qloc) = tdb.Qr;
	QElectron(tdb.Qloc) = tdb.Qe;
	QIon(tdb.Qloc) = tdb.Qi;
    }

    setBoundary(alpha, beta, bSrc);
    
    std::cerr << "Made it before diffSolver ctor" << endl;

    SP<MatrixSolver> spMatrixSolver = new MatrixSolver(spMesh, pcg_db);
    spDiffSolver = new DS(diffdb, spMesh, spMatrixSolver);
    
    std::cerr << "Made it after diffSolver ctor" << endl;

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
    spTsMin->set_fixed_value(tdb.dtMin);

    // Set up a max timestep

    SP<fixed_ts_advisor> spTsMax =
	new fixed_ts_advisor("Maximum",
			     ts_advisor::max, 
			     ts_advisor::small());
    spTsManager->add_advisor(spTsMax);
    spTsMax->set_fixed_value(tdb.dtMax);

    // Set up a lower limit on the timestep rate of change

    SP<ratio_ts_advisor> spTsLowRateChange = 
	new ratio_ts_advisor("Rate of Change Lower Limit",
			     ts_advisor::min, 0.8);
    spTsManager->add_advisor(spTsLowRateChange);

    // Set up an upper limit on the time-step rate of change

    SP<ratio_ts_advisor>  spTsHighRateChange =
	new ratio_ts_advisor("Rate of Change Upper Limit");
    spTsManager->add_advisor(spTsHighRateChange);

    const int ncycles = tdb.nsteps;
    double dt = tdb.dt;
    double time = 0.0;
    
    for (int cycle = 1; cycle <= ncycles; cycle++)
    {
	timestep(time, dt, cycle, matStateCC, matStateFC, radState,
		 electEnergyDep, ionEnergyDep, momentumDeposition,
		 QRad, QElectron, QIon, alpha, beta, bSrc);
    }
}

template<class UMCMP>
void testFullP13T<UMCMP>::timestep(double &time, double &dt, int &cycle,
				   MatStateCC &matStateCC,
				   MatStateFC &matStateFC,
				   RadiationStateField &radState,
				   ccsf &electEnergyDep, ccsf &ionEnergyDep,
				   ncvsf &momentumDeposition,
                                   const ccsf &QRad, const ccsf &QElectron,
				   const ccsf &QIon, const bssf &alpha,
				   const bssf &beta, const bssf &bSrc) const
{
    // end of cycle time
    
    time += dt;

    // advance the timestep manager

    spTsCurrent->set_fixed_value(dt);
    spTsManager->set_cycle_data(dt, cycle, time);
    
    ccsf TElec(spMesh);
    ccsf TIon(spMesh);
    ccsf CvElec(spMesh);
    ccsf CvIon(spMesh);
    ccsf sigAbs(spMesh);
    ccsf coupleEI(spMesh);
    ccsf kappaElec(spMesh);
    ccsf kappaIon(spMesh);
    ncvsf velocity(spMesh);
    
    matStateCC.getElectronTemperature(TElec);
    matStateCC.getIonTemperature(TIon);
    matStateCC.getElectronSpecificHeat(CvElec);
    matStateCC.getIonSpecificHeat(CvIon);
    matStateCC.getSigmaAbsorption(1,sigAbs);
    matStateCC.getElectronIonCoupling(coupleEI);
    matStateCC.getElectronConductionCoeff(kappaElec);
    matStateCC.getIonConductionCoeff(kappaIon);
    ncvsf::value_type onevec;
    onevec = 1.;
    velocity = onevec;

    if (tdb.verbose)
    {
	cout << "TElectron: " << TElec << endl;
	cout << "TIon: " << TIon << endl;
	cout << "Cv Electron: " << CvElec << endl;
	cout << "Cv Ion: " << CvIon << endl;
	cout << "Absorption: " << sigAbs << endl;
	cout << "coupleEI: " << coupleEI << endl;
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
	cout << "coupleEI (min), (max): " << min(coupleEI, opAbs) << " " <<
	    max(coupleEI, opAbs) << endl;
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
    
    std::cerr << "Made it before solve3T" << endl;
    
    RadiationStateField newRadState(spMesh);
    ccsf QEEM(spMesh);
    ccsf REEM(spMesh);
	
    spP13T->solve3T(newRadState, QEEM, REEM,
		    electEnergyDep, ionEnergyDep, momentumDeposition, 
                    TElec, TIon,
		    *spDiffSolver, dt, matStateCC, matStateFC, velocity,
                    radState,
		    QRad, QElectron, QIon,
		    alpha, beta, bSrc);

    std::cerr << "Made it after solve3T" << endl;

    Assert(isContinuous<UMCMP>(newRadState.F));

    if (cycle % tdb.dumpcycles == 0)
    {
	static int dumpno = 0;
	gmvDump(newRadState, TElec, TIon, dumpno++, cycle, time);
    }

    ccvsf density(spMesh);
    matStateCC.getDensity(density);

    ccvsf volfrac(spMesh);
    matStateCC.getVolumeFraction(volfrac);

    ccvif matid(spMesh);
    matStateCC.getMatId(matid);

    // Map electron and ion temperatures back to the individual materials
    // within each cell.
    
    ccvsf TIonByMat(spMesh);
    ccvsf TElecByMat(spMesh);
    matStateCC.mapAvgIonTemp(TIonByMat, TIon);
    matStateCC.mapAvgElectronTemp(TElecByMat, TElec);
    
    MatStateCC newMatStateCC = getMatStateCC(TElecByMat, TIonByMat, 
					     density, volfrac, matid);
    MatStateFC newMatStateFC = getMatStateFC(newMatStateCC);
	
#if 1
    if (cycle % tdb.dumpcycles == 0)
    {
	static int dumpno = 0;
    
	std::ostrstream oss;
	oss << "testFullP13T."
	    << std::setw(5) << std::setfill('0') << dumpno++
	    << ".dat"
	    << std::ends;
	std::ofstream ofs(oss.str());
	ofs << "# testFullP13T cycle=" << cycle << " time=" << time << endl;
	ofs << "# z \t phi(z)" << endl;
	ccsf TElect(spMesh);
	newMatStateCC.getElectronTemperature(TElect);
	const double dz = spMesh->get_dz();
	for (int m=0; m<nz; m++)
	    ofs << dz*(m+.5) << '\t' << TElect(0,0,m)
		<< endl;
    }
#endif

    postProcess(radState, newRadState, matStateCC, newMatStateCC,
		electEnergyDep, ionEnergyDep, QRad, QElectron, QIon,
		QEEM, REEM, dt);

    dt = spTsManager->compute_new_timestep();
    spTsManager->print_summary();

    matStateCC = newMatStateCC;
    matStateFC = newMatStateFC;
    radState = newRadState;
}

template<class UMCMP>
void testFullP13T<UMCMP>::setBoundary(bssf &alpha, bssf &beta, bssf &bSrc) const
{
    for (int j=0; j<ny; j++)
    {
	for (int k=0; k<nz; k++)
	{
	    alpha(0,    j, k, LEFT) = diffdb.alpha_left;
	    beta (0   , j, k, LEFT) = diffdb.beta_left;
	    bSrc (0   , j, k, LEFT) = tdb.src_left;
	    alpha(nx-1, j, k, RIGHT) = diffdb.alpha_right;
	    beta (nx-1, j, k, RIGHT) = diffdb.beta_right;
	    bSrc (nx-1, j, k, RIGHT) = tdb.src_right;
	}
    }

    for (int i=0; i<nx; i++)
    {
	for (int k=0; k<nz; k++)
	{
	    alpha(i, 0   , k, FRONT) = diffdb.alpha_front;
	    beta (i, 0   , k, FRONT) = diffdb.beta_front;
	    bSrc (i, 0   , k, FRONT) = tdb.src_front;
	    alpha(i, ny-1, k, BACK) = diffdb.alpha_back;
	    beta (i, ny-1, k, BACK) = diffdb.beta_back;
	    bSrc (i, ny-1, k, BACK) = tdb.src_back;
	}
    }
    
    if (tdb.Te_bottom > 0.0)
    {
	const RadiationPhysics radphys(units);
	double phi_bottom = 0.0;
	radphys.getPlanck(tdb.Te_bottom, phi_bottom);

	phi_bottom *= 4.0*PhysicalConstants::pi;

	for (int i=0; i<nx; i++)
	{
	    for (int j=0; j<ny; j++)
	    {
		alpha(i, j, 0   , BOTTOM) = diffdb.alpha_bottom;
		beta (i, j, 0   , BOTTOM) = diffdb.beta_bottom;
		bSrc (i, j, 0   , BOTTOM) = phi_bottom*diffdb.alpha_bottom;
		alpha(i, j, nz-1, TOP) = diffdb.alpha_top;
		beta (i, j, nz-1, TOP) = diffdb.beta_top;
		bSrc (i, j, nz-1, TOP) = tdb.src_top;
	    }
	}
    }
    else
    {
	for (int i=0; i<nx; i++)
	{
	    for (int j=0; j<ny; j++)
	    {
		alpha(i, j, 0   , BOTTOM) = diffdb.alpha_bottom;
		beta (i, j, 0   , BOTTOM) = diffdb.beta_bottom;
		bSrc (i, j, 0   , BOTTOM) = tdb.src_bottom;
		alpha(i, j, nz-1, TOP) = diffdb.alpha_top;
		beta (i, j, nz-1, TOP) = diffdb.beta_top;
		bSrc (i, j, nz-1, TOP) = tdb.src_top;
	    }
	}
    }
}    

template<class UMCMP>
void testFullP13T<UMCMP>::postProcess(const RadiationStateField &radState,
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

    if (tdb.verbose)
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

    const double deltaRadEnergy = sum<UMCMP>(deltaRadEnergyCC)*volpcell;
    cout << "deltaRadEnergy: " << deltaRadEnergy
	 << "\trate: " << deltaRadEnergy/dt << endl;

    const double electEnergyDep = sum<UMCMP>(electEnergyDepCC)*volpcell;
    cout << "electEnergyDep: " << electEnergyDep
	 << "\trate: " << electEnergyDep/dt << endl;

    ccsf CvElec(spMesh);
    ccsf CvIon(spMesh);
    matStateCC.getElectronSpecificHeat(CvElec);
    matStateCC.getIonSpecificHeat(CvIon);
    
    // ccsf temp(spMesh);
    // temp = (TElec - TElec0)*CvElec;
    // delta = sum<UMCMP>(temp)*volpcell;
    // cout << "electEnergyDep (recalc): " << delta
    //      << "\trate: " << delta/dt << endl;

    const double ionEnergyDep = sum<UMCMP>(ionEnergyDepCC)*volpcell;
    cout << "  ionEnergyDep: " << ionEnergyDep
	 << "\trate: " << ionEnergyDep/dt << endl;

    const double energyDep = deltaRadEnergy + electEnergyDep + ionEnergyDep;

    cout << "energyDep: " << energyDep
	 << "\trate: " << energyDep/dt << endl;

    const double inhomosrc = sum<UMCMP>(QRad) * volpcell;
    const double qsrc = inhomosrc + (sum<UMCMP>(QElectron) +
				     sum<UMCMP>(QIon))*volpcell;

    cout << "Volume src: " << qsrc << endl;

    double bndsrc;
    bndsrc  = (tdb.src_left + tdb.src_right) * ny*dy * nz*dz ;
    bndsrc += (tdb.src_front + tdb.src_back) * nz*dz * nx*dx;
    bndsrc += (tdb.src_bottom + tdb.src_top) * nx*dx * ny*dy;

    cout << "Boundary src: " << bndsrc << endl;

    const double externsrc = bndsrc + qsrc;
    cout << "Total external src: " << externsrc << endl;

    const double timedepsrc = sum<UMCMP>(radState.phi) * volpcell / (c*dt);
    cout << "time dep src: " << timedepsrc << endl;
    
    const double timedeprem = sum<UMCMP>(newRadState.phi) * volpcell / (c*dt);
    cout << "time dep removal: " << timedeprem << endl;
    
    ccsf sigAbs(spMesh);
    matStateCC.getSigmaAbsorption(1,sigAbs);

    ccsf temp(spMesh);
    temp = newRadState.phi*sigAbs;
    const double absorption = sum<UMCMP>(temp) * volpcell;
    cout << "absorption: " << absorption << endl;

    temp = newRadState.phi*REEM;
    const double emisrem = sum<UMCMP>(temp) * volpcell;
    cout << "emissive removal: " << emisrem << endl;
    
    const double emission = sum<UMCMP>(QEEM) *	volpcell;
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

    cout << endl;
    cout.flags(oldOptions);
}

template<class UMCMP>
double
testFullP13T<UMCMP>::calcLeakage(const RadiationStateField &radstate) const
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
		double src_bottom = tdb.src_bottom;

		double F_bottom = radstate.F(k,l,0,BOTTOM);
		
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
		double src_top = tdb.src_top;

		double F_top = radstate.F(k,l,nz-1,TOP);
		
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
