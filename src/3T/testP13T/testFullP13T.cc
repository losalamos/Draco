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
#include "3T/OpCrossSectionMapper.hh"
#include "units/Units.hh"
#include "radphys/RadiationPhysics.hh"
#include "nml/Group.hh"
#include "timestep/ts_manager.hh"
#include "timestep/fixed_ts_advisor.hh"
#include "timestep/ratio_ts_advisor.hh"
#include "matprops/TempMapper.hh"
#include "3T/testP13T/GmvDump.hh"
#include "3T/testP13T/testMaterialProps.hh"
#include "3T/testP13T/utils.hh"
#include "c4/global.hh"

#include <functional>
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

#include "3T/testP13T/Mesh_XYZ_IO.hh"

using dsxx::SP;

using namespace XTM;

namespace
{

 using rtt_3T_testP13T::operator<<;
 
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


    // Create an object that will map cross sections from faces
    // to vertices.
    // You can create your mapper class own as long as it inherits
    // from rtt_3T::CrossSectionMapper<DS>.
    
    SP<const rtt_3T::OpCrossSectionMapper<DS,MT::OpMinAssign> >
	spOpCrossSectionMapper
	= new rtt_3T::OpCrossSectionMapper<DS,MT::OpMinAssign>;
    
    spP13T = new P13T(options, spMesh, spOpCrossSectionMapper, spTsManager);

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
	using rtt_3T_testP13T::BOTTOM;
	using rtt_3T_testP13T::setBoundary;

	MT::bstf<vector<double> > bndTemp(spMesh);

	MT::gather(bndTemp, TElectFC, MT::OpAssign());
	setBoundary(bndTemp, tdb.Te_bottom, BOTTOM);
	MT::gather(TElectFC, bndTemp, MT::OpAssign());

	MT::gather(bndTemp, TIonFC, MT::OpAssign());
	setBoundary(bndTemp, tdb.Te_bottom, BOTTOM);
	MT::gather(TIonFC, bndTemp, MT::OpAssign());
    }

    if (tdb.verbose)
    {
	cout << "TElectFC: " << TElectFC << endl;
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
    std::string fname = rtt_3T_testP13T::getFileName("testFullP13T.gmvout.",
						     "", dumpno);
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

    rtt_3T_testP13T::GmvDump<MT> gmv(fname, spMesh, cycle, time);
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
	rtt_3T_testP13T::setTempFromFile(TElect0_in, TFile, tdb.Te);
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
    std::cerr << C4::node() << " Made it before first cpCell2Mats" << endl;
    cpCell2Mats<ccvsf,ccsf>(TElect0, TElect0_in, nmat_prob);
    cpCell2Mats<ccvsf,ccsf>(TIon0, TIon0_in, nmat_prob);
    cpCell2Mats<ccvsf,ccsf>(density, density_in, nmat_prob);
    cpCell2Mats<ccvsf,ccsf>(VolFrac, VolFrac_in, nmat_prob);
    cpCell2Mats<ccvif,ccif>(matid, matid_in, nmat_prob);
    std::cerr << C4::node() << " Made it after last cpCell2Mats" << endl;

    std::cerr << C4::node() << " Made it before getMatStateCC" << endl;
    MatStateCC matStateCC = getMatStateCC(TElect0, TIon0, density, 
					  VolFrac, matid);
    std::cerr << C4::node() << " Made it after getMatStateCC" << endl;

    std::cerr << C4::node() << " Made it before getMatStateFC" << endl;
    MatStateFC matStateFC = getMatStateFC(matStateCC);
    std::cerr << C4::node() << " Made it after getMatStateFC" << endl;
   
    RadiationStateField radState(spMesh);

    spP13T->initializeRadiationState(testMaterialProps(matStateCC, matStateFC),
				     radState);

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

    std::cerr << C4::node() << " Made it before setBoundary" << endl;

    setBoundary(alpha, beta, bSrc);

    std::cerr << C4::node() << " Made it after setBoundary" << endl;

    if (tdb.verbose)
    {
	if (C4::node() == 0)
	    cout << "alpha: " << endl;
	cout << alpha << endl;
    
	if (C4::node() == 0)
	    cout << "beta: " << endl;
	cout << beta << endl;
    
	if (C4::node() == 0)
	    cout << "bSrc: " << endl;
	cout << bSrc << endl;
    
    }

    if (C4::node() == 0)
	std::cerr << "Made it before diffSolver ctor" << endl;

    SP<MatrixSolver> spMatrixSolver = new MatrixSolver(spMesh, pcg_db);
    spDiffSolver = new DS(diffdb, spMesh, spMatrixSolver);
    
    if (C4::node() == 0)
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
	using rtt_3T_testP13T::min;
	using rtt_3T_testP13T::max;
	
	double (*pabs)(double) = std::abs;
	std::pointer_to_unary_function<double,double> opAbs =
	    std::ptr_fun(pabs);

	double minval;
	double maxval;

	minval = min(TElec, opAbs);
	maxval = max(TElec, opAbs);
	if (C4::node() == 0)
	    cout << "TElec (min), (max): " << minval << " " <<
		maxval << endl;

	minval = min(TIon, opAbs);
	maxval = max(TIon, opAbs);
	if (C4::node() == 0)
	    cout << "TIon (min), (max): " << minval << " " <<
		maxval << endl;
	
	minval = min(radState.phi, opAbs);
	maxval = max(radState.phi, opAbs);
	if (C4::node() == 0)
	    cout << "radState.phi (min), (max): " << minval
		 << " " << maxval << endl;

	minval = min(radState.F, opAbs);
	maxval = max(radState.F, opAbs);
	if (C4::node() == 0)
	    cout << "radState.F (min), (max): " << minval
		 << " " << maxval << endl;
    }

    std::cerr << C4::node() << " Made it before solve3T" << endl;
    
    RadiationStateField newRadState(spMesh);
    ccsf QEEM(spMesh);
    ccsf REEM(spMesh);
	
    spP13T->solve3T(newRadState, QEEM, REEM,
		    electEnergyDep, ionEnergyDep, momentumDeposition,
		    TElec, TIon, *spDiffSolver, dt,
		    testMaterialProps(matStateCC, matStateFC), velocity,
		    radState, QRad, QElectron, QIon, alpha, beta, bSrc);

    std::cerr << C4::node() << " Made it after solve3T" << endl;

    Assert(rtt_3T_testP13T::isContinuous<MT>(newRadState.F, spMesh));

    if (cycle % tdb.dumpcycles == 0)
    {
	static int dumpno = 1;
	std::cerr << "before gmvDump()" << std::endl;
	gmvDump(newRadState, TElec, TIon, dumpno++, cycle, time);
	std::cerr << "after gmvDump()" << std::endl;
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
	
    if (cycle % tdb.dumpcycles == 0)
    {
	int dumpno = cycle / tdb.dumpcycles;
    
	std::cerr << "before creating dump file name" << std::endl;

	std::string fname = rtt_3T_testP13T::getFileName("testFullP13T.",
						  ".dat", dumpno);
	
	std::cerr << "after creating dump file name" << std::endl;

	ccsf TElect(spMesh);
	newMatStateCC.getElectronTemperature(TElect);

	std::cerr << "before dumpInZ()" << std::endl;
	rtt_3T_testP13T::dumpInZ(fname, cycle, time, TElect);
	std::cerr << "after dumpInZ()" << std::endl;
    }

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
    using rtt_3T_testP13T::setBoundary;
    
    setBoundary(alpha, diffdb.alpha_left, diffdb.alpha_right,
		diffdb.alpha_front, diffdb.alpha_back,
		diffdb.alpha_bottom, diffdb.alpha_top);
    setBoundary(beta, diffdb.beta_left, diffdb.beta_right,
		diffdb.beta_front, diffdb.beta_back,
		diffdb.beta_bottom, diffdb.beta_top);

    if (tdb.Te_bottom > 0.0)
    {
	const RadiationPhysics radphys(units);
	double phi_bottom = 0.0;
	radphys.getPlanck(tdb.Te_bottom, phi_bottom);

	phi_bottom *= 4.0*PhysicalConstants::pi;

	setBoundary(bSrc, tdb.src_left, tdb.src_right,
		    tdb.src_front, tdb.src_back,
		    phi_bottom*diffdb.alpha_bottom, tdb.src_top);
    }
    else
    {
	setBoundary(bSrc, tdb.src_left, tdb.src_right,
		    tdb.src_front, tdb.src_back,
		    tdb.src_bottom, tdb.src_top);
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

    using rtt_3T_testP13T::sum;

    const double deltaRadEnergy = sum(deltaRadEnergyCC)*volpcell;
    if (C4::node() == 0)
	cout << "deltaRadEnergy: " << deltaRadEnergy
	     << "\trate: " << deltaRadEnergy/dt << endl;

    const double electEnergyDep = sum(electEnergyDepCC)*volpcell;
    if (C4::node() == 0)
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

    const double ionEnergyDep = sum(ionEnergyDepCC)*volpcell;
    if (C4::node() == 0)
	cout << "  ionEnergyDep: " << ionEnergyDep
	     << "\trate: " << ionEnergyDep/dt << endl;

    const double energyDep = deltaRadEnergy + electEnergyDep + ionEnergyDep;
    if (C4::node() == 0)
	cout << "energyDep: " << energyDep
	     << "\trate: " << energyDep/dt << endl;

    const double inhomosrc = sum(QRad) * volpcell;
    const double qsrc = inhomosrc + (sum(QElectron) +
				     sum(QIon))*volpcell;

    if (C4::node() == 0)
	cout << "Volume src: " << qsrc << endl;

    double bndsrc;
    bndsrc  = (tdb.src_left + tdb.src_right) * ny*dy * nz*dz ;
    bndsrc += (tdb.src_front + tdb.src_back) * nz*dz * nx*dx;
    bndsrc += (tdb.src_bottom + tdb.src_top) * nx*dx * ny*dy;

    if (C4::node() == 0)
	cout << "Boundary src: " << bndsrc << endl;

    const double externsrc = bndsrc + qsrc;
    if (C4::node() == 0)
	cout << "Total external src: " << externsrc << endl;

    const double timedepsrc = sum(radState.phi) * volpcell / (c*dt);
    if (C4::node() == 0)
	cout << "time dep src: " << timedepsrc << endl;
    
    const double timedeprem = sum(newRadState.phi) * volpcell / (c*dt);
    if (C4::node() == 0)
	cout << "time dep removal: " << timedeprem << endl;
    
    ccsf sigAbs(spMesh);
    matStateCC.getSigmaAbsorption(1,sigAbs);

    ccsf temp(spMesh);
    temp = newRadState.phi*sigAbs;
    const double absorption = sum(temp) * volpcell;
    if (C4::node() == 0)
	cout << "absorption: " << absorption << endl;

    temp = newRadState.phi*REEM;
    const double emisrem = sum(temp) * volpcell;
    if (C4::node() == 0)
	cout << "emissive removal: " << emisrem << endl;
    
    const double emission = sum(QEEM) *	volpcell;
    if (C4::node() == 0)
	cout << "emission: " << emission << endl;

    // cout << "emission-absorption: " << emission - absorption << endl;

    const double leakage = calcLeakage(newRadState);
    
    if (C4::node() == 0)
	cout << "leakage (x-y): " << leakage << endl;

    const double totsrc = inhomosrc + bndsrc + timedepsrc + emission;
    if (C4::node() == 0)
	cout << "total radiation sources: " << totsrc << endl;

    const double totrem = leakage + absorption + emisrem + timedeprem;
    if (C4::node() == 0)
	cout << "total radiation removal: " << totrem << endl;
    
    const double balance = totsrc - totrem;
    if (C4::node() == 0)
	cout << "Rad only Balance: " << balance << endl;

    const double relbal = balance / totsrc;
    if (C4::node() == 0)
	cout << "Relative Balance: " << relbal << endl;

    if (C4::node() == 0)
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
	    if (C4::node() == 0 && diffdb.alpha_bottom != 0)
	    {
		double alpha_bottom = diffdb.alpha_bottom;
		double beta_bottom = diffdb.beta_bottom;
		double src_bottom = tdb.src_bottom;

		using rtt_3T_testP13T::BOTTOM;
		double F_bottom = radstate.F(k,l,0,BOTTOM);
		
		// calculate phi along bottom face
		double phi_bottom = (src_bottom - beta_bottom*F_bottom) /
		    alpha_bottom;
		
		leakage += phi_bottom/4.0 + F_bottom/2.0;
		// leakage += alpha_bottom*phi_bottom - beta_bottom*F_bottom;
	    }

	    if (C4::node() == C4::nodes() - 1 && diffdb.alpha_top != 0)
	    {
		double alpha_top = diffdb.alpha_top;
		double beta_top = diffdb.beta_top;
		double src_top = tdb.src_top;

		using rtt_3T_testP13T::TOP;
		double F_top = radstate.F(k,l,nz-1,TOP);
		
		// calculate phi along top face
		double phi_top = (src_top - beta_top*F_top) / alpha_top;

		leakage += phi_top/4.0 + F_top/2.0;
		// leakage += alpha_top*phi_top - beta_top*F_top;
	    }
	}
    }

    leakage *= dx*dy;

    C4::gsum(leakage);
    
    return leakage;
}

//---------------------------------------------------------------------------//
//                              end of testFullP13T.cc
//---------------------------------------------------------------------------//
