//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstIMC_Objects_Builder.cc
 * \author Thomas M. Evans
 * \date   Sat Aug 23 11:17:44 2003
 * \brief  IMC_Objects_Builder unit test.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Soft_Equivalence.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "rng/Random.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/RZWedge_Mesh.hh"
#include "mc/RZWedge_Builder.hh"
#include "mc/Rep_Topology.hh"
#include "mc/Comm_Patterns.hh"
#include "mc/Global_Mesh_Data.hh"
#include "../Release.hh"
#include "../IMC_Objects_Builder.hh"
#include "../Frequency.hh"
#include "../Multigroup_Particle.hh"
#include "imc_test.hh"
#include "IMC_Test.hh"

using namespace std;

using rtt_imc_test::Parser;

using rtt_imc::IMC_Objects_Builder;

using rtt_mc::OS_Mesh;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::RZWedge_Builder;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::Comm_Patterns;
using rtt_mc::Global_Mesh_Data;

using rtt_rng::Rnd_Control;

using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh                            MT;
typedef rtt_mc::RZWedge_Mesh                       RZ;
typedef rtt_imc::Multigroup_Frequency              MG;
typedef rtt_imc::Multigroup_Particle<OS_Mesh>      MGP;
typedef rtt_imc::Multigroup_Particle<RZWedge_Mesh> MGPRZ;
typedef rtt_imc_test::IMC_CDI_Interface<MGPRZ>     CDI_Interface;
typedef rtt_imc_test::IMC_Flat_Interface<MGP>      Flat_Interface;

typedef rtt_dsxx::SP<MG>                                 SP_Frequency;
typedef rtt_dsxx::SP<rtt_imc::Mat_State<MT> >            SP_Mat_State;
typedef rtt_dsxx::SP<rtt_imc::Opacity<MT,MG> >           SP_Opacity;
typedef rtt_dsxx::SP<rtt_imc::Diffusion_Opacity<MT> >    SP_Diff_Opacity;

typedef rtt_dsxx::SP<rtt_imc::Source<MT,MG,MGP> >         SP_Source;
typedef rtt_dsxx::SP<rtt_imc::Tally<MT> >                 SP_Tally;
typedef rtt_dsxx::SP<rtt_imc::Random_Walk<MT> >           SP_Random_Walk;
typedef rtt_dsxx::SP<rtt_imc::Extrinsic_Surface_Tracker>  SP_Tracker;
typedef rtt_dsxx::SP<rtt_imc::Source_Builder<MT,MG,MGP> > SP_Source_Builder;

typedef rtt_dsxx::SP<MG>                                 RZ_Frequency;
typedef rtt_dsxx::SP<rtt_imc::Mat_State<RZ> >            RZ_Mat_State;
typedef rtt_dsxx::SP<rtt_imc::Opacity<RZ,MG> >           RZ_Opacity;
typedef rtt_dsxx::SP<rtt_imc::Diffusion_Opacity<RZ> >    RZ_Diff_Opacity;

typedef rtt_dsxx::SP<rtt_imc::Source<RZ,MG,MGPRZ> >         RZ_Source;
typedef rtt_dsxx::SP<rtt_imc::Tally<RZ> >                   RZ_Tally;
typedef rtt_dsxx::SP<rtt_imc::Random_Walk<RZ> >             RZ_Random_Walk;
typedef rtt_dsxx::SP<rtt_imc::Extrinsic_Surface_Tracker>    RZ_Tracker;
typedef rtt_dsxx::SP<rtt_imc::Source_Builder<RZ,MG,MGPRZ> > RZ_Source_Builder; 
 

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// derived IMC_Objects_Builder class

template<class IT, class M, class P>
class DBuilder : public IMC_Objects_Builder<IT, M, MG, P>
{
  private:
    typedef IMC_Objects_Builder<IT,M,MG,P> Base;

  public:
    // Constructor.
    DBuilder(SP<IT> intface) : IMC_Objects_Builder<IT,M,MG,P>(intface)
    {
	Ensure (Base::interface);
    }

    // Builder.
    void build_IMC_objects(SP<M>                    mesh,
			   SP<Topology>             topology,
			   SP<Global_Mesh_Data<M> > mesh_data,
			   SP<Comm_Patterns>        comm_patterns,
			   SP<Rnd_Control>          rnd_control)
    {
	Require (mesh);
	Require (topology);
	Require (mesh_data);
	Require (comm_patterns);
	Require (rnd_control);

	Require (!Base::frequency);
	Require (!Base::mat_state);
	Require (!Base::opacity);
	Require (!Base::diff_opacity);
	Require (!Base::source);
	Require (!Base::tally);
	Require (!Base::random_walk);
	Require (!Base::tracker);
	Require (!Base::source_builder);

	Require (mesh->num_cells() == topology->num_cells(rtt_c4::node()));

	// first build the material state objects, the appropriate internal
	// function will be called depending upon the base class of the
	// interface (CDI_Data_interface or Flat_Data_Interface)
	Base::build_mat_state_objects(mesh);
	Check (Base::frequency);
	Check (Base::mat_state);
	Check (Base::opacity);

	// build the source
	Base::build_source_builder_object(mesh, topology);
	Base::build_source_object(mesh, rnd_control, comm_patterns);
	Check (Base::source);
	Check (Base::source_builder);

	// build the tally
	Base::build_tally_object(mesh);
	Check (Base::tally);

	// build the particle transport objects
	Base::build_time_dependent_particle_objects(mesh);
	Base::build_time_independent_particle_objects(mesh, mesh_data);
    }

    void reset()
    {
	Base::reset(Base::frequency);
	Base::reset(Base::mat_state);
	Base::reset(Base::opacity);
	Base::reset(Base::diff_opacity);
	Base::reset(Base::source);
	Base::reset(Base::tally);
	Base::reset(Base::random_walk);
	Base::reset(Base::tracker);
	Base::reset(Base::source_builder);
    }
};

//---------------------------------------------------------------------------//

void cdi_objects_test()
{
    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("RZWedge_Input"));
    SP<RZWedge_Builder> mb(new RZWedge_Builder(parser));
    SP<RZ> mesh = mb->build_Mesh();

    if (mesh->num_cells() != 6) ITFAILS;

    // build a Topology:  a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a Global_Mesh_Data
    SP<Global_Mesh_Data<RZ> > mesh_data(
	new Global_Mesh_Data<RZ>(topology, *mesh));

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // make the interface (with hybrid on)
    SP<CDI_Interface> interface(new CDI_Interface(1));

    // make the objects builder
    DBuilder<CDI_Interface, RZ, MGPRZ> builder(interface);

    builder.build_IMC_objects(mesh, topology, mesh_data, patterns, rcon); 

    // get objects

    RZ_Frequency    frequency = builder.get_Frequency();
    RZ_Mat_State    mat       = builder.get_Mat_State();
    RZ_Opacity      opacity   = builder.get_Opacity();
    RZ_Diff_Opacity diff      = builder.get_Diffusion_Opacity();

    // some simple mat objects tests (these are tested thoroughly---with the
    // same interface---in tstCDI_Mat_State_Builder)
    {
	if (!frequency) ITFAILS;
	if (!mat)       ITFAILS;
	if (!opacity)   ITFAILS;
	if (!diff)      ITFAILS;

	if (!frequency->is_multigroup()) ITFAILS;

	if (mat->num_cells() != mesh->num_cells())   ITFAILS;
	if (!soft_equiv(mat->get_rho(3), 1.0))       ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(4), 0.2)) ITFAILS;

	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 0.5 * 1.0)) ITFAILS;

	// from tstCDI_Mat_State_Builder
	double ros2 = 0.600756358394604;
	double ros4 = 1.20140466622632;

	if (!soft_equiv(diff->get_Rosseland_opacity(2), ros2)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(4), ros4)) ITFAILS;
    }

    RZ_Source source                 = builder.get_Source();
    RZ_Source_Builder source_builder = builder.get_Source_Builder(); 

    // some simple source tests (more extensive in the flat case)
    {
	if (!source)         ITFAILS;
	if (!source_builder) ITFAILS;
    }

    RZ_Tally tally = builder.get_Tally();

    // some simple tests on tally
    {
	if (!tally) ITFAILS;

	if (!tally->get_RW_Sub_Tally())      ITFAILS;
	if (!tally->get_Surface_Sub_Tally()) ITFAILS;
    }

    // get particle objects
    RZ_Random_Walk rwalk   = builder.get_Random_Walk();
    RZ_Tracker     tracker = builder.get_Surface_Tracker();
    {
	if (!rwalk)  ITFAILS;
	if (!tracker) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Object builder ok for CDI_Data_Interface.");  
}

//---------------------------------------------------------------------------//

void flat_objects_test()
{
    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));
    rtt_rng::rn_stream = 0;

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build a Topology:  a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a Global_Mesh_Data
    SP<Global_Mesh_Data<OS_Mesh> > mesh_data(
	new Global_Mesh_Data<OS_Mesh>(topology, *mesh));

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // make the interface (with hybrid off)
    SP<Flat_Interface> interface(new Flat_Interface(mb, true));

    // make the objects builder
    DBuilder<Flat_Interface, MT, MGP> builder(interface);

    builder.build_IMC_objects(mesh, topology, mesh_data, patterns, rcon);  

    // get objects

    SP_Frequency    frequency = builder.get_Frequency();
    SP_Mat_State    mat       = builder.get_Mat_State();
    SP_Opacity      opacity   = builder.get_Opacity();
    SP_Diff_Opacity diff      = builder.get_Diffusion_Opacity();

    // some simple mat objects tests (these are tested thoroughly---with the
    // same interface---in tstFlat_Mat_State_Builder)
    {
	if (!frequency) ITFAILS;
	if (!mat)       ITFAILS;
	if (!opacity)   ITFAILS;
	if (diff)       ITFAILS;

	if (!frequency->is_multigroup()) ITFAILS;

	if (mat->num_cells() != mesh->num_cells())   ITFAILS;
	if (!soft_equiv(mat->get_rho(3), 1.0))       ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(4), 0.2)) ITFAILS;

	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 0.1))  ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 2), 0.1))  ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 3), 0.1))  ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 1), 0.02)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 2), 0.02)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 0.02)) ITFAILS;
    }

    SP_Source source                 = builder.get_Source();
    SP_Source_Builder source_builder = builder.get_Source_Builder(); 

    // some simple source tests--> see tstSource_Builder.cc for detailed
    // tests; we simply copied some simple ones here for a sanity check
    {
	if (!source)         ITFAILS;
	if (!source_builder) ITFAILS;

	// get the global numbers for each species
	int global_nsstot  = source_builder->get_nsstot();
	int global_ncentot = source_builder->get_ncentot();
	int global_nvoltot = source_builder->get_nvoltot();

	// check to make sure globals and locals match
	{
	    int local_nsstot  = source->get_nsstot();
	    int local_ncentot = source->get_ncentot();
	    int local_nvoltot = source->get_nvoltot();
	    rtt_c4::global_sum(local_nsstot);
	    rtt_c4::global_sum(local_ncentot);
	    rtt_c4::global_sum(local_nvoltot);

	    if (local_nsstot != global_nsstot)   ITFAILS;
	    if (local_ncentot != global_ncentot) ITFAILS;
	    if (local_nvoltot != global_nvoltot) ITFAILS;
	}
    
	// get deterministically calculated global values of total energy
	// from the source builder; NOTE: this is the initial census
	// energy
	double vol_ew = source_builder->get_evoltot();
	double ss_ew  = source_builder->get_esstot();
	double cen_ew = source_builder->get_initial_census_energy();

	// check to hand calculations of same energies (from xess
	// spreadsheet file source.xs4)
	double hand_vol = 12.747095;
	double hand_ss  = 248.280745;
	double hand_cen = 1630.676880;

	if (!soft_equiv(vol_ew, hand_vol, 1.e-4)) ITFAILS;
	if (!soft_equiv(ss_ew, hand_ss, 1.e-4))   ITFAILS;
	if (!soft_equiv(cen_ew, hand_cen, 1.e-4)) ITFAILS;
    }

    SP_Tally tally = builder.get_Tally();

    // some simple tests on tally
    {
	if (!tally) ITFAILS;

	if (tally->get_RW_Sub_Tally())      ITFAILS;
	if (tally->get_Surface_Sub_Tally()) ITFAILS;
    }

    // get particle objects
    SP_Random_Walk rwalk   = builder.get_Random_Walk();
    SP_Tracker     tracker = builder.get_Surface_Tracker();
    {
	if (rwalk)   ITFAILS;
	if (tracker) ITFAILS;
    }

    // make sure it resets
    builder.reset();
    if (builder.get_Random_Walk())       ITFAILS;
    if (builder.get_Surface_Tracker())   ITFAILS;
    if (builder.get_Tally())             ITFAILS;
    if (builder.get_Source_Builder())    ITFAILS;
    if (builder.get_Source())            ITFAILS;
    if (builder.get_Frequency())         ITFAILS;
    if (builder.get_Mat_State())         ITFAILS;
    if (builder.get_Opacity())           ITFAILS;         
    if (builder.get_Diffusion_Opacity()) ITFAILS;

    // the objects should still be here though
    if (!tally)          ITFAILS;
    if (rwalk)           ITFAILS;
    if (tracker)         ITFAILS;
    if (!source)         ITFAILS;
    if (!source_builder) ITFAILS;
    if (!frequency)      ITFAILS;
    if (!mat)            ITFAILS;
    if (!opacity)        ITFAILS;
    if (diff)            ITFAILS;
	
    if (rtt_imc_test::passed)
	PASSMSG("Object builder ok for Flat_Data_Interface.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	cdi_objects_test();
	flat_objects_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstIMC_Objects_Builder, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstIMC_Objects_Builder Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstIMC_Objects_Builder on " << rtt_c4::node()
	 << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstIMC_Objects_Builder.cc
//---------------------------------------------------------------------------//
