//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstTally.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:11:40 2001
 * \brief  Tally test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Release.hh"
#include "../Tally.hh"
#include "../Tally_Builder.hh"
#include "../Surface_Tracking_Interface.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/Surface_Descriptor.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "mc/Math.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

using rtt_imc_test::Parser;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::Surface_Descriptor;
using rtt_imc::Tally;
using rtt_imc::Tally_Builder;
using rtt_imc::Random_Walk_Sub_Tally;
using rtt_imc::Surface_Sub_Tally;
using rtt_dsxx::SP;
using std::pow;
using rtt_mc::global::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace tally_test
{

class Interface : public rtt_imc::Surface_Tracking_Interface
{
  private:
    vector<Surface_Descriptor> surfaces;
    vector<double>             cosines;

    int hybrid;

  public:
   
    Interface(int h, int s) : surfaces(s), hybrid(h), cosines(5)
    {
	cosines[0] = -1.0;
	cosines[1] = -0.5;
	cosines[2] =  0.0;
	cosines[3] =  0.5;
	cosines[4] =  1.0;
    }

    int get_hybrid_diffusion_method() const { return hybrid; }
    int number_of_surfaces() const { return surfaces.size(); }

    const vector<Surface_Descriptor>& get_surface_data() const
    {
	return surfaces;
    }

    const vector<double>& get_bin_cosines() const
    {
	return cosines;
    }  
};

}

//---------------------------------------------------------------------------//

void Tally_Test()
{
    // make a mesh and tally
    SP<Parser> parser(new Parser("OS_Input"));
    OS_Builder mb(parser);
    SP<OS_Mesh> mesh = mb.build_Mesh();
    Tally<OS_Mesh> t(mesh);

    if (t.num_faces() != 4) ITFAILS;

    // there should be no random walk sub tally
    if (t.get_RW_Sub_Tally()) ITFAILS;

    // now create a random walk sub_tally
    SP<Random_Walk_Sub_Tally> rwsub(new Random_Walk_Sub_Tally);
    t.assign_RW_Sub_Tally(rwsub);
    if (!t.get_RW_Sub_Tally()) ITFAILS;
    
    if (t.num_cells() != 6) ITFAILS;

    std::vector<double> omega(2, 0.0);

    // running and constant cell-number-dependent sum
    int sum_i = 0;
    int const_sum_i = 0;
    for (int i = 1; i <= mesh->num_cells(); i++)
	const_sum_i += i;

    // particles to "run" per cell
    int ppcell = 2;

    // particle count (j-th particle in cell i)
    int pcount = 0;

    SP<Random_Walk_Sub_Tally> st = t.get_RW_Sub_Tally();
	
    // add some stuff and check the tally
    for (int j = 1; j <= 2; j++)
    {
	sum_i = 0;
	for (int i = 1; i <= mesh->num_cells(); i++)
	{
	    // increment particle count (j-th particle in cell i)
	    pcount = (j-1)*mesh->num_cells() + i;

	    // increment running cell-number-dependent sum
            sum_i += i;

	    t.deposit_energy(i, i * 10.0);
	    t.accumulate_ewpl(i, i * 20.0);
	    t.accumulate_cen_info(i, i * 30);
	    t.accum_n_effscat();
	    t.accum_n_thomscat();
	    t.accum_n_killed();
	    t.accum_ew_killed(i * 2.0);
	    t.accum_n_escaped();
	    t.accum_ew_escaped(i * 1.0, 1 + (i % 4));  // Hit all four faces
	    t.accum_n_bndcross();
	    t.accum_n_reflections();
	    t.get_RW_Sub_Tally()->accum_n_random_walks();
	    st->accum_sphere_radii(1.1);
	    st->accum_step_length(0.9);

	    // momentum checks - hardwired for 2 dimensions, 2 depositions
	    omega[0] = (1.0 - pow(-1.0,j))/2.0;
	    omega[1] = (1.0 - pow(-1.0,j+1))/2.0;
	    t.accumulate_momentum(i,static_cast<double>(i),omega);

	    // checks
	    if (t.get_energy_dep(i) != j * i * 10)     ITFAILS;
	    if (t.get_accum_ewpl(i) != j * i * 20)     ITFAILS;
	    if (t.get_new_ecen(i) != j * i * 30)       ITFAILS;
	    if (t.get_new_ncen(i) != j * 1)            ITFAILS;

	    if (t.get_accum_ew_killed() != ((j-1)*const_sum_i+sum_i)*2.0) 
		ITFAILS;

	    if (t.get_accum_n_effscat() != pcount)         ITFAILS;
	    if (t.get_accum_n_thomscat() != pcount)        ITFAILS;
	    if (t.get_accum_n_killed() != pcount)          ITFAILS;
	    if (t.get_accum_n_escaped() != pcount)         ITFAILS;
	    if (t.get_accum_n_bndcross() != pcount)        ITFAILS;
	    if (t.get_accum_n_reflections() != pcount)     ITFAILS;
	    if (st->get_accum_n_random_walks() != pcount)  ITFAILS;

	}
    }

    // now test accumulator with non-default (default=1) values
    t.accum_n_effscat(5);
    t.accum_n_thomscat(20);
    t.accum_n_killed(10);
    t.accum_n_escaped(2);
    t.accum_n_bndcross(3);
    t.accum_n_reflections(6);
    st->accum_n_random_walks(11);

    if (t.get_accum_n_effscat() != pcount + 5)        ITFAILS;
    if (t.get_accum_n_thomscat() != pcount + 20)      ITFAILS;
    if (t.get_accum_n_killed() != pcount + 10)        ITFAILS;
    if (t.get_accum_n_escaped() != pcount + 2)        ITFAILS;
    if (t.get_accum_n_bndcross() != pcount + 3)       ITFAILS;
    if (t.get_accum_n_reflections() != pcount + 6)    ITFAILS;

    if (t.get_RW_Sub_Tally()->get_accum_n_random_walks() != pcount + 11)
	ITFAILS;

    if (!soft_equiv(t.get_RW_Sub_Tally()->get_accum_step_lengths(),
		    pcount * 0.9)) 
	ITFAILS;

    if (t.get_RW_Sub_Tally()->get_accum_n_spheres() != pcount) ITFAILS;
    

    if (!soft_equiv(t.get_RW_Sub_Tally()->get_accum_sphere_radii(),
		    pcount * 1.1)) 
	ITFAILS;

    for (int i = 1; i <= mesh->num_cells(); i++)
	t.accumulate_cen_info(i, i*1.0, i);

    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	if (t.get_new_ecen(i) != ppcell*i*30 + i) ITFAILS;
	if (t.get_new_ncen(i) != ppcell + i)      ITFAILS;
    }
    
    // momentum deposition checks
    std::vector<double> mom(2, 0.0);
    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	mom = t.get_momentum_dep(i);
	if (mom.size() != 2)                            ITFAILS;
	if (!soft_equiv(mom[0],static_cast<double>(i))) ITFAILS;
	if (!soft_equiv(mom[1],static_cast<double>(i))) ITFAILS;
    }
    
    // totals test
    if (t.get_energy_dep_tot() != 10.0 * sum_i*ppcell)       ITFAILS;
    if (t.get_new_ecen_tot() != (30.0*sum_i*ppcell) + sum_i) ITFAILS;
    if (t.get_new_ncen_tot() != pcount + sum_i)              ITFAILS;
    if (t.get_ew_escaped() != sum_i*ppcell)                  ITFAILS;
    if (t.num_cells() != 6)                                  ITFAILS;

    double total_ew_escaped = 0.0;
    for (int face = 1; face <=4; ++face) total_ew_escaped += t.get_ew_escaped(face);

    if (!soft_equiv(total_ew_escaped, t.get_ew_escaped())) ITFAILS;
    
    if (rtt_imc_test::passed)
	PASSMSG("Tally tests ok.")
}

//---------------------------------------------------------------------------//

void Tally_Builder_Test()
{
    const int groups = 5;

    // make a mesh and tally
    SP<Parser> parser(new Parser("OS_Input"));
    OS_Builder mb(parser);
    SP<OS_Mesh> mesh = mb.build_Mesh();

    // make an interface that has no random walk or surfaces
    SP<tally_test::Interface> null(new tally_test::Interface(0, 0));

    Tally_Builder<OS_Mesh> nullb(null, groups);
    SP<Tally<OS_Mesh> > nullt = nullb.build_Tally(mesh);

    if (!nullt)                  ITFAILS;
    if (nullt->num_cells() != 6) ITFAILS;

    // all tallies should be null
    for (int i = 1; i <= 6; i++)
    {
	if (nullt->get_energy_dep(i) != 0) ITFAILS;
	if (nullt->get_accum_ewpl(i) != 0) ITFAILS;
	if (nullt->get_new_ecen(i) != 0)   ITFAILS;
	if (nullt->get_new_ncen(i) != 0)   ITFAILS;
    }

    if (nullt->get_energy_dep_tot() != 0)      ITFAILS;
    if (nullt->get_new_ecen_tot() != 0)        ITFAILS;
    if (nullt->get_new_ncen_tot() != 0)        ITFAILS;
    if (nullt->get_accum_n_effscat() != 0)     ITFAILS;
    if (nullt->get_accum_n_thomscat() != 0)    ITFAILS;
    if (nullt->get_accum_n_killed() != 0)      ITFAILS;
    if (nullt->get_accum_ew_killed() != 0)     ITFAILS;
    if (nullt->get_ew_escaped() != 0)          ITFAILS;
    if (nullt->get_accum_n_escaped() != 0)     ITFAILS;
    if (nullt->get_accum_n_bndcross() != 0)    ITFAILS;
    if (nullt->get_accum_n_reflections() != 0) ITFAILS;

    if (nullt->get_RW_Sub_Tally())      ITFAILS;
    if (nullt->get_Surface_Sub_Tally()) ITFAILS;

    if (nullt->get_num_surface_tallies() != 0) ITFAILS;

    // make an interface that has random walk and surfaces
    SP<tally_test::Interface> has(new tally_test::Interface(1, 1));

    Tally_Builder<OS_Mesh> hasb(has, groups);
    SP<Tally<OS_Mesh> > hast = hasb.build_Tally(mesh);

    if (!hast->get_RW_Sub_Tally())            ITFAILS;
    if (hast->get_num_surface_tallies() != 5) ITFAILS;

    SP<Random_Walk_Sub_Tally> rw_tally = hast->get_RW_Sub_Tally();

    if (rw_tally->get_accum_n_random_walks() != 0) ITFAILS;
    if (rw_tally->get_accum_n_spheres() != 0)      ITFAILS;
    if (rw_tally->get_accum_sphere_radii() != 0.0) ITFAILS;
    if (rw_tally->get_accum_step_lengths() != 0.0) ITFAILS;

    for (int group = 1; group <= groups; ++group)
    {
	SP<Surface_Sub_Tally> s_tally = hast->get_Surface_Sub_Tally(group);

        if (!s_tally) ITFAILS;
	
	if (s_tally->get_number_surfaces() != 1)              ITFAILS;
	if (s_tally->get_outward_weight_tally(1).size() != 4) ITFAILS;
	if (s_tally->get_mesh_size() != 4 )                   ITFAILS;
    }
	
    if (rtt_imc_test::passed)
	PASSMSG("Tally_Builder tests ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is essentially a scalar test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_imc::release() 
		 << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// Tally tests
	Tally_Test();
	Tally_Builder_Test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstTally, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstTally Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstTally on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstTally.cc
//---------------------------------------------------------------------------//
